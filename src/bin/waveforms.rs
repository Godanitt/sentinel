use std::io:: Result;
use std::path::Path;
use std::rc::Rc;
use std::time::Instant;

use itertools::Itertools;
use clap::Parser;
use rayon::prelude::*;
use hdf5_metno as hdf5;

use sentinel::invalid_input;
use sentinel::configure::waveforms::Configure;
use sentinel::electronics::{add_dark_noise, add_pedestal, add_white_noise, apply_pulse_shape, digitize_sensor_hits, digitize_voltage, to_adc, Waveforms};
use sentinel::io::{Writer, ArrayHdf5Writer, ScalarHdf5Writer, read_hdf5};
use sentinel::io::hdf5_types::SensorHit;
use sentinel::medium::Medium;
use sentinel::progress::MaybeProgressBar;
use sentinel::sensor_db::SensorDB;

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct CLI {

    #[arg(short, long)]
    conf: String,

    #[arg(short, long)]
    input: String,

    #[arg(short, long)]
    output: String,

    #[arg(long, value_enum, default_value_t=Medium::Xenon)]
    medium: Medium,

    #[arg(long, default_value_t=false)]
    overwrite: bool,

    #[arg(short, long, default_value_t=10)]
    threads: usize,
}

pub fn main() -> Result<()> {
    let timer = Instant::now();

    let args = CLI::parse();
    let conf = Configure::new(&args.conf).unwrap();
    let path = Path::new(&args.output);
    if path.exists() & !args.overwrite {
        return invalid_input!("Outfile already exists!");
    }

    let ofile  = hdf5::File::create(&path.to_str().unwrap()).unwrap();
    let ofile  = Rc::new(ofile);
    let nsens  = 125_usize;
    let ntime  = (conf.buffer_size_mus / conf.sampling_time_mus).ceil() as usize;
    let shape  = vec![nsens, ntime];

    let evt_writer = ScalarHdf5Writer::<u64>::new(Rc::clone(&ofile), "/Waveforms/event_id"  , 1024       ).unwrap();
    let chn_writer = ScalarHdf5Writer::<u16>::new(Rc::clone(&ofile), "/Waveforms/channel_id", 1024       ).unwrap();
    let  wf_writer =  ArrayHdf5Writer::<u16>::new(Rc::clone(&ofile), "/Waveforms/adcs"      ,  128, shape).unwrap();

    let sensor_db       = SensorDB::pcolina();
    let pulse_shape     = args.medium.pulse_shape(conf.sampling_time_mus);
    let dark_noise_prob = args.medium.dark_noise_prob(conf.sampling_time_mus, conf.sipm_size_mm2);
    println!("Initialization time: {:?}", timer.elapsed().as_secs_f64());


    let shits   = read_hdf5::<SensorHit>(&args.input, "/MC/sensor_hits")?;
    let events  = shits.into_iter()
                       .chunk_by(|sh| sh.event)
                       .into_iter()
                       .map(|(_, shs)| shs.collect::<Vec<SensorHit>>())
                       .collect::<Vec<Vec<SensorHit>>>();
    let n_evt   = events.len();

    let pb = MaybeProgressBar::new(n_evt, false);

    rayon::ThreadPoolBuilder::new().num_threads(args.threads).build_global().unwrap();

    let timer = Instant::now();

    let mut waveforms : Vec<Waveforms<u16>> =
        events.into_iter()
              .inspect(|_| pb.inc(1))
              .par_bridge()
              .map(|shs| digitize_sensor_hits(shs, nsens, ntime, conf.sampling_time_mus, conf.pretrigger_mus, &sensor_db) )
              .map(|wfs| add_dark_noise(wfs, dark_noise_prob))
              .map(|wfs| apply_pulse_shape(wfs, &pulse_shape) )
              .map(|wfs| to_adc(wfs, conf.adc_per_pe, conf.pe_sigma))
              .map(|wfs| add_white_noise(wfs, conf.noise_std))
              .map(|wfs| digitize_voltage(wfs))
              .map(|wfs| add_pedestal(wfs, conf.baseline))
              .collect();

    pb.finish();

    waveforms.sort_by(|a, b| b.event.cmp(&a.event));
    for wf in waveforms {
        evt_writer.write(wf.event)?;
         wf_writer.write(wf.data)?;
    }
    for channel in sensor_db.sensor_id {
        chn_writer.write(channel)?;
    }


    let exe_time = timer.elapsed().as_secs_f64();
    println!( "Execution time for {} events: {:.2} s => {:.1} evt/s or {:.8} s/evt"
            , n_evt
            , exe_time
            , n_evt as f64 / exe_time
            , exe_time / n_evt as f64
            );
    Ok(())
}
