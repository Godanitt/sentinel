use std::io:: Result;
use std::path::Path;
use std::rc::Rc;
use std::time::Instant;

use clap::Parser;
use itertools::Itertools;
use hdf5_metno as hdf5;
use nalgebra::point;
use rayon::prelude::*;

use sentinel::configure::detsim::Configure;
use sentinel::field::Field;
use sentinel::geometry::Cone;
use sentinel::invalid_input;
use sentinel::io::{Writer, ScalarHdf5Writer, read_hdf5};
use sentinel::io::hdf5_types::{IonizationHit, SensorHit};
use sentinel::lt::{S1LightTable, S2LightTable};
use sentinel::medium::Medium;
use sentinel::progress::MaybeProgressBar;
use sentinel::random::{poisson, multiexpo, uniform, exp_survival};
use sentinel::tracker::Tracker;

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

    #[arg(short, long)]
    maps: String,

    #[arg(long, default_value_t=false)]
    overwrite: bool,

    #[arg(short, long, default_value_t=10)]
    threads: usize,

    #[arg(long, default_value_t=false)]
    batch: bool,
}

pub fn main() -> Result<()> {
    let timer = Instant::now();

    let args = CLI::parse();
    let conf = Configure::new(&args.conf).unwrap();
    let path = Path::new(&args.output);
    if path.exists() & !args.overwrite {
        return invalid_input!("Outfile already exists!");
    }

    let s1_lt = S1LightTable::from_file(&args.maps)?;
    let s2_lt = S2LightTable::from_file(&args.maps)?;

    let geometry = Cone::new(conf.rmin, conf.form_factor, conf.neck_length, conf.drift_length);
    let field    = Field::from_file(&conf.field_file, conf.field_to_mm, conf.field_to_Vpercm, conf.neck_length, true);
    let tracker  = Tracker::new(field, geometry.clone(), conf.t_step);

    let (fractions, constants) = args.medium.time_constants();

    let ofile  = hdf5::File::create(&path.to_str().unwrap()).unwrap();
    let writer = ScalarHdf5Writer::<SensorHit>::new(Rc::new(ofile), "/MC/sensor_hits", 1024).unwrap();

    println!("Initialization time: {:?}", timer.elapsed().as_secs_f64());

    let events = read_hdf5::<IonizationHit>(&args.input, "/MC/ionization_hits")?
        .into_iter()
        .chunk_by(|sh| sh.event)
        .into_iter()
        .map(|(_, shs)| shs.collect::<Vec<IonizationHit>>())
        .collect::<Vec<Vec<IonizationHit>>>();

    let n_evt = events.len();

    let pb = MaybeProgressBar::new(n_evt, args.batch);

    rayon::ThreadPoolBuilder::new().num_threads(args.threads).build_global().unwrap();

    let     timer = Instant::now();
    let mut n_ok  = 0;

    for ihits in events.into_iter() {
        let mut sensor_hits : Vec<SensorHit> =
            ihits.par_iter()
                 .flat_map_iter(|ihit| {
                     let event = ihit.event;
                     let n_ave = ihit.e as f64 * 1e6 / args.medium.w_s();// hits come in MeV, w_s is in eV
                     let n_ph  = poisson(n_ave) as f64;
                     s1_lt.get(ihit.z, ihit.x, ihit.y)
                          .iter() // over sensors
                          .map(move |p| (p.sensor_id, poisson(p.prob*n_ph)))
                          .filter(|(_, n)| *n>0)
                          .flat_map(move |(sid, n)| (0..n).map(move |_| (event, sid)))
                 })
                 .map(|(event, sensor_id)| {
                     let t = multiexpo(&fractions, &constants);
                     SensorHit{event, sensor_id, t: t as f32}
                 })
                 .collect();

        let mut s2_sensor_hits : Vec<SensorHit> =
            ihits.into_par_iter()
                 .flat_map_iter(|ihit| {
                     let event  = ihit.event;
                     let n_ave  = ihit.e as f64 * 1e6 / args.medium.w_i(); // hits come in MeV, w_i is in eV
                     let n_ie   = poisson(n_ave);
                     let origin = point![ihit.x as f64, ihit.y as f64, ihit.z as f64];
                     (0..n_ie).map(move |_| origin.clone())
                              .map(|origin| tracker.propagate_from(origin))
                              .filter_map(move |trajectory| {
                                  let tdrift   = trajectory.len() as f64 * conf.t_step;
                                  let last     = trajectory.last().unwrap();
                                  let endpoint = point![last.x as f32, last.y as f32, last.z as f32];
                                  if exp_survival(tdrift, conf.lifetime) {
                                      Some((event, tdrift, endpoint))
                                  }
                                  else {
                                      None
                                  }
                              })
                 })
                 .flat_map_iter(|(evt, ie_time, ep)| {
                     s2_lt.get(ep.x, ep.y)
                          .iter() // over sensors
                          .map(|p| (p.sensor_id, poisson(p.prob*args.medium.light_yield())))
                          .filter(|(_, n)| *n>0)
                          .flat_map(move |(sid, n)| (0..n).map(move |_| (evt, ie_time, sid)))
                 })
                 .map(|(event, ie_time, sensor_id)| {
                     // HACK: uniform emission in time within el emission range
                     let t = ie_time + uniform(0.0, 0.1) + multiexpo(&fractions, &constants);
                     SensorHit{event, sensor_id, t: t as f32}
                 })
                 .collect();

        sensor_hits.append(&mut s2_sensor_hits);
        sensor_hits.sort_by(|a, b| b.event.cmp(&a.event));
        if !sensor_hits.is_empty() { n_ok += 1; }

        writer.write_many(sensor_hits)?;
        pb.inc(1);
    }
    pb.finish();

    let exe_time = timer.elapsed().as_secs_f64();
    println!( "{}/{} = {:.1}% events produced signal", n_ok, n_evt, n_ok as f32 / n_evt as f32 * 100.0);
    println!( "Execution time for {} events: {:.2} s => {:.1} event/s or {:.8} s/event"
            , n_evt
            , exe_time
            , n_evt as f64 / exe_time
            , exe_time / n_evt as f64
            );
    Ok(())
}
