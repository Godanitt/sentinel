use std::io:: Result;
use std::path::Path;
use std::rc::Rc;
use std::time::Instant;

use clap::Parser;
use rayon::prelude::*;
use hdf5_metno as hdf5;

use sentinel::configure::trajectories::Configure;
use sentinel::geometry::Cone;
use sentinel::field::Field;
use sentinel::tracker::Tracker;
use sentinel::io::{Writer, CsvWriter, Hdf5Writer, DatasetShape};
use sentinel::io::hdf5_types::Endpoint;
use sentinel::invalid_input;
use sentinel::progress::MaybeProgressBar;

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct CLI {

    #[arg(short, long)]
    conf: String,

    #[arg(short, long)]
    n_events: usize,

    #[arg(short, long, default_value_t=1000)]
    batch_size: usize,

    #[arg(short, long)]
    output: String,

    #[arg(long, default_value_t=false)]
    overwrite: bool,

    #[arg(short, long, default_value_t=10)]
    threads: usize,

    #[arg(long, default_value_t=false)]
    csv: bool,

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

    rayon::ThreadPoolBuilder::new().num_threads(args.threads).build_global().unwrap();

    let geometry   = Cone::new(conf.rmin, conf.form_factor, conf.neck_length, conf.drift_length);
    let field      = Field::from_file(&conf.field_file, conf.field_to_mm, conf.field_to_Vpercm, conf.neck_length, true);
    let tracker    = Tracker::new(field, geometry.clone(), conf.t_step);
    let sampler    = conf.generator.sampler(&geometry);

    let writer : Box<dyn Writer<Endpoint>> =
        if args.csv {
            let header = "event x0 y0 z0 x1 y1 z1 t".split(" ")
                                                    .map(|v| v.to_string())
                                                    .collect::<Vec<String>>();

            Box::new(CsvWriter::new(&path.to_str().unwrap(), " ", header).unwrap())
        }
        else {
            let file = hdf5::File::create(&path.to_str().unwrap()).unwrap();
            Box::new(Hdf5Writer::<Endpoint>::new(Rc::new(file), "endpoints", 1024, DatasetShape::Scalar).unwrap())
        };

    println!("Initialization time: {:?}", timer.elapsed().as_secs_f64());

    let nbatch = args.n_events.div_ceil(args.batch_size);
    let pb = MaybeProgressBar::new(nbatch, args.batch);

    let timer = Instant::now();
    for batch in 0..nbatch {
        let data : Vec<Endpoint> =
        (0..args.batch_size).into_par_iter()
                            .filter_map( |evt| {
                                let evt = evt + batch * args.batch_size;

                                sampler.sample(&geometry).map(|starting_pos| {
                                    let trajectory = tracker.propagate_from(starting_pos.clone());
                                    let last = trajectory.last().unwrap();
                                    Endpoint{ event: evt as u32
                                            , x0   : starting_pos.x as f32
                                            , y0   : starting_pos.y as f32
                                            , z0   : starting_pos.z as f32
                                            , x1   :         last.x as f32
                                            , y1   :         last.y as f32
                                            , z1   :         last.z as f32
                                            , t    : trajectory.len() as f32 * conf.t_step as f32
                                            }
                                })
                            })
                            .collect();
        writer.write_many(data)?;
        pb.inc(1);
    }
    pb.finish();

    let exe_time = timer.elapsed().as_secs_f64();
    println!( "Execution time for {} events: {:.2} s => {:.1} evt/s or {:.8} s/evt"
            , args.n_events
            , exe_time
            , args.n_events as f64 / exe_time
            , exe_time / args.n_events as f64
            );
    Ok(())
}
