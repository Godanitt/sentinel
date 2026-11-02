//! Batched transport of explicitly identified electrons prepared by Python.
use std::{collections::{HashMap, HashSet}, error::Error, io, path::PathBuf, rc::Rc, time::{Duration, Instant}};
use clap::Parser;
use hdf5_metno as hdf5;
use nalgebra::Point3;
use rayon::prelude::*;
use sentinel::{
    configure::trajectories::Configure,
    field::Field,
    geometry::{Cone, Geometry},
    io::{ScalarHdf5Writer, Writer, hdf5_types::{Endpoint, TrajectoryPoint}},
    tracker::Tracker,
};

#[derive(Parser)]
struct Args {
    #[arg(long)] input: PathBuf,
    #[arg(long)] field: String,
    #[arg(long)] output: PathBuf,
    #[arg(long, value_parser = ["vector", "legacy"])] mode: String,
    #[arg(long)] endpoints: bool,
    #[arg(long)] trajectories: bool,
}

#[derive(hdf5::H5Type, Clone)]
#[repr(C)]
struct Initial {
    event: u32,
    cell_id: u64,
    x0: f64,
    y0: f64,
    z0: f64,
}
#[derive(hdf5::H5Type, Clone)]
#[repr(C)]
struct EventCell {
    event_id: u32,
    cell_id: u64,
}

/// Numerical map support, not detector walls or material surfaces.
struct MapSupport { rmin: f64, rmax: f64, zmin: f64, zmax: f64 }
impl Geometry for MapSupport {
    fn is_within(&self, p: &Point3<f64>) -> bool {
        p.coords.iter().all(|v| v.is_finite()) &&
        p.x.hypot(p.y) >= self.rmin && p.x.hypot(p.y) <= self.rmax && p.z >= self.zmin && p.z <= self.zmax
    }
    fn cathode_z(&self) -> f64 { self.zmax }
}
enum TransportVolume { Legacy(Cone), Map(MapSupport) }
impl Geometry for TransportVolume {
    fn is_within(&self, p: &Point3<f64>) -> bool {
        match self { Self::Legacy(g) => g.is_within(p), Self::Map(g) => g.is_within(p) }
    }
    fn cathode_z(&self) -> f64 {
        match self { Self::Legacy(g) => g.cathode_z(), Self::Map(g) => g.cathode_z() }
    }
}
fn invalid(message: impl Into<String>) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidInput, message.into())
}

fn progress(phase: &str, done: usize, total: usize, details: &str) {
    eprintln!("SENTINEL_PROGRESS\t{phase}\t{done}\t{total}\t{details}");
}

fn main() -> Result<(), Box<dyn Error>> {
    let args = Args::parse();
    let (field, geometry, timestep) = if args.mode == "vector" {
        let field = Field::from_vector_file(&args.field)?;
        let rmin = field.field_points.iter().map(|p| p.r()).fold(f64::INFINITY, f64::min);
        let rmax = field.field_points.iter().map(|p| p.r()).fold(0., f64::max);
        let zmin = field.field_points.iter().map(|p| p.z()).fold(f64::INFINITY, f64::min);
        let zmax = field.field_points.iter().map(|p| p.z()).fold(f64::NEG_INFINITY, f64::max);
        (field, TransportVolume::Map(MapSupport { rmin, rmax, zmin, zmax }), 0.01)
    } else {
        let conf = Configure::new("conf/test.toml")?;
        let field = Field::from_file(&args.field, conf.field_to_mm, conf.field_to_Vpercm,
                                     conf.neck_length, conf.field_invert_z);
        let geometry = Cone::new(conf.rmin, conf.form_factor, conf.neck_length, conf.drift_length);
        (field, TransportVolume::Legacy(geometry), conf.t_step)
    };
    let tracker = Tracker::new(field, geometry, timestep);
    let inputs = hdf5::File::open(&args.input)?;
    let dataset = inputs.dataset("initial_positions")?;
    let total = dataset.size();
    let total_cells = inputs.dataset("cells")?.size();
    if total == 0 { return Err(invalid("No initial electrons").into()); }

    // If any sampled electron in a cell lacks field support, discard that
    // entire cell. Keep validating the full input before starting transport.
    let mut invalid_cells = HashSet::new();
    let mut electrons_per_cell = HashMap::<u64, usize>::new();
    let mut last_validation_report = Instant::now();
    for start in (0..total).step_by(8192) {
        let batch = dataset.read_slice_1d::<Initial, _>(start..(start+8192).min(total))?;
        for (offset, row) in batch.iter().enumerate() {
            if row.event as usize != start + offset {
                return Err(invalid("Initial electron IDs must be unique, contiguous and ordered").into());
            }
            *electrons_per_cell.entry(row.cell_id).or_default() += 1;
            let p = Point3::new(row.x0, row.y0, row.z0);
            if !tracker.valid_start(&p) {
                invalid_cells.insert(row.cell_id);
            }
        }
        let scanned = (start + batch.len()).min(total);
        if last_validation_report.elapsed() >= Duration::from_millis(500) || scanned == total {
            progress("validation", scanned, total,
                     &format!("electrones {scanned}/{total} | celdas descartadas {}", invalid_cells.len()));
            last_validation_report = Instant::now();
        }
    }
    println!();
    let cancelled_electrons: usize = invalid_cells.iter()
        .map(|cell| electrons_per_cell.get(cell).copied().unwrap_or(0)).sum();
    let expected_electrons = total - cancelled_electrons;
    if !invalid_cells.is_empty() {
        println!("Se descartan {} celdas ({} electrones) porque al menos una posición inicial no tiene campo válido.",
                 invalid_cells.len(), cancelled_electrons);
    }
    let output = Rc::new(hdf5::File::create_excl(&args.output)?);
    let processed = ScalarHdf5Writer::<EventCell>::new(Rc::clone(&output), "processed", 1024)?;
    let endpoints = args.endpoints.then(|| ScalarHdf5Writer::<Endpoint>::new(Rc::clone(&output), "endpoints", 1024)).transpose()?;
    let trajectories = args.trajectories.then(|| ScalarHdf5Writer::<TrajectoryPoint>::new(Rc::clone(&output), "trajectories", 4096)).transpose()?;

    // At most one complete trajectory per worker; endpoints require constant memory.
    let workers = rayon::current_num_threads().max(1);
    let batch_size = if args.trajectories { workers } else { 1024 };
    let mut processed_electrons = 0usize;
    let mut last_progress_report = Instant::now();
    for start in (0..total).step_by(batch_size) {
        let stop = (start + batch_size).min(total);
        let batch = dataset.read_slice_1d::<Initial, _>(start..stop)?.to_vec()
            .into_iter().filter(|row| !invalid_cells.contains(&row.cell_id)).collect::<Vec<_>>();
        let batch_len = batch.len();
        let results: Result<Vec<_>, io::Error> = batch.into_par_iter().map(|row| {
            let mut points = Vec::new();
            let pos = Point3::new(row.x0, row.y0, row.z0);
            let (last, steps) = tracker.propagate_checked(pos, 2_000_000, |p, i| {
                if args.trajectories {
                    points.push(TrajectoryPoint { event: row.event, x: p.x as f32, y: p.y as f32,
                                                z: p.z as f32, t: (i as f64 * timestep) as f32 });
                }
            })?;
            let endpoint = Endpoint { event: row.event, x0: row.x0 as f32, y0: row.y0 as f32, z0: row.z0 as f32,
                x1: last.x as f32, y1: last.y as f32, z1: last.z as f32, t: (steps as f64 * timestep) as f32 };
            Ok((EventCell { event_id: row.event, cell_id: row.cell_id }, endpoint, points))
        }).collect();
        for (id, endpoint, points) in results? {
            processed.write(id)?;
            if let Some(writer) = &endpoints { writer.write(endpoint)?; }
            if let Some(writer) = &trajectories { writer.write_many(points)?; }
        }
        processed_electrons += batch_len;
        let completed_cells = if stop == total {
            total_cells
        } else {
            dataset.read_slice_1d::<Initial, _>(stop..stop + 1)?[0].cell_id as usize
        };
        if last_progress_report.elapsed() >= Duration::from_millis(500) || stop == total {
            let cancelled_completed = invalid_cells.iter().filter(|&&cell| (cell as usize) < completed_cells).count();
            progress("transport", stop, total,
                     &format!("celdas {completed_cells}/{total_cells} | restan {} | descartadas {cancelled_completed} | electrones {processed_electrons}/{expected_electrons}",
                              total_cells - completed_cells));
            last_progress_report = Instant::now();
        }
    }
    println!();
    processed.flush()?;
    if let Some(writer) = endpoints { writer.flush()?; }
    if let Some(writer) = trajectories { writer.flush()?; }
    output.flush()?;
    println!("Rust: processed {total} electrons.");
    Ok(())
}
