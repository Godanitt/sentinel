use ndarray::Array2;

use crate::io::hdf5_types::SensorHit;
use crate::random::poisson;
use crate::sensor_db::SensorDB;


#[derive(Clone, Debug)]
pub struct Waveforms<T> {
    pub event: u64,
    pub data : Array2<T>,
}


pub fn digitize_sensor_hits(sensor_hits: Vec<SensorHit>, nsens: usize, ntime: usize, sampling_time: f64, sensor_db: &SensorDB) -> Waveforms<u16> {
    if sensor_hits.len() == 0 {
        panic!("There should be at least one sensor hit here")
    }

    let     event = sensor_hits.first().unwrap().event;
    let mut data  = Array2::<u16>::zeros((nsens, ntime));

    for sh in sensor_hits.into_iter() {
        let id    = sensor_db.get_index(sh.sensor_id).expect("Unknown sensor_id");
        let index = (sh.t as f64 / sampling_time).floor() as usize;

        // micro optimization to avoid double lookup
        // data[[id, index]] = data[[id, index]].saturating_add(1);
        let item = data.get_mut((id, index)).unwrap();
        *item = item.saturating_add(1);
    }

    Waveforms{event, data}
}

pub fn apply_pulse_shape(wfs: Waveforms<u16>, pulse_shape: &Vec<f32>) -> Waveforms<f32> {
    let Waveforms{event, data} = wfs;
    let mut wfs = Array2::<f32>::zeros(data.raw_dim());

    let imax = data.shape()[0];
    let jmax = data.shape()[1];

    for i in 0..imax {
        for j in 0..jmax {
            let a = data[[i,j]] as f32;
            pulse_shape.iter()
                       .enumerate()
                       .filter(|(k, _)| j+k < jmax)
                       .for_each(|(k, s)| {
                           wfs[[i, j+k]] += s * a;
                       })
        }
    }
    Waveforms{event, data:wfs}
}

pub fn add_noise(wfs: Waveforms<u16>, p: f64) -> Waveforms<u16> {
    let Waveforms{event, mut data} = wfs;

    data.iter_mut()
        .for_each(|v| *v += v.saturating_add(poisson(p) as u16));
    Waveforms{event, data}
}
