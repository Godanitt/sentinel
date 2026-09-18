use ndarray::Array2;

use crate::io::hdf5_types::SensorHit;
use crate::random::{normal, poisson};
use crate::sensor_db::SensorDB;


#[derive(Clone, Debug)]
pub struct Waveforms<T> {
    pub event: u64,
    pub data : Array2<T>,
}


pub fn digitize_sensor_hits(sensor_hits: Vec<SensorHit>, nsens: usize, ntime: usize, sampling_time: f64, pretrigger: f64, sensor_db: &SensorDB) -> Waveforms<u32> {
    if sensor_hits.len() == 0 {
        panic!("There should be at least one sensor hit here")
    }

    let pretrigger = (pretrigger / sampling_time).floor() as usize;
    let     event = sensor_hits.first().unwrap().event;
    let mut data  = Array2::<u32>::zeros((nsens, ntime));

    for sh in sensor_hits.into_iter() {
        let id    = sensor_db.get_index(sh.sensor_id).expect("Unknown sensor_id");
        let index = (sh.t as f64 / sampling_time).floor() as usize + pretrigger;

        // micro optimization to avoid double lookup
        // data[[id, index]] = data[[id, index]].saturating_add(1);
        let item = data.get_mut((id, index)).unwrap();
        *item = item.saturating_add(1);
    }

    Waveforms{event, data}
}


pub fn add_dark_noise(wfs: Waveforms<u32>, p: f64) -> Waveforms<u32> {
    let Waveforms{event, mut data} = wfs;

    // It's unrealistic to saturate, but better safe than sorry
    data.iter_mut()
        .for_each(|v| *v = v.saturating_add(poisson(p) as u32));
    Waveforms{event, data}
}


pub fn apply_pulse_shape(wfs: Waveforms<u32>, pulse_shape: &Vec<f32>) -> Waveforms<f32> {
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


pub fn to_adc(wfs: Waveforms<f32>, adc_per_pe: f64, pe_sigma: f64) -> Waveforms<f32> {
    let Waveforms{event, mut data} = wfs;

    data.iter_mut()
        .for_each(|a| *a *= (normal(1.0, pe_sigma) * adc_per_pe) as f32);
    Waveforms{event, data}

}

pub fn add_white_noise(wfs: Waveforms<f32>, std: f64) -> Waveforms<f32> {
    let Waveforms{event, mut data} = wfs;

    data.iter_mut()
        .for_each(|v| *v += normal(0.0, std) as f32);
    Waveforms{event, data}
}


pub fn digitize_voltage(wfs: Waveforms<f32>) -> Waveforms<u16> {
    let Waveforms{event, data} = wfs;
    let daq_max    = u16::MAX as f32;
    let mut daq_wf = Array2::<u16>::zeros(data.raw_dim());

    daq_wf.iter_mut()
          .zip(data.iter())
          .for_each(|(out, in_)| {
              *out = in_.floor().min(daq_max) as u16;
          });

    Waveforms{event, data: daq_wf}
}


pub fn add_pedestal(wfs: Waveforms<u16>, ped: u16) -> Waveforms<u16> {
    let Waveforms{event, mut data} = wfs;

    data.iter_mut()
        .for_each(|v| *v = ped.saturating_sub(*v));
    Waveforms{event, data}
}
