

pub struct SensorDB {
    pub sensor_id: Vec<u16>,
    pub gain     : Vec<f64>,
}

impl SensorDB {
    pub fn pcolina() -> SensorDB {
        let sensor_id = (0..5).flat_map(|b| (0..25).map(move |s| 1000*b + s)).collect();
        let gain      = vec![1f64; 125];
        Self{sensor_id, gain}
    }

    pub fn get_index(&self, id: u16) -> Option<usize> {
        self.sensor_id
            .iter()
            .position(|&sid| sid==id)
    }
}
