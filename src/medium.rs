use crate::io::read_csv;
use crate::utils::interpolate_1d;

#[derive(clap::ValueEnum, Clone, Debug)]
pub enum Medium {
    Xenon,
    Argon,
    Krypton,
}

impl Medium {
    pub fn light_yield(&self) -> f64 {
        30.0
    }

    pub fn w_i(&self) -> f64 {
        match self {
            Medium::Xenon   => 15.6,
            Medium::Argon   => unimplemented!(),
            Medium::Krypton => unimplemented!(),
        }
    }
    pub fn w_s(&self) -> f64 {
        match self {
            Medium::Xenon   => 17.9,
            Medium::Argon   => unimplemented!(),
            Medium::Krypton => unimplemented!(),
        }
    }

    //                   cumulative fractions, constants
    pub fn time_constants(&self) -> (Vec<f64>, Vec<f64>) {
        match self {
            Medium::Xenon   => (vec![0.03, 1.00], vec![2.0, 42.5]),
            Medium::Argon   => unimplemented!(),
            Medium::Krypton => unimplemented!(),
        }
    }

    pub fn pulse_shape(&self, sampling_time_mus: f64) -> Vec<f32> {
        let (time, amplitude) = match self {
            Medium::Xenon => {
                let mut data = read_csv("data/pulse_shape_vuv4_lxe.csv", " ", 0).expect("Could not read file");
                (data.remove(0), data.remove(0))
            },
            Medium::Argon => {
                let mut data = read_csv("data/pulse_shape_vuv4_lar.csv", " ", 0).expect("Could not read file");
                (data.remove(0), data.remove(0))
            },
            Medium::Krypton => unimplemented!(),
        };
        let time = time.into_iter().map(|t| t*1e-3).collect();
        let out = (0..200).map(|i| (i as f64) * sampling_time_mus)
                          .filter_map(|t| interpolate_1d(&time, &amplitude, t).ok())
                          .collect::<Vec<_>>();

        let norm = out.iter().sum::<f64>();
        out.into_iter()
           .map(|x| (x / norm) as f32)
           .collect()
    }

    pub fn dark_noise_prob(&self, sampling_time_mus: f64, sipm_size_mm2: f64) -> f64 {
        let freq = match self { // Hz
            Medium::Xenon   => 0.137   , // https://www.sciencedirect.com/science/article/pii/S0168900224002730
            Medium::Argon   => 0.978e-5, // extrapolated from above
            Medium::Krypton => 0.134e-3, // extrapolated from above
        };
        freq * 1e-6 * sipm_size_mm2 * sampling_time_mus // per sampling time
    }
}

impl<T> From<T> for Medium where T: ToString {
    fn from(value: T) -> Self {
        let value = value.to_string().to_lowercase();
        match &value[..] {
            "xenon"   => Medium::Xenon,
            "argon"   => Medium::Argon,
            "krypton" => Medium::Krypton,
            _         => panic!("Invalid medium {}", value)
        }
    }
}


impl Into<String> for Medium  {
    fn into(self) -> String {
        match self {
            Medium::Xenon   => "Xenon",
            Medium::Argon   => "Argon",
            Medium::Krypton => "Krypton",
        }.into()
    }
}
