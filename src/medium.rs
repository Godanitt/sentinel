use crate::io::read_csv;

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

    pub fn pulse_shape(medium: Medium) -> (Vec<f64>, Vec<f64>) {
        match medium {
            Medium::Xenon => {
                let mut data = read_csv("data/pulse_shape_vuv4_lxe.csv", " ", 0).expect("Could not read file");
                (data.remove(0), data.remove(0))
            },
            Medium::Argon => {
                let mut data = read_csv("data/pulse_shape_vuv4_lar.csv", " ", 0).expect("Could not read file");
                (data.remove(0), data.remove(0))
            },
            Medium::Krypton => unimplemented!(),
        }
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
