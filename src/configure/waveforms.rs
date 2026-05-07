use config::{Config, ConfigError, File};
use float_eq::float_eq;
use serde::Deserialize;

#[allow(non_snake_case)]
#[derive(Debug, Deserialize, Clone)]
pub struct Configure {
    pub buffer_size_mus  : f64,
    pub sampling_time_mus: f64,
    pub sipm_size_mm2    : f64,
    pub noise_std        : f64,
    pub adc_per_pe       : f64,
    pub pe_sigma         : f64,
    pub baseline         : u16,
}

impl Configure {
    pub fn new(filename: &str) -> Result<Self, ConfigError> {
        let s = Config::builder()
            .add_source(File::with_name(filename))
            .build()?;

        // You can deserialize (and thus freeze) the entire configuration as
        let out : Self = s.try_deserialize()?;

        let off = out.buffer_size_mus % out.sampling_time_mus;
        if ! (float_eq!(off,                   0.0, abs<=1e-6) ||
              float_eq!(off, out.sampling_time_mus, abs<=1e-6)) {
            panic!("Buffer size must be a multiple of the sampling time.")
        }
        Ok(out)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use pretty_assertions::assert_eq;
    use float_eq::assert_float_eq;

    #[test]
    fn test_new() {
        let c = Configure::new("conf/electronics.toml").unwrap();
        assert_float_eq!(c.buffer_size_mus  ,   504.0, ulps<=2);
        assert_float_eq!(c.sampling_time_mus,   16e-3, ulps<=2);
        assert_float_eq!(c.sipm_size_mm2    , 6.0*6.0, ulps<=2);
        assert_float_eq!(c.noise_std        ,     0.1, ulps<=2);
        assert_float_eq!(c.adc_per_pe       ,    20.0, ulps<=2);
        assert_float_eq!(c.pe_sigma         ,     0.2, ulps<=2);
        assert_eq!      (c.baseline         ,      30);
    }

    #[test]
    fn test_new_err() {
        assert!(Configure::new("does_not_exist").is_err());
    }
}
