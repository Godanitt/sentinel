use nalgebra::{Point3, Vector3, Rotation3};

use crate::io::read_csv;
use crate::utils::interpolate_1d;
use crate::geometry::Geometry;
use crate::field::Field;
use crate::field_point::FieldPoint;
use crate::random;

#[derive(Debug)]
pub struct Tracker<G: Geometry> {
    field     : Field,
    geometry  : G,
    drift     : (Vec<f64>, Vec<f64>),
    tstep     : f64
}

impl<G: Geometry> Tracker<G> {
    pub fn new(field: Field, geometry: G, tstep: f64) -> Self {
        let data        = read_csv("data/drift_velocity_gushchin.dat", ";", 1).unwrap();
        let mut drift_e = data.iter().map(|x| *x.first().unwrap()).collect::<Vec<f64>>(); // drift field
        let mut drift_v = data.iter().map(|x| *x. last().unwrap()).collect::<Vec<f64>>(); // drift velocity

        // Extra points to avoid going under/over interpolation range
        drift_e.insert(0, 0f64);
        drift_v.insert(0, 0f64);
        drift_e.push(1e6);
        drift_v.push(*drift_v.last().unwrap());

        Self{ field
            , geometry
            , drift : (drift_e, drift_v)
            , tstep
            }
    }

    pub fn propagate_from(&self, pos: Point3<f64>) -> Vec<Point3<f64>> {
        let mut trajectory : Vec<Point3<f64>> = Vec::new();
        trajectory.push(pos);

        while self.geometry.is_within(trajectory.last().unwrap()) {
            let current_pos = trajectory.last().unwrap();
            let direction   = self.compute_step(&current_pos);
            let next_pos = current_pos + direction;
            trajectory.push(next_pos);
        }
        trajectory
    }

    /// Same transport as propagate_from, without storing the trajectory.
    pub fn propagate_endpoint(&self, mut pos: Point3<f64>) -> (Point3<f64>, usize) {
        let mut steps = 0;
        while self.geometry.is_within(&pos) {
            pos += self.compute_step(&pos);
            steps += 1;
        }
        (pos, steps)
    }

    pub fn valid_start(&self, pos: &Point3<f64>) -> bool {
        if !pos.coords.iter().all(|v| v.is_finite()) || !self.geometry.is_within(pos) {
            return false;
        }
        let field = self.field.find_nearest(pos);
        field.mag.is_finite() && field.mag > 0. &&
        field.dir.iter().all(|v| v.is_finite()) && field.dir.norm() > 0.
    }

    /// Bounded integration for the Python interface. Domain/invalid-field exits
    /// are numerical endpoints, without any detector/material classification.
    pub fn propagate_checked<F>(&self, mut pos: Point3<f64>, max_steps: usize, mut visit: F)
        -> std::io::Result<(Point3<f64>, usize)>
    where F: FnMut(&Point3<f64>, usize) {
        if !self.valid_start(&pos) {
            return Err(std::io::Error::new(std::io::ErrorKind::InvalidInput, "Invalid starting position"));
        }
        visit(&pos, 0);
        let mut steps = 0;
        while self.valid_start(&pos) {
            if steps >= max_steps {
                return Err(std::io::Error::other("Transport step limit reached; check field and numerical domain"));
            }
            pos += self.compute_step(&pos);
            if !pos.coords.iter().all(|v| v.is_finite()) {
                return Err(std::io::Error::other("Non-finite transport step"));
            }
            steps += 1;
            visit(&pos, steps);
        }
        Ok((pos, steps))
    }

    fn compute_step(&self, pos: &Point3<f64>) -> Vector3<f64> {
        let closest = self.field.find_nearest(pos);

        let dv          = self.drift_velocity(closest.mag);
        let step_length = dv * self.tstep;
        let field_dir   = closest.to_vec3(&pos);                           // direction set by e field
        let smeared_dir = self.randomize_direction(&closest, step_length); // diffusion smeared

        let rot = Rotation3::face_towards(&field_dir, &Vector3::y());
        let step = rot * smeared_dir;

        step
    }

    fn randomize_direction(&self, fp: &FieldPoint, step_length: f64) -> Vector3<f64> {
        let dt = self.  transverse_diffusion(fp.mag);
        let dl = self.longitudinal_diffusion(fp.mag);
        let sigmat = (2f64 * dt * self.tstep * 1e-4).sqrt(); // 1e-4 = cm2/s to mm2/us
        let sigmal = (2f64 * dl * self.tstep * 1e-4).sqrt();

        let smearx = random::normal(0., sigmat);
        let smeary = random::normal(0., sigmat);
        let smearz = random::normal(0., sigmal);

        Vector3::new(smearx, smeary, step_length + smearz)
    }

    fn drift_velocity(&self, edrift: f64) -> f64 {
        interpolate_1d(&self.drift.0, &self.drift.1 , edrift).unwrap()
    }

    fn transverse_diffusion(&self, edrift: f64) -> f64 {
        // from NEST, parametrization of red data points in
        // https://arxiv.org/pdf/1609.04467
        37.368 * edrift.powf(0.093452) * (-8.1651e-5 * edrift).exp()
    }

    fn longitudinal_diffusion(&self, edrift: f64) -> f64 {
        // from NEST, parametrization of magenta data points in
        // https://arxiv.org/pdf/1911.11580
        66.350 * edrift.powf(-0.24855) + 36.85 * (-edrift / 35.661).exp()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use float_eq::assert_float_eq;
    use nalgebra::{Point2, Vector2};
    use crate::geometry::Cone;

    fn homogeneous_field() -> Field {
        let points   = vec![
            FieldPoint::new(Point2::<f64>::new(0., 11.), 1e3, 0, Vector2::new(0., -1.)), // startpoint
            FieldPoint::new(Point2::<f64>::new(0.,  1.), 1e3, 0, Vector2::new(0., -1.)), //
            FieldPoint::new(Point2::<f64>::new(0.,  0.), 1e3, 0, Vector2::new(0., -1.)), //   endpoint
        ];
        Field::new(points)
    }
    fn default_geometry() -> Cone          { Cone::new(1., 1., 1.0, 10.) }
    fn default_tracker () -> Tracker<Cone> { Tracker ::new(homogeneous_field(), default_geometry(), 1.) }

    #[test]
    fn test_tracker_new() {
        let tracker = default_tracker();
        assert!(tracker.drift.0.len() > 1);
        assert!(tracker.drift.1.len() > 1);
    }

    #[test]
    fn test_drift_velocity() {
        let tracker = default_tracker();
        assert_float_eq!(tracker.drift_velocity( 100.), 1.4938051013386355, rel <= 1e-7);
        assert_float_eq!(tracker.drift_velocity(1000.), 2.2215991423542030, rel <= 1e-7);
        assert_float_eq!(tracker.drift_velocity(2000.), 2.5070676356146624, rel <= 1e-7);
    }

    #[test]
    fn test_propagate_thisone() {
        let tracker  = Tracker::new(homogeneous_field(), default_geometry(), 1e-2);
        let t = tracker.propagate_from(Point3::new(0., 0., 9.));
        assert!(t.len() > 100, "Track length too short: {}", t.len());
        assert!(t.last().unwrap().z < 2e-2); // close to exit
    }

    #[test]
    fn test_propagate_slow() {
        let tracker  = Tracker::new(homogeneous_field(), default_geometry(), 5e-4);
        let t = tracker.propagate_from(Point3::new(0., 0., 9.));

        assert!(t.len() > 1000, "Track length too short: {}", t.len());
        assert!(t.last().unwrap().z < 1e-3); // close to exit
    }

    #[test]
    fn checked_transport_rejects_invalid_start_and_enforces_limit() {
        let tracker = default_tracker();
        assert!(tracker.propagate_checked(Point3::new(100., 0., 9.), 10, |_, _| {}).is_err());
        assert!(tracker.propagate_checked(Point3::new(0., 0., 9.), 0, |_, _| {}).is_err());
    }

    #[test]
    fn endpoint_only_transport_leaves_geometry() {
        let tracker = default_tracker();
        let (end, steps) = tracker.propagate_endpoint(Point3::new(0., 0., 9.));
        assert!(steps > 0);
        assert!(!tracker.geometry.is_within(&end));
    }
}
