use std::num::NonZero;

use nalgebra::{Point2, Point3};
use kiddo::{ImmutableKdTree, SquaredEuclidean};

use crate::io::read_csv;
use crate::field_point::FieldPoint;

#[derive(Clone, Debug)]
pub struct Field{
    pub tree        : ImmutableKdTree<f64, 2>,
    pub field_points: Vec<FieldPoint>,
}

impl Field {
    pub fn new(field_points: Vec<FieldPoint>) -> Self {
        let points : Vec<[f64; 2]> =
            field_points.iter()
                        .map(|fp| [fp.r(), fp.z()])
                        .collect();

        let tree = ImmutableKdTree::new_from_slice(&points);
        Field{tree, field_points}
    }

    #[allow(non_snake_case)]
    pub fn from_file(filename: &str, to_mm: f64, to_Vpercm: f64, offset_z: f64, invert_z: bool) -> Self {
        let data = read_csv(filename, " ", 8).expect("Could not read field file");

        let mut previous_pos  = Point2::<f64>::origin();
        let mut  current_line = 999_999_999_999_999usize;

        let sign = if invert_z {-1.0} else {1.0};
        let field_points : Vec<FieldPoint> =
            data.into_iter()
                .rev()
                .map(|row| vec![ row[0] * to_mm
                               , row[1] * to_mm * sign + offset_z
                               , row[2]
                               , row[3] * to_Vpercm
                ])
                .filter_map(|row| FieldPoint::from_csv_row(row, &mut previous_pos, &mut current_line))
                .collect();

        Field::new(field_points)
    }

    /// Vector map in mm and V/cm, prepared from a checked COMSOL header by Python.
    /// Keep invalid samples so nearest-neighbour lookup cannot bridge NaN regions.
    pub fn from_vector_file(filename: &str) -> std::io::Result<Self> {
        use std::io::{BufRead, BufReader, Error, ErrorKind};
        use nalgebra::Vector2;
        let mut points = Vec::new();
        for (index, line) in BufReader::new(std::fs::File::open(filename)?).lines().enumerate() {
            let line = line?;
            let values: Result<Vec<f64>, _> = line.split_whitespace().map(str::parse).collect();
            let row = values.map_err(|_| Error::new(ErrorKind::InvalidData, format!("Invalid field row {}", index+1)))?;
            if row.len() != 4 || !row[0].is_finite() || !row[1].is_finite() || row[0] < 0. {
                return Err(Error::new(ErrorKind::InvalidData, "Expected finite r>=0, z, Er, Ez"));
            }
            let direction = Vector2::new(-row[2], -row[3]); // electrons drift opposite E
            let mag = direction.norm();
            points.push(FieldPoint::new(Point2::new(row[0], row[1]), mag, 0, direction));
        }
        if points.is_empty() || !points.iter().any(|p| p.mag.is_finite() && p.mag > 0.) {
            return Err(Error::new(ErrorKind::InvalidData, "No valid nonzero field vectors"));
        }
        Ok(Self::new(points))
    }

    pub fn find_nearest(&self, pos: &Point3<f64>) -> &FieldPoint {
        let r   = (pos.x*pos.x + pos.y*pos.y).sqrt();
        let z   = pos.z;
        let one = NonZero::new(1).unwrap();
        let index = self.tree.nearest_n::<SquaredEuclidean>(&[r, z], one)[0].item;
        self.field_points.get(index as usize).unwrap()
    }
}



#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::Vector2;
    use pretty_assertions::assert_eq;
    use float_eq::assert_float_eq;
    use crate::assert_vector2_eq;

    #[test]
    fn test_from_file() {
        use std::io::Write;
        use tempfile::NamedTempFile;

        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, "1\n2\n3\n4\n5\n6\n7\n8").unwrap(); // 8 rows to skip
        writeln!(file, "1.2 -24.5 6 7.89").unwrap();       // so it's easy to test
        writeln!(file, "1.2 -14.5 6 7.89").unwrap();       // repeat the same values
        writeln!(file, "1.2  -4.5 6 7.89").unwrap();       // the actual rows

        let field = Field::from_file(file.path().to_str().unwrap(), 1., 1., 0., true);

        assert_eq!(field.field_points.len(), 2); // first one is skipped for being the end of the line

        let dir = Vector2::new(0.0, -1.0);
        for (i, fp) in field.field_points.iter().enumerate() {
            let z = 14.5 + (10*i) as f64;
            assert_float_eq!  (fp.r(),  1.2, ulps <= 2);
            assert_float_eq!  (fp.z(),    z, ulps <= 2);
            assert_float_eq!  (fp.mag, 7.89, ulps <= 2);
            assert_eq!        (fp.line, 6usize);
            assert_vector2_eq!(fp.dir,  dir, ulps <= 2);
        }
    }

    #[test]
    fn test_find_nearest() {
        let points = vec![
            FieldPoint::new(Point2::<f64>::new(0., 0.), 1., 0, Vector2::zeros()), // 0
            FieldPoint::new(Point2::<f64>::new(1., 1.), 1., 0, Vector2::zeros()), // 1
            FieldPoint::new(Point2::<f64>::new(2., 2.), 1., 0, Vector2::zeros()), // 2
            FieldPoint::new(Point2::<f64>::new(3., 3.), 1., 0, Vector2::zeros()), // 3
        ];

        let field = Field::new(points.clone());
        let p0 = Point3::new(0.1, 0.0, 0.1); // closest is 0
        let p1 = Point3::new(0.9, 0.0, 0.9); // closest is 1
        let p2 = Point3::new(1.9, 0.0, 2.1); // closest is 2
        let p3 = Point3::new(3.1, 0.0, 2.9); // closest is 3

        let c0 = field.find_nearest(&p0);
        let c1 = field.find_nearest(&p1);
        let c2 = field.find_nearest(&p2);
        let c3 = field.find_nearest(&p3);

        assert_float_eq!(c0.r(), points[0].r(), ulps <= 2);
        assert_float_eq!(c1.r(), points[1].r(), ulps <= 2);
        assert_float_eq!(c2.r(), points[2].r(), ulps <= 2);
        assert_float_eq!(c3.r(), points[3].r(), ulps <= 2);
    }


    #[test]
    fn vector_field_preserves_z_and_electron_sign() {
        use std::io::Write;
        let mut file = tempfile::NamedTempFile::new().unwrap();
        writeln!(file, "1 2 3 -4").unwrap();
        writeln!(file, "2 2 NaN NaN").unwrap();
        let field = Field::from_vector_file(file.path().to_str().unwrap()).unwrap();
        let point = field.find_nearest(&Point3::new(1., 0., 2.));
        assert_float_eq!(point.mag, 5., abs <= 1e-12);
        assert_float_eq!(point.dir.x, -0.6, abs <= 1e-12);
        assert_float_eq!(point.dir.y, 0.8, abs <= 1e-12);
        assert_float_eq!(point.z(), 2., abs <= 1e-12);
        assert!(field.find_nearest(&Point3::new(2., 0., 2.)).mag.is_nan());
    }
}
