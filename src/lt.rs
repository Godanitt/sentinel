use std::collections::HashMap;

use kiddo::{ImmutableKdTree, SquaredEuclidean};

use nalgebra::{point, Point2, Rotation2};
use hdf5_metno as hdf5;

use crate::io::read_hdf5;
use crate::io::hdf5_types::{S1LightTablePoint, S2LightTablePoint, DriftMapPoint};
use crate::random::{exp_survival_n, choice_nonexhaustive};

#[derive(Debug)]
pub struct LTprob {
    pub sensor_id: u16,
    pub prob     : f64,
}

#[allow(non_camel_case_types)]
pub struct S1LightTable {
    probs: Vec<Vec<LTprob>>,
    tree : ImmutableKdTree<f32, 3>,
}

#[allow(non_camel_case_types)]
pub struct S2LightTable {
    probs: Vec<Vec<LTprob>>,
    tree : ImmutableKdTree<f32, 2>,
}

pub struct DriftTable {
    points: Vec<Vec<Point2<f32>>>,
    t_aves: Vec<f64>,
    probs : Vec<Vec<f64>>,
    tree  : ImmutableKdTree<f32, 2>,
}


impl S1LightTable {
    pub fn from_file(filename: &str) -> hdf5::Result<Self> {
        let data = read_hdf5::<S1LightTablePoint>(filename, "/S1/LT")?;

        let mut map : HashMap<(u32, u32, u32), Vec<LTprob>> = HashMap::new();
        for p in data {
            map.entry( (p.z.to_bits(), p.x.to_bits(), p.y.to_bits() ) )
                .or_insert(Vec::with_capacity(125))
                .push(LTprob {
                    sensor_id: p.sensor_id,
                    prob     : p.prob as f64,
                });
        }
        let mut points = Vec::with_capacity(map.len());
        let mut probs  = Vec::with_capacity(map.len());
        for ((z, x, y), v) in map {
            points.push([f32::from_bits(z), f32::from_bits(x), f32::from_bits(y)]);
            probs .push(v);
        }
        let tree = ImmutableKdTree::new_from_slice(&points);
        Ok(Self{probs, tree})
    }

    pub fn get(&self, x: f32, y: f32, z: f32) -> &Vec<LTprob> {
        let index = self.tree.nearest_one::<SquaredEuclidean>(&[z, x, y]).item;
        self.probs.get(index as usize).unwrap()
    }
}

impl S2LightTable {
    pub fn from_file(filename: &str) -> hdf5::Result<Self> {
        let data = read_hdf5::<S2LightTablePoint>(filename, "/S2/LT")?;

        let mut map : HashMap<(u32, u32), Vec<LTprob>> = HashMap::new();
        for p in data {
            map.entry( (p.x.to_bits(), p.y.to_bits() ) )
                .or_insert(Vec::with_capacity(125))
                .push(LTprob {
                    sensor_id: p.sensor_id,
                    prob     : p.prob as f64,
                });
        }
        let mut points = Vec::with_capacity(map.len());
        let mut probs  = Vec::with_capacity(map.len());
        for ((x, y), v) in map {
            points.push([f32::from_bits(x), f32::from_bits(y)]);
            probs .push(v);
        }
        let tree = ImmutableKdTree::new_from_slice(&points);
        Ok(Self{probs, tree})
    }

    pub fn get(&self, x: f32, y: f32) -> &Vec<LTprob> {
        let index = self.tree.nearest_one::<SquaredEuclidean>(&[x, y]).item;
        self.probs.get(index as usize).unwrap()
    }
}

impl DriftTable {
    pub fn from_file(filename: &str) -> hdf5::Result<Self> {
        let data = read_hdf5::<DriftMapPoint>(filename, "/Drift/Map")?;

        let mut map : HashMap<(u32, u32, u32), Vec<(Point2<f32>, f64)>> = HashMap::new();
        for p in data {
            map.entry( (p.z0.to_bits(), p.r0.to_bits(), p.t_ave.to_bits()) )
               .or_insert(Vec::with_capacity(100))
               .push( ( point![p.x1, p.y1]
                      , (p.n_dst as f64) / (p.n_src as f64)
                      )
                    );
        }

        let mut sources = Vec::with_capacity(map.len());
        let mut points  = Vec::with_capacity(map.len());
        let mut t_aves  = Vec::with_capacity(map.len());
        let mut probs   = Vec::with_capacity(map.len());

        for ((r, z, t_ave), v) in map {
            sources.push([ f32::from_bits(r), f32::from_bits(z) ]);
            t_aves.push(f32::from_bits(t_ave) as f64);
            let (pts, pbs) = v.into_iter().unzip();
            points.push(pts);
            probs .push(pbs);
        }

        let tree = ImmutableKdTree::new_from_slice(&sources);
        Ok(Self{points, t_aves, probs, tree})
    }

    pub fn get(&self, x: f32, y: f32, z: f32, n: usize, tau: f64) -> (Vec<Point2<f32>>, f64) {
        let rot = Rotation2::new(y.atan2(x));

        let r      = (x*x + y*y).sqrt();
        let index  = self.tree.nearest_one::<SquaredEuclidean>(&[z, r]).item as usize;
        let points = self.points.get(index).unwrap();
        let probs  = self.probs .get(index).unwrap();
        let t_ave  = self.t_aves.get(index).unwrap();

        let n = exp_survival_n(*t_ave, tau, n);

        let points = (0..n).filter_map(|_| choice_nonexhaustive(&points, &probs))
                           .map(|xy| rot * xy)
                           .collect();
        (points, *t_ave)
    }

}
