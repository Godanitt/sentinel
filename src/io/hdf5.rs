use std::cell::RefCell;
use std::io::Result;
use std::rc::Rc;

use hdf5_metno::{self as hdf5, Error, Extent};
use hdf5_metno::filters::BloscShuffle;

use ndarray::s;

use crate::io::Writer;

pub fn read_hdf5<T: hdf5::H5Type + Clone>(filename: &str, dataset : &str) -> hdf5::Result<Vec<T>> {
    let file    = hdf5::File::open(filename)?;
    let dataset = file.dataset(dataset)?;
    dataset.read_slice_1d::<T,_>(s![..]).map(|v| v.to_vec())
}

pub enum DatasetShape {
    Scalar,
    Vector(usize),
    Matrix(usize, usize),
}

pub struct Hdf5Writer<T: hdf5::H5Type>{
    #[allow(dead_code)]
    _file     : Rc<hdf5::File>, // file needs to live while dataset lives, this is a way of ensuring that
    dataset   : hdf5::Dataset,
    chunk_size: usize,
    cache     : RefCell<Vec<T>>,
    shape     : DatasetShape,
}


impl<T: hdf5::H5Type> Hdf5Writer<T>{
    pub fn new(file: Rc<hdf5::File>, dataset: &str, chunk_size: usize, shape: DatasetShape) -> hdf5::Result<Self> {
        if chunk_size == 0 {
            return Err(Error::Internal("Hdf5Writer::new: invalid chunk size 0".to_owned()));
        }
        let dataset = match shape {
            DatasetShape::Scalar      => file.new_dataset::<T>()
                                             .chunk((chunk_size,))
                                             .blosc_zlib(4, BloscShuffle::Byte)
                                             .shape(Extent::resizable(0))
                                             .create(dataset)?,
            DatasetShape::Vector(n)   => file.new_dataset::<T>()
                                             .chunk((chunk_size, n))
                                             .blosc_zlib(4, BloscShuffle::Byte)
                                             .shape((Extent::resizable(0), n))
                                             .create(dataset)?,
            DatasetShape::Matrix(r,c) => file.new_dataset::<T>()
                                             .chunk((chunk_size, r, c))
                                             .blosc_zlib(4, BloscShuffle::Byte)
                                             .shape((Extent::resizable(0), r, c))
                                             .create(dataset)?,

        };

        let cache = Vec::with_capacity(chunk_size);
        Ok(Hdf5Writer{_file: file, dataset, chunk_size, cache: RefCell::new(cache), shape})
    }


    fn dump_cache(&self) -> Result<()> {
        let n_write = self.cache.borrow().len();
        if n_write > 0 {
            let size_before = self.dataset.shape()[0];
            match self.shape {
                DatasetShape::Scalar      => {
                    self.dataset.resize( size_before + n_write       )?;
                    self.dataset.write_slice(&self.cache.borrow()[..], size_before..)?;
                },
                DatasetShape::Vector(n)   => {
                    self.dataset.resize((size_before + n_write, n   ))?;
                    self.dataset.write_slice(&self.cache.borrow()[..], (size_before.., ..))?;
                },
                DatasetShape::Matrix(r,c) => {
                    self.dataset.resize((size_before + n_write, r, c))?;
                    self.dataset.write_slice(&self.cache.borrow()[..], (size_before.., .., ..))?;
                },
            }
            self.cache.borrow_mut().clear();
        }
        Ok(())
    }

    pub fn write(&self, item: T) -> Result<()> {
        self.cache.borrow_mut().push(item);

        if self.cache.borrow().len() == self.chunk_size {
            self.dump_cache()
        }
        else {
            Ok(())
        }
    }

    pub fn flush(&self) -> Result<()> {
        self.dump_cache()
    }
}

impl<T: hdf5::H5Type> Writer<T> for Hdf5Writer<T> {
    fn write(&self, value: T) -> Result<()> { self.write(value) }
}

impl<T: hdf5::H5Type> Drop for Hdf5Writer<T> {
    fn drop(&mut self) {
        self.flush().unwrap()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use pretty_assertions::assert_eq;
    use float_eq::assert_float_eq;
    use hdf5;
    use crate::utils::tempfile;
    use crate::io::hdf5_types::{Endpoint, IonizationHit};

    fn extremes() -> (Endpoint, Endpoint) {
        let e0 = Endpoint{ event:  0, x0:  1.0, y0:  2.0, z0:  3.0, x1:  4.0, y1:  5.0, z1:  6.0, t :  7.0};
        let e1 = Endpoint{ event: 10, x0: 11.0, y0: 12.0, z0: 13.0, x1: 14.0, y1: 15.0, z1: 16.0, t : 17.0};
        (e0, e1)
    }

    #[test]
    fn test_read_hdf5() {
        let filename = "data/kr_hits.h5";
        let data_read = read_hdf5::<IonizationHit>(&filename, "/MC/ionization_hits").unwrap();
        assert_eq!(data_read.len(), 150);

        let first = data_read.first().unwrap();
        let last  = data_read.last ().unwrap();

        assert_eq!      (first.event   , 0);
        assert_eq!      (first.track_id, 1);
        assert_float_eq!(first.x       , -  2.1049984 , abs<=1e-6);
        assert_float_eq!(first.y       , -  0.6166766 , abs<=1e-6);
        assert_float_eq!(first.z       ,  115.36533   , abs<=1e-6);
        assert_float_eq!(first.e       ,    0.00100068, abs<=1e-6);

        assert_eq!      (last.event   , 9);
        assert_eq!      (last.track_id, 1);
        assert_float_eq!(last.x       ,  26.793226  , abs<=1e-6);
        assert_float_eq!(last.y       , -29.236149  , abs<=1e-6);
        assert_float_eq!(last.z       ,  83.36896   , abs<=1e-6);
        assert_float_eq!(last.e       ,   0.00191707, abs<=1e-6);
    }

    #[test]
    fn dataset_writer_new() {
        let (_dir, filename) = tempfile("dataset_writer_new");
        let file             = hdf5::File::create(filename).unwrap();
        let writer           = Hdf5Writer::<Endpoint>::new(Rc::new(file), "a_dataset", 123);
        assert!(writer.is_ok());
    }

    #[test]
    fn dataset_writer_new_invalid_filename() {
        let (_dir, filename) = tempfile("dataset_writer_new");
        let file             = hdf5::File::create(filename).unwrap();
        let writer           = Hdf5Writer::<Endpoint>::new(Rc::new(file), "", 123);
        assert!(writer.is_err());
    }

    #[test]
    fn dataset_writer_new_invalid_chunksize() {
        let (_dir, filename) = tempfile("dataset_writer_new");
        let file             = hdf5::File::create(filename).unwrap();
        let writer           = Hdf5Writer::<Endpoint>::new(Rc::new(file), "a_dataset", 0);
        assert!(writer.is_err());
    }

    #[test]
    fn dataset_writer_round_trip_single() {
        let (_dir, filename) = tempfile("round_trip_single");
        let file             = hdf5::File::create(filename.clone()).unwrap();
        let writer           = Hdf5Writer::<Endpoint>::new(Rc::new(file), "a_dataset", 1).unwrap();

        let data_written = extremes().0;
        writer.write(data_written).unwrap();

        let data_read = read_hdf5::<Endpoint>(&filename, "a_dataset").unwrap();

        assert_eq!(data_read.len(), 1);

        let extremes = data_read.first().unwrap();
        assert_eq!      (extremes.event, 0);
        assert_float_eq!(extremes.x0, 1.0, ulps<=2);
        assert_float_eq!(extremes.y0, 2.0, ulps<=2);
        assert_float_eq!(extremes.z0, 3.0, ulps<=2);
        assert_float_eq!(extremes.x1, 4.0, ulps<=2);
        assert_float_eq!(extremes.y1, 5.0, ulps<=2);
        assert_float_eq!(extremes.z1, 6.0, ulps<=2);
        assert_float_eq!(extremes.t , 7.0, ulps<=2);
    }

    #[test]
    fn dataset_writer_flushes_on_drop() {
        let (_dir, filename) = tempfile("flush_on_drop");
        let file             = hdf5::File::create(filename.clone()).unwrap();
        let writer           = Hdf5Writer::<Endpoint>::new(Rc::new(file), "a_dataset", 10).unwrap();

        let data_written = extremes().0;
        writer.write(data_written).unwrap();
        drop(writer);

        let data_read = read_hdf5::<Endpoint>(&filename, "a_dataset").unwrap();

        assert_eq!(data_read.len(), 1);

        let extremes = data_read.first().unwrap();
        assert_eq!      (extremes.event, 0);
        assert_float_eq!(extremes.x0, 1.0, ulps<=2);
        assert_float_eq!(extremes.y0, 2.0, ulps<=2);
        assert_float_eq!(extremes.z0, 3.0, ulps<=2);
        assert_float_eq!(extremes.x1, 4.0, ulps<=2);
        assert_float_eq!(extremes.y1, 5.0, ulps<=2);
        assert_float_eq!(extremes.z1, 6.0, ulps<=2);
        assert_float_eq!(extremes.t , 7.0, ulps<=2);
    }

}
