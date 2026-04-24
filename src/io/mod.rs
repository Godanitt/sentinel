mod csv;
mod hdf5;
pub mod hdf5_types;

pub use csv::read_csv;
pub use csv::CsvWriter;

pub use hdf5::read_hdf5;
pub use hdf5::ScalarHdf5Writer;
pub use hdf5:: ArrayHdf5Writer;

use std::io;
pub trait Writer<T> {
    fn write(&self, value: T) -> io::Result<()>;

    fn write_many(&self, values: Vec<T>) -> io::Result<()> {
        for v in values.into_iter() {
            self.write(v)?;
        }
        Ok(())
    }
}
