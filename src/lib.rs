mod atom;
mod autodock_map;
mod fld_parser;
mod force_fields_parser;
mod hydrogen_bonds;
mod map_parser;
mod params_parser;
mod run_waterkit;
mod sampling;
mod utils;
mod water;
mod waterbox;
pub mod waterkit;

#[pyfunction]
fn process_atom(input: (i32, String, String, i32, String, [f32; 3], f32, String)) -> PyResult<()> {
    println!("got an atom from python!: {:?}", input);
    Ok(())
}

#[pyfunction]
fn process_atom_list(
    input: Vec<(i32, String, String, i32, String, [f32; 3], f32, String)>,
) -> PyResult<()> {
    println!("got atom list from python!: {:?}", input);
    Ok(())
}

#[pyfunction]
fn process_hydrogen_bonds(input: Vec<(i32, [f32; 3], String, String)>) -> PyResult<()> {
    println!("got hydrogen bond list from python!: {:?}", input);
    Ok(())
}

#[pyfunction]
fn process_rotatable_bonds(
    input: Vec<(
        i32,
        i32,
        [f32; 3],
        [f32; 3],
        [f32; 3],
        [f32; 3],
        Vec<i32>,
        String,
    )>,
) -> PyResult<()> {
    println!("got rotatable bond list from python!: {:?}", input);
    Ok(())
}

pub fn add(left: u64, right: u64) -> u64 {
    left + right
}

use pyo3::{prelude::*, types::PyTuple};

#[pyfunction]
fn hello_rust(name: &str) -> PyResult<String> {
    Ok(format!("hello {:?}!!", name))
}

/// Formats the sum of two numbers as string.
#[pyfunction]
fn sum_as_string(a: usize, b: usize) -> PyResult<String> {
    Ok((a + b).to_string())
}

/// A Python module implemented in Rust. The name of this function must match
/// the `lib.name` setting in the `Cargo.toml`, else Python will not be able to
/// import the module.
#[pymodule]
fn string_sum(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(sum_as_string, m)?)?;
    m.add_function(wrap_pyfunction!(hello_rust, m)?)?;
    m.add_function(wrap_pyfunction!(process_atom, m)?)?;
    m.add_function(wrap_pyfunction!(process_atom_list, m)?)?;
    m.add_function(wrap_pyfunction!(process_hydrogen_bonds, m)?)?;
    m.add_function(wrap_pyfunction!(process_rotatable_bonds, m)?)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_works() {
        let result = add(2, 2);
        assert_eq!(result, 4);
    }
}
