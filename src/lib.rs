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
use serde::{Deserialize, Serialize};

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

// python:
// receptor_atoms_tuples = [tuple(atom) for atom in receptor.atoms]
// hydrogen_bond_tuples = list(receptor.hydrogen_bonds.itertuples(index=False, name=None))
// rotatable_bond_tuples = list(receptor.rotatable_bonds.itertuples(index=False, name=None))
// molecule_for_rust = {
//     "atoms": receptor_atoms_tuples,
//     "hydrogen_bonds": hydrogen_bond_tuples,
//     "rotatable_bonds": rotatable_bond_tuples
// }
// string_sum.process_molecule(molecule_for_rust)
#[pyfunction]
fn process_molecule(input: PythonMolecule) -> PyResult<()> {
    println!("got molecule from python!: {:?}", input);
    let res = save_molecule_json(&input, "./myreceptor.json");

    if let Err(e) = res {
        println!("e: {:?}", e);
    }

    Ok(())
}

#[derive(Debug, FromPyObject, Serialize, Deserialize)]
struct PythonMolecule {
    #[pyo3(item)]
    atoms: Vec<(i32, String, String, i32, String, [f32; 3], f32, String)>,
    #[pyo3(item)]
    hydrogen_bonds: Vec<(i32, [f32; 3], String, String)>,
    #[pyo3(item)]
    rotatable_bonds: Vec<(
        i32,
        i32,
        [f32; 3],
        [f32; 3],
        [f32; 3],
        [f32; 3],
        Vec<i32>,
        String,
    )>,
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

use std::{
    fs::File,
    io::{Read, Write},
};

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
    m.add_function(wrap_pyfunction!(process_molecule, m)?)?;
    Ok(())
}

fn save_molecule_json(molecule: &PythonMolecule, filename: &str) -> std::io::Result<()> {
    let json_string = serde_json::to_string_pretty(molecule).unwrap();
    let mut file = File::create(filename)?;
    file.write_all(json_string.as_bytes())?;
    Ok(())
}

fn load_molecule_json(filename: &str) -> std::io::Result<PythonMolecule> {
    let mut file = File::open(filename)?;
    let mut buffer = String::new();
    file.read_to_string(&mut buffer)?;
    let molecule: PythonMolecule = serde_json::from_str(&buffer).unwrap();
    Ok(molecule)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn run_waterkit() {
        // to not have to always start rust with python during development, we serialize some parameters and use that
        let receptor_res = load_molecule_json("./dev_pars/myreceptor.json");
        assert!(receptor_res.is_ok());
        let receptor = receptor_res.unwrap();
        println!("atoms: {}", receptor.atoms.len());
        println!("hbs: {}", receptor.hydrogen_bonds.len());
        println!("rbs: {}", receptor.rotatable_bonds.len());
    }
}
