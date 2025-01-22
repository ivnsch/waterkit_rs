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
use anyhow::Result;
use atom::{Atom, Bond, Molecule};
use autodock_map::Map;
use hydrogen_bonds::HydrogenBond;
use run_waterkit::hydrate_rust;
use serde::{Deserialize, Serialize};
use vek::Vec3;

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
fn process_ad_map(input: PythonAdMap) -> PyResult<()> {
    println!("got ad map from python!: {:?}", input);
    Ok(())
}

// python:
// receptor_atoms_tuples = [tuple(atom) for atom in receptor.atoms]
// hydrogen_bond_tuples = list(receptor.hydrogen_bonds.itertuples(index=False, name=None))
// rotatable_bond_tuples = list(receptor.rotatable_bonds.itertuples(index=False, name=None))
// receptor_for_rust = {
//     "atoms": receptor_atoms_tuples,
//     "hydrogen_bonds": hydrogen_bond_tuples,
//     "rotatable_bonds": rotatable_bond_tuples
// }
// (same for water)
// string_sum.process_molecules(receptor_for_rust, water_for_rust)
#[pyfunction]
fn save_parameters(
    receptor: PythonMolecule,
    water: PythonMolecule,
    ad_map: PythonAdMap,
) -> PyResult<()> {
    let receptor_res = save_molecule_json(&receptor, "./myreceptor.json");
    let water_res = save_molecule_json(&water, "./mywater_molecule.json");
    let map_res = save_ad_map_json(&ad_map, "./my_ad_map.json");

    println!("water from python: {:?}", water);

    if let Err(e) = receptor_res {
        println!("e with receptor: {:?}", e);
    }
    if let Err(e) = water_res {
        println!("e with water: {:?}", e);
    }
    if let Err(e) = map_res {
        println!("e with ad map: {:?}", e);
    }

    Ok(())
}

///
///
#[pyfunction]
async fn run(
    receptor_python: PythonMolecule,
    water_python: PythonMolecule,
    ad_map_python: PythonAdMap,
) -> PyResult<()> {
    let receptor = to_molecule(receptor_python);
    let water = to_molecule(water_python);
    let ad_map = to_map(ad_map_python).await.unwrap();
    let res = hydrate_rust(receptor, water, ad_map, "traj", 10000, 3, 20.).await;
    Ok(())
}

#[derive(Debug, FromPyObject, Serialize, Deserialize)]
struct PythonMolecule {
    #[pyo3(item)]
    atoms: Vec<(usize, String, String, u32, String, [f32; 3], f32, String)>,
    #[pyo3(item)]
    hydrogen_bonds: Vec<(usize, [f32; 3], String, String)>,
    #[pyo3(item)]
    rotatable_bonds: Vec<(
        usize,
        usize,
        [f32; 3],
        [f32; 3],
        [f32; 3],
        [f32; 3],
        Vec<i32>,
        String,
    )>,
}

#[derive(Debug, FromPyObject, Serialize, Deserialize)]
struct PythonAdMap {
    #[pyo3(item)]
    box_center: [f32; 3],
    #[pyo3(item)]
    box_size: [usize; 3],
    #[pyo3(item)]
    box_spacing: f32,
    #[pyo3(item)]
    labels: Vec<String>,
    #[pyo3(item)]
    files: Vec<String>,
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
fn waterkit_rs(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(sum_as_string, m)?)?;
    m.add_function(wrap_pyfunction!(hello_rust, m)?)?;
    m.add_function(wrap_pyfunction!(process_atom, m)?)?;
    m.add_function(wrap_pyfunction!(process_atom_list, m)?)?;
    m.add_function(wrap_pyfunction!(process_hydrogen_bonds, m)?)?;
    m.add_function(wrap_pyfunction!(process_rotatable_bonds, m)?)?;
    m.add_function(wrap_pyfunction!(save_parameters, m)?)?;
    m.add_function(wrap_pyfunction!(process_ad_map, m)?)?;
    m.add_function(wrap_pyfunction!(run, m)?)?;
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

fn save_ad_map_json(ad_map: &PythonAdMap, filename: &str) -> std::io::Result<()> {
    let json_string = serde_json::to_string_pretty(ad_map).unwrap();
    let mut file = File::create(filename)?;
    file.write_all(json_string.as_bytes())?;
    Ok(())
}

fn load_ad_map_json(filename: &str) -> std::io::Result<PythonAdMap> {
    let mut file = File::open(filename)?;
    let mut buffer = String::new();
    file.read_to_string(&mut buffer)?;
    let molecule: PythonAdMap = serde_json::from_str(&buffer).unwrap();
    Ok(molecule)
}

pub fn to_molecule(python_mol: PythonMolecule) -> Molecule {
    let atoms = python_mol
        .atoms
        .into_iter()
        .map(|tuple| Atom {
            index: tuple.0,
            name: tuple.1.clone(),
            resname: tuple.2.clone(),
            resnum: tuple.3,
            // tuple.4 not used?
            coords: tuple.5.into(),
            q: tuple.6,
            t: tuple.7.clone(),
        })
        .collect::<Vec<Atom>>();

    let hydrogen_bonds = Some(
        python_mol
            .hydrogen_bonds
            .into_iter()
            .map(|tuple| {
                // (2, array([17.14364591, 39.8412896 , -4.20977717]), 'donor', 'H_N_000')
                HydrogenBond {
                    atom_i: tuple.0,
                    vector_xyz: tuple.1.into(),
                    anchor_type: tuple.2,
                    atom_types: tuple.3,
                }
            })
            .collect::<Vec<HydrogenBond>>(),
    );

    let rotatable_bonds = python_mol
        .rotatable_bonds
        .into_iter()
        .map(|tuple| Bond {
            atom_i: tuple.0,
            atom_j: Some(tuple.1),
            // these 2 fields seem to be set later, in specific context
            // TODO port review them and whether initialization to 0 is right
            molecule_i: 0,
            molecule_j: 0,
            atom_i_xyz: tuple.2.into(),
            atom_j_xyz: tuple.3.into(),
            atom_k_xyz: tuple.4.into(),
            atom_l_xyz: tuple.5.into(),
        })
        .collect::<Vec<Bond>>();
    Molecule {
        atoms,
        hydrogen_bonds,
        rotatable_bonds,
        // TODO port: where to these come from, probably have to be passed from python too now
        hb_anchor: Vec3::zero(),
        hb_vector: Vec3::zero(),
        hb_type: "".to_string(),
    }
}

async fn to_map(python_map: PythonAdMap) -> Result<Map> {
    Map::new(python_map.files, python_map.labels).await
}

#[cfg(test)]
mod tests {
    use run_waterkit::hydrate_rust;

    use super::*;
    #[tokio::test]
    async fn run_waterkit() {
        // to not have to always start rust with python during development, we serialize some parameters and use that
        let receptor_res: Result<PythonMolecule, std::io::Error> =
            load_molecule_json("./dev_pars/myreceptor.json");
        assert!(receptor_res.is_ok());
        let receptor_python = receptor_res.unwrap();
        println!("receptor atoms: {}", receptor_python.atoms.len());
        println!("receptor hbs: {}", receptor_python.hydrogen_bonds.len());
        println!("receptor rbs: {}", receptor_python.rotatable_bonds.len());
        let receptor = to_molecule(receptor_python);

        let water_res: Result<PythonMolecule, std::io::Error> =
            load_molecule_json("./dev_pars/mywater_molecule.json");
        assert!(water_res.is_ok());
        let water_python = water_res.unwrap();
        println!("water atoms: {}", water_python.atoms.len());
        println!("water hbs: {}", water_python.hydrogen_bonds.len());
        println!("water rbs: {}", water_python.rotatable_bonds.len());
        let water = to_molecule(water_python);

        let ad_map_res = load_ad_map_json("./dev_pars/my_ad_map.json");
        assert!(ad_map_res.is_ok());
        let ad_map_python = ad_map_res.unwrap();
        let ad_map = to_map(ad_map_python).await.unwrap();

        let res = hydrate_rust(receptor, water, ad_map, "traj", 10000, 3, 20.).await;
        println!("res: {:?}", res);
    }
}
