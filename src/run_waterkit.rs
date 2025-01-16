use crate::{
    atom::{self, Molecule},
    autodock_map::Map,
    fld_parser::parse_fld,
    waterbox::WaterBox,
    waterkit::hydrate_single,
};
use anyhow::Result;
use pdbtbx::PDB;
use pyo3::pyfunction;
use vek::Vec3;

// TODO async + pyfunction: "experimental feature"
// #[pyfunction]
pub async fn hydrate_rust(
    pdbqt_file: &str,
    fld_file: &str,
    output_dir: &str,
    n_frames: usize,
    n_layer: usize,
    temperature: f32,
) -> Result<()> {
    println!(
        "called hydrate: pdbqt_file: {}, fld_path: {}, output_dir: {}, n_frames: {}, n_layer: {}",
        pdbqt_file, fld_file, output_dir, n_frames, n_layer
    );

    // let (mut pdbqt, _errors) = pdbtbx::open("protein_prepared_amber.pdbqt").unwrap();
    let (mut pdbqt, _errors) = pdbtbx::open(pdbqt_file).unwrap();
    let molecule = to_molecule(pdbqt);

    let parsed_fld = parse_fld(fld_file).await?;
    let map = Map::new(parsed_fld.map_files, parsed_fld.labels).await?;

    // TODO port: spherical_water_map

    let mut water_box_res = WaterBox::new(map, temperature, "tip3p", None).await;

    println!("will iterate over {} frames", n_frames);
    match water_box_res {
        Ok(mut water_box) => {
            // for now single threaded
            for i in 0..n_frames {
                hydrate_single(&mut water_box, n_layer, i, output_dir);
            }
        }
        Err(e) => {
            println!("Error initializing water box: {:?}", e)
        }
    }

    Ok(())
}

pub fn to_molecule(pdbqt: PDB) -> Molecule {
    let mol_atoms = pdbqt
        .atoms()
        .map(|a| atom::Atom {
            index: a.serial_number(),
            name: a.name().to_string(),
            // TODO port what field is this?
            resname: "".to_string(),
            // TODO port what field is this?
            resnum: 0,
            coords: Vec3::new(a.x() as f32, a.y() as f32, a.z() as f32),
            q: a.charge(),
            t: a.autodock_type(),
        })
        .collect();
    Molecule {
        atoms: mol_atoms,
        hydrogen_bonds: None,
        rotatable_bonds: vec![],
        coordinates: vec![],
        // TODO port
        hb_anchor: Vec3::zero(),
        // TODO port
        hb_vector: Vec3::zero(),
        // TODO port
        hb_type: "".to_string(),
    }
}

#[cfg(test)]
mod tests {
    use super::hydrate_rust;

    #[tokio::test]
    async fn run() {
        let res = hydrate_rust(
            "protein_prepared_amber.pdbqt",
            "receptor_maps.fld",
            ".",
            100,
            123,
            20.,
        )
        .await;
        println!("res: {:?}", res);
    }
}
