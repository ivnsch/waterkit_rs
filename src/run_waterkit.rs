use crate::{
    atom::{self, Molecule},
    autodock_map::Map,
    waterbox::WaterBox,
    waterkit::hydrate_single,
};
use anyhow::Result;
use pdbtbx::PDB;
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use vek::Vec3;

// TODO async + pyfunction: "experimental feature"
// #[pyfunction]
pub async fn hydrate_rust(
    receptor: Molecule,
    water_molecule: Molecule,
    map: Map,
    output_dir: &str,
    n_frames: usize,
    n_layer: usize,
    temperature: f32,
) -> Result<()> {
    println!(
        "called hydrate: output_dir: {}, n_frames: {}, n_layer: {}",
        output_dir, n_frames, n_layer
    );

    let water_box_res =
        WaterBox::new(receptor, water_molecule, map, temperature, "tip3p", None).await;

    println!("will iterate over {} frames", n_frames);
    match water_box_res {
        Ok(water_box) => {
            (0..n_frames).into_par_iter().for_each(|i| {
                let mut water_box = water_box.clone();
                hydrate_single(&mut water_box, n_layer, i, output_dir);
            });
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
        // TODO port
        hb_anchor: Vec3::zero(),
        // TODO port
        hb_vector: Vec3::zero(),
        // TODO port
        hb_type: "".to_string(),
    }
}
