use crate::{atom::Molecule, autodock_map::Map, waterbox::WaterBox, waterkit::hydrate_single};
use anyhow::Result;
use rayon::iter::{IntoParallelIterator, ParallelIterator};

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
