use crate::waterbox::WaterBox;

pub fn hydrate_single(water_box: &mut WaterBox, n_layer: usize, frame_id: usize, output_dir: &str) {
    // Single job waterkit hydration.
    let mut i = 1;

    loop {
        // build_next_shell returns True if
        // it was able to put water molecules,
        // otherwise it returns False and we break

        if !water_box.build_next_shell() {
            break;
        }
        // Stop after building n_layer or when the number of
        // layers is equal to 26. In the PDB format, the segid
        // is encoded with only one caracter, so the maximum number
        // of layers is 26, which corresponds to the letter Z.
        if i == n_layer || i == 26 {
            break;
        }

        i += 1;
    }

    let output_filename = format!("{}/water_{:06}.pdb", output_dir, frame_id + 1);
    water_box.to_pdb(&output_filename, "tip3p", false);
}
