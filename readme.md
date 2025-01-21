run_waterkit (hardcoded parameters in JSON/test)

```
cargo test --release -- --nocapture
```

Python example (this doesn't yet start run_waterkit)

```
receptor_for_rust = {
    "atoms": receptor_atoms_tuples,
    "hydrogen_bonds": receptor_hydrogen_bond_tuples,
    "rotatable_bonds": receptor_rotatable_bond_tuples
}
water_for_rust = {
    "atoms": water_receptor_atoms_tuples,
    "hydrogen_bonds": water_hydrogen_bond_tuples,
    "rotatable_bonds": water_rotatable_bond_tuples
}
ad_map_for_rust = {
    "box_center": ad_map._center,
    "box_size": ad_map._npts,
    "box_spacing": ad_map._spacing,
    "labels": list(ad_map._maps.keys()),
    "files": ad_map._files,
}

string_sum.save_parameters(receptor_for_rust, water_for_rust, ad_map_for_rust)
```
