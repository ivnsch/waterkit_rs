#### Rust test (standalone)

The parameters are hardcoded in JSON / the test. This is just during development, to not have start from Python.

```
cargo test --release -- --nocapture
```

#### Build Rust to be used by python

Install [maturin](https://github.com/PyO3/maturin?tab=readme-ov-file#usage)

Activate the waterkit environment, e.g:

```
conda activate waterkit
```

Then run:

```
maturin develop
```

#### Run Rust from python

```
conda activate waterkit
```

```
import waterkit_rs

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

print(str(await waterkit_rs.run(receptor_for_rust, water_for_rust, ad_map_for_rust)))
```

Currently requires waterkit/data directory to be directly under root, and files receptor_maps.fld, receptor.d.map, receptor.e.map and receptor.OW.map to be under root as well.
