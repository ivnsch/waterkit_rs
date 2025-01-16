use vek::Vec3;

use crate::atom::{Atom, HydrogenBond, Molecule};

#[derive(Debug, Clone)]
pub struct Water {
    pub atoms: Vec<Atom>,
    pub hb_anchor: Vec3<f32>,
    pub hb_vector: Vec3<f32>,
    pub hb_type: String,

    pub coordinates: Vec<Vec3<f32>>,
    pub molecule: Molecule,
}

impl Water {
    // xyz (array_like): 3d coordinates of the oxygen water molecule
    // atom_type (str): atom types of the spherical water molecule (default: OW)
    // partial_charge (float): partial charge of the spherical water molecule (default: -0.834)
    // hb_anchor (array_like): 3d coordinates of the HB anchor (Default: [0, 0, 0])
    // hb_vector (array_like): 3d coordinates of the HB vector (Default: xyz variable)
    // hb_type (str): type of the HB anchor (acceptor or donor)
    pub fn new(
        xyz: &Vec3<f32>,
        atom_type: &str,
        partial_charge: f32,
        hb_anchor: Option<Vec3<f32>>,
        hb_vector: Option<Vec3<f32>>,
        hb_type: Option<String>,
    ) -> Water {
        let hb_anchor = hb_anchor.unwrap_or_else(|| Vec3 {
            x: 0.,
            y: 0.,
            z: 0.,
        });
        let hb_vector = hb_vector.unwrap_or_else(|| xyz.clone());
        let hb_type = hb_type.unwrap_or_else(|| "acceptor".to_string());

        // Add the oxygen atom

        Water {
            atoms: vec![create_atom(0, xyz, atom_type, partial_charge)],
            hb_anchor,
            hb_vector,
            hb_type: hb_type.clone(),
            // TODO port: inheritance: trait? this composition not looking good
            coordinates: vec![*xyz],
            molecule: Molecule {
                atoms: vec![],
                hydrogen_bonds: Some(vec![]),
                rotatable_bonds: vec![],
                coordinates: vec![*xyz],
                hb_anchor,
                hb_vector,
                hb_type,
            },
        }
    }

    fn add_atom(&mut self, xyz: &Vec3<f32>, atom_type: &String, partial_charge: f32) {
        let atom = create_atom(self.atoms.len(), xyz, atom_type, partial_charge);
        self.atoms.push(atom);
    }

    pub fn build_explicit_water(&mut self, water_model: &str) -> bool {
        todo!()
    }

    pub fn update_coordinates(&mut self, xyz: &Vec3<f32>, atom_id: usize) -> bool {
        self.molecule.update_coordinates(xyz, atom_id)
    }

    pub fn atom_types(&self) -> Vec<String> {
        self.molecule.atom_types()
    }

    pub fn translate(&mut self, vector: &Vec3<f32>) {
        self.molecule.translate(vector)
    }

    /// port: this is really just "atoms" as we don't separate the fields (for now?)
    /// Get atom informations (xyz, q, type).
    // Args:
    //     atom_ids (int, list): index of one or multiple atoms
    // Returns:
    //     ndarray: atom information (i, xyz, q, t)
    pub fn atom_informations(&self, atom_ids: Option<Vec<usize>>) -> Vec<Atom> {
        self.molecule.atom_informations(atom_ids)
    }

    pub fn hydrogen_bonds(&self) -> Option<Vec<HydrogenBond>> {
        self.molecule.hydrogen_bonds.clone()
    }
}

fn create_atom(index: usize, xyz: &Vec3<f32>, atom_type: &str, partial_charge: f32) -> Atom {
    Atom {
        index,
        name: atom_type.to_string(),
        resname: "HOH".to_string(),
        resnum: 1,
        coords: *xyz,
        q: partial_charge,
        t: atom_type.to_string(),
    }
}
