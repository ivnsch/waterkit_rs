use vek::Vec3;

use crate::{hydrogen_bonds::HydrogenBond, water::Water};

#[derive(Debug, Clone)]
pub struct Atom {
    pub index: usize,
    pub name: String,
    pub resname: String,
    pub resnum: u32,
    pub coords: Vec3<f32>,
    pub q: f32,
    pub t: String,
}

// is this really needed? using Atom instead for now
// #[derive(Debug, Clone)]
// pub struct AtomInformations {
//     pub index: usize,
//     pub coords: Vec3<f32>,
//     pub q: f32,
//     pub t: String,
// }

#[derive(Debug, Clone)]
pub enum MoleculeType {
    Receptor(Molecule),
    Water(Water),
}

impl MoleculeType {
    pub fn atoms(&self) -> Vec<Atom> {
        match self {
            MoleculeType::Receptor(Molecule { atoms, .. }) => atoms.clone(),
            MoleculeType::Water(water) => water.atoms().to_owned(),
        }
    }

    pub fn hydrogen_bonds(&self) -> Option<Vec<HydrogenBond>> {
        match self {
            MoleculeType::Receptor(m) => m.hydrogen_bonds.clone(),
            MoleculeType::Water(w) => w.hydrogen_bonds(),
        }
    }

    pub fn coordinates(&self, atom_ids: Option<&[usize]>) -> Vec<Vec3<f32>> {
        match self {
            MoleculeType::Receptor(m) => m.coordinates(atom_ids).clone(),
            MoleculeType::Water(w) => w.coordinates(atom_ids).clone(),
        }
    }
}

#[derive(Debug, Clone)]
pub struct Molecule {
    pub atoms: Vec<Atom>,
    pub hydrogen_bonds: Option<Vec<HydrogenBond>>,
    pub rotatable_bonds: Vec<Bond>,

    pub hb_anchor: Vec3<f32>,
    pub hb_vector: Vec3<f32>,

    pub hb_type: String,
}

impl Molecule {
    /// Return coordinates of all atoms or a certain atom
    /// Args:
    ///     atom_ids (int, list): index of one or multiple atoms
    /// Returns:
    ///     ndarray: 2d ndarray of 3d coordinates
    pub fn coordinates(&self, atom_ids: Option<&[usize]>) -> Vec<Vec3<f32>> {
        if let Some(atom_ids) = atom_ids {
            // -1 because we use 0 based indices
            let atom_ids_minus_1: Vec<usize> = atom_ids.iter().map(|a| a - 1).collect();
            self.atoms
                .iter()
                .enumerate()
                .filter(|(index, _)| atom_ids_minus_1.contains(&index))
                .map(|(_, a)| a.coords)
                .collect()
        } else {
            self.atoms.iter().map(|a| a.coords).collect()
        }
    }

    /// Update the coordinates of an atom.
    /// Args:
    ///     xyz (array_like): 3d coordinates of the new atomic position
    ///     atom_id (int): atom id
    /// Returns:
    ///     bool: True if successfull, False otherwise
    pub fn update_coordinates(&mut self, xyz: &Vec3<f32>, atom_id: usize) -> bool {
        if atom_id <= self.atoms.len() {
            if self.atoms.len() > 1 {
                self.atoms[atom_id - 1].coords = *xyz;
            } else if let Some(atom) = self.atoms.first_mut() {
                atom.coords = *xyz;
            }
            return true;
        }
        false
    }

    /// Return atom types of all atoms or a certain atom.
    /// Args:
    ///     atom_ids (int, list): index of one or multiple atoms
    /// Returns:
    ///     list: atom types
    pub fn atom_types(&self, atom_ids: Option<&[usize]>) -> Vec<String> {
        if let Some(atom_ids) = atom_ids {
            if self.atoms.len() > 1 {
                // -1 because we use 0 based indices
                let atom_ids_minus_1: Vec<usize> = atom_ids.iter().map(|a| a - 1).collect();
                self.atoms
                    .iter()
                    .enumerate()
                    .filter(|(index, _)| atom_ids_minus_1.contains(&index))
                    .map(|(_, a)| a.t.clone())
                    .collect()
            } else {
                self.atoms.iter().map(|a: &Atom| a.t.clone()).collect()
            }
        } else {
            self.atoms.iter().map(|a: &Atom| a.t.clone()).collect()
        }
    }

    /// Translate the water molecule.
    /// Args:
    ///     vector (array_like): 3d vector
    /// Returns:
    ///     None
    pub fn translate(&mut self, vector: &Vec3<f32>) {
        let water_xyz: Vec<Vec3<f32>> = self.coordinates(None).iter().map(|c| c + vector).collect();

        for (atom_id, coord) in water_xyz.iter().enumerate() {
            // +1, because atom ids are 1-based
            self.update_coordinates(coord, atom_id + 1);
        }
        // We have also to translate the hydrogen bond vectors if present
        if let Some(hydrogen_bonds) = self.hydrogen_bonds.as_mut() {
            for bond in hydrogen_bonds {
                bond.vector_xyz += *vector;
            }
        }
    }

    /// port: this is really just "atoms" as we don't separate the fields (for now?)
    /// Get atom informations (xyz, q, type).
    // Args:
    //     atom_ids (int, list): index of one or multiple atoms
    // Returns:
    //     ndarray: atom information (i, xyz, q, t)
    pub fn atom_informations(&self, atom_ids: Option<Vec<usize>>) -> Vec<Atom> {
        if let Some(atom_ids) = atom_ids {
            // # -1 because numpy array is 0-based
            let atom_ids = atom_ids.iter().map(|a| a - 1).collect::<Vec<usize>>();
            atom_ids.iter().map(|&i| self.atoms[i].clone()).collect()
        } else {
            self.atoms.clone()
        }
    }
}

#[derive(Debug, Clone)]
pub struct Bond {
    pub atom_i: usize,
    pub atom_j: Option<usize>,
    pub molecule_i: usize,
    pub molecule_j: usize,
    // these should be options too, sometimes using defaults (zero) instead
    pub atom_i_xyz: Vec3<f32>,
    pub atom_j_xyz: Vec3<f32>,
    pub atom_k_xyz: Vec3<f32>,
    pub atom_l_xyz: Vec3<f32>,
}

impl Bond {
    // TODO port: better name for this
    pub fn dihedral_pars(&self) -> DihedralPars {
        DihedralPars {
            atom_i_xyz: self.atom_i_xyz,
            atom_j_xyz: self.atom_j_xyz,
            atom_k_xyz: self.atom_k_xyz,
            atom_l_xyz: self.atom_l_xyz,
        }
    }
}

// TODO port: better name for this
pub struct DihedralPars {
    pub atom_i_xyz: Vec3<f32>,
    pub atom_j_xyz: Vec3<f32>,
    pub atom_k_xyz: Vec3<f32>,
    pub atom_l_xyz: Vec3<f32>,
}
