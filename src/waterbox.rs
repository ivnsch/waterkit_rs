use anyhow::Result;
use kd_tree::{KdPoint, KdTree};
use vek::Vec3;

use crate::{
    atom::{Bond, Molecule, MoleculeType},
    autodock_map::Map,
    sampling::WaterSampler,
    water::Water,
};

#[derive(Debug, Clone)]
pub struct WaterBox {
    pub molecules: Vec<MoleculeType>,
    // pub receptor: Molecule,
    pub map: Map,
    pub kd_tree: KdTree<[f32; 3]>,
    pub water_model: String,
    pub temperature: f32,
    pub wopt: WaterSampler,

    // df
    pub data: BoxData,
}

#[derive(Debug, Clone)]
pub struct BoxData {
    pub connections: Vec<Bond>,
    pub shells: Vec<Shell>,
    pub kdtree_relations: Vec<KdTreeRelation>,
}

impl WaterBox {
    pub async fn new(
        receptor: Molecule,
        map: Map,
        temperature: f32,
        water_model: &str,
        spherical_water_map: Option<&str>,
    ) -> Result<WaterBox> {
        let mut wopt = WaterSampler::new(
            map.clone(),
            temperature,
            spherical_water_map,
            2.5,
            3.6,
            90.,
            0.,
            water_model,
        )
        .await?;

        let water_box = WaterBox {
            molecules: vec![MoleculeType::Receptor(receptor)],
            map,
            kd_tree: KdTree::build_by_ordered_float(vec![]),
            water_model: water_model.to_string(),
            temperature: 0.,
            wopt: wopt.clone(),
            data: BoxData {
                connections: vec![],
                shells: vec![],
                kdtree_relations: vec![],
            },
        };

        wopt.set_water_box(water_box.clone());

        Ok(water_box)
    }

    fn add_receptor(&mut self, receptor: Molecule) {
        self.add_molecules(&vec![MoleculeType::Receptor(receptor)], true);
    }

    // TODO port: not sure why python autoincrementing dictionary key? using just a vec for now
    // TODO port: add_kd_tree not used?
    fn add_molecules(&mut self, molecules: &[MoleculeType], add_kd_tree: bool) {
        self.molecules.extend(molecules.iter().cloned());
    }

    /// Add connections between molecules
    fn add_connections(&mut self, connections: &mut [Bond]) {
        let (last_molecule_i, last_molecule_j) =
            if let Some(last_connection) = self.data.connections.last() {
                (
                    last_connection.molecule_i as i32,
                    last_connection.molecule_j as i32,
                )
            } else {
                (-1, 0)
            };

        for connection in connections.iter_mut() {
            connection.molecule_i += (last_molecule_i + 1) as usize;
            connection.molecule_j += (last_molecule_j + 1) as usize;
        }
        // TODO port: confirm, is this what add informations is doing?
        // self._add_informations(connections, "connections")
        self.data.connections.extend(connections.iter().cloned());
    }

    fn add_molecules_to_kd_tree(&mut self, molecules: &[MoleculeType]) {
        let mut points = self.current_kd_tree_points();
        let new_points = molecules_as_points(molecules);
        points.extend(new_points);

        let points_arr = points.into_iter().map(|p| p.into_array()).collect();
        let new_kd_tree = KdTree::build_by_ordered_float(points_arr);

        self.kd_tree = new_kd_tree;
    }

    fn closest_atoms(&self, point: &Vec3<f32>, radius: f32) -> Vec<Vec3<f32>> {
        let arr = self.kd_tree.within_radius(&point.into_array(), radius);
        arr.into_iter().map(|r| Vec3::from(*r)).collect()
    }

    fn current_kd_tree_points(&self) -> Vec<Vec3<f32>> {
        // NOTE: assumes self.kd_tree reflects self.molecules, i.e. self.molecules are in self.kd_tree
        // hmm is there really nothing to get the points in a kd tree or otherwise something to rebuild it?
        molecules_as_points(&self.molecules)
    }

    /// Place spherical water molecules in the optimal position.
    /// Args:
    ///     molecules (Molecule): molecules on which spherical water molecules will be placed
    ///     atom_type (str): atom type of the spherical water molecules (default: OW)
    ///     partial_charges (float): partial charges of the spherical water molecules (default: -0.834)
    /// Returns:
    ///     list: list of water molecules
    ///     DataFrame: contains connections between molecules
    /// TODO port: return connections
    fn place_optimal_spherical_waters(
        &self,
        molecules: &[MoleculeType],
        atom_type: &str,
        partial_charge: f32,
    ) -> Vec<Water> {
        let mut waters = vec![];

        for molecule in molecules {
            if let Some(bonds) = &molecule.hydrogen_bonds() {
                for bond in bonds {
                    // Add water molecule only if it's in the map
                    // TODO what structure is this? for now removing the [0]
                    // let anchor_xyz = molecule.coordinates[bond.atom_i][0];
                    let anchor_xyz = molecule.coordinates()[bond.atom_i];

                    if self.map.is_in_map(anchor_xyz) {
                        let w = Water::new(
                            &bond.vector_xyz,
                            atom_type,
                            partial_charge,
                            Some(anchor_xyz),
                            Some(bond.vector_xyz),
                            Some(bond.anchor_type.clone()),
                        );
                        waters.push(w);
                    }
                }
            }
        }
        waters
    }

    /// Build the next hydration shell.
    /// Returns:
    ///     bool: True if water molecules were added or False otherwise
    pub fn build_next_shell(&mut self) -> bool {
        let sw_type = "SW";
        let partial_charge = 0.0;

        let shell_id = self.data.shells.len();
        let molecules = self.molecules_in_shell(&[]);

        let mut waters = self.place_optimal_spherical_waters(&molecules, sw_type, partial_charge);

        // Only the receptor contains disordered hydrogens
        let (w, df) = if shell_id == 0 {
            // TODO port: connections
            self.wopt.sample_grid(&mut waters, &mut vec![], true)
        } else {
            // After the first hydration layer, we don't care anymore about
            // connections. It was only useful for the disordered hydrogen atoms.
            // TODO port: connections: note that here (opposed to if above) we actually don't want to pass any, should we make it an option
            self.wopt.sample_grid(&mut waters, &mut vec![], true)
        };
        waters = w;

        if !waters.is_empty() {
            // TODO port: wrap in enum: performance?
            let molecules = waters
                .into_iter()
                .map(|w| MoleculeType::Water(w))
                .collect::<Vec<MoleculeType>>();
            self.add_molecules(&molecules, false);
            // Add informations about the new shell
            // TODO port: review, is this the same?
            if !self.data.shells.is_empty() {
                self.data.shells.extend(df.shells);
            }
            // if "shells" in df.keys():
            // self._add_informations(df["shells"], "shells")
            true
        } else {
            false
        }
    }

    /// Get all the molecule in shell.
    /// Args:
    ///     shell_ids (list): ids of the shell(s) (default: None)
    /// Returns:
    ///     list: list of all the molecules in the selected shell(s)
    pub fn molecules_in_shell(&self, shell_ids: &[usize]) -> Vec<MoleculeType> {
        // TODO port: filter by ids
        // println!("self.data.shells: {:?}", self.data.shells);
        // self.data
        //     .shells
        //     .iter()
        //     .flat_map(|s| s.molecules.clone())
        //     .collect()

        self.molecules.clone()
    }

    pub fn number_of_shells(&self) -> usize {
        self.data.shells.len()
    }

    /// Write all the content of the water box in a PDB file.
    /// We cannot use OpenBabel to write the PDB file of water molecules
    /// because it is always changing the atom types...
    /// Args:
    ///     fname (str): name of the output file
    ///     water_model (str): water model we want to use to export the water (default: tip3p)
    ///     include_receptor (bool): include or exclude the receptor (default: True)
    /// Returns:
    ///     None
    pub fn to_pdb(&self, fname: &str, water_model: &str, include_receptor: bool) {
        let mut i = 0;
        let mut j = 0;
        let mut output_str = "".to_string();
        let pdb_str = "ATOM  %5d  %-4s%-3s%2s%4d    %8.3f%8.3f%8.3f  1.00  1.00          %2s\n";

        let shell_id = self.number_of_shells();
        let water_shells = (1..=shell_id).map(|i| self.molecules_in_shell(&vec![shell_id]));

        if include_receptor {
            let atoms = self.molecules[0].atoms();
            for atom in atoms {
                // TODO port: hmm why "xyz", probably different to coords? python: x, y, z = atom["xyz"]
                // let [x, y, z] = atom.coords;
                // TODO port what's this?
                // output_str += pdb_str % (atom["i"], atom["name"], atom["resname"], "A", atom["resnum"], x, y, z, atom["t"])
                // output_str = format!("{}{}", output_str, pdb_str);
            }

            // TODO port iterate through this
            // let ascii_uppercase = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
            // for (water_shell, chain) in water_shells.zip(ascii_uppercase.into) {

            // TODO port: mostly more pdb_str % ...
            // }

            // TODO port: what is this?
            // Get the index of the nex atom and residue
            // i = atom["i"] + 1
            // j = atom["resnum"] + 1
        }
    }
}

#[derive(Debug, Clone)]
pub struct Shell {
    pub molecules: Vec<Molecule>,
    pub shell_id: usize,
    pub energy_position: f32,
    pub energy_orientation: Vec<f32>,
}

#[derive(Debug, Clone)]
pub struct KdTreeRelation {
    pub atom_i: usize,
    pub molecule_i: usize,
}

#[cfg(test)]
mod test {
    #[test]
    fn load_pdb() {
        let (mut pdb, _errors) = pdbtbx::open("protein_prepared_amber.pdb").unwrap();
    }
}

fn molecules_as_points(molecules: &[MoleculeType]) -> Vec<Vec3<f32>> {
    let mut points = vec![];
    for m in molecules {
        for a in &m.atoms() {
            points.push(a.coords);
        }
    }
    points
}
