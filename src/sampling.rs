use anyhow::Result;
use nalgebra::{UnitQuaternion, Vector3};
use ndarray::{Array2, Array3, Axis};
use rand::Rng;
use vek::Vec3;

use crate::{
    atom::{Bond, Molecule, MoleculeType},
    autodock_map::Map,
    utils,
    water::Water,
    waterbox::{BoxData, Shell, WaterBox},
};

#[derive(Debug, Clone)]
pub struct WaterSampler {
    pub water_box: Option<Box<WaterBox>>,
    pub ad_map: Map,
    pub receptor: Option<Molecule>,
    pub rotation: f32,
    pub temperature: f32,
    pub max_distance: f32,
    pub min_distance: f32,
    pub energy_cutoff: f32,
    pub water_model: String,
    pub angle: f32,
    pub water_orientations: Vec<Vec<f32>>,
    pub water_ref: Molecule,
    pub water_map: Map,
}

impl WaterSampler {
    pub async fn new(
        ad_map: Map,
        temperature: f32,
        spherical_water_map: Option<&str>,
        min_distance: f32,
        max_distance: f32,
        angle: f32,
        energy_cutoff: f32,
        water_model: &str,
        water_molecule: Molecule,
    ) -> Result<WaterSampler> {
        // Water grids
        let mut map_list = vec![];
        let spherical_water_map =
            spherical_water_map.unwrap_or_else(|| "data/water/spherical/water_SW.map");
        map_list.push(spherical_water_map);

        // Explicit water grids
        let (atom_types, map_files) = match water_model {
            "tip3p" => (
                vec!["SW", "OW", "HW"],
                vec![
                    "data/water/tip3p/water_OW.map",
                    "data/water/tip3p/water_HW.map",
                ],
            ),
            "tip5p" => (
                vec!["SW", "OT", "HT", "LP"],
                vec![
                    "data/water/tip5p/water_OT.map",
                    "data/water/tip5p/water_HT.map",
                    "data/water/tip5p/water_LP.map",
                    "data/water/tip3p/water_HW.map",
                ],
            ),
            // TODO enum
            _ => panic!("not supported water model"),
        };
        map_list.extend(map_files);

        Ok(WaterSampler {
            // initialized via a setter, since circular reference
            water_box: None,
            ad_map,
            receptor: None,
            // TODO review these defaults
            rotation: 0.,
            temperature,
            max_distance,
            min_distance,
            energy_cutoff,
            water_model: water_model.to_string(),
            angle,
            water_orientations: vec![],
            water_ref: water_molecule,
            water_map: Map::new(
                map_list.into_iter().map(|s| s.to_string()).collect(),
                atom_types.into_iter().map(|s| s.to_string()).collect(),
            )
            .await?,
        })
    }

    pub fn set_water_box(&mut self, water_box: WaterBox) {
        let mols = water_box.molecules_in_shell(Some(&vec![0]));
        let receptor = match &mols[0] {
            MoleculeType::Receptor(molecule) => molecule,
            // TODO look for better way to structure this,
            // probably simply store separately and merge when they're needed in a list
            MoleculeType::Water(_) => panic!("we assume mol[0] to be the receptor.."),
        };

        self.receptor = Some(receptor.clone());

        self.water_model = water_box.water_model.clone();

        self.water_box = Some(Box::new(water_box));
    }

    pub fn sample_grid(
        &mut self,
        waters: &[Water],
        connections: Option<&mut [Bond]>,
        opt_disordered: bool,
    ) -> (Vec<Water>, BoxData) {
        let mut waters = waters.to_vec();
        // Optimize position of water molecules.
        let mut unfavorable_water_indices: Vec<usize> = vec![];

        let mut df = BoxData {
            connections: vec![],
            shells: vec![],
            kdtree_relations: vec![],
        };

        let mut shells = vec![];

        if let Some(water_box) = &self.water_box {
            let shell_id = water_box.number_of_shells();
            let add_noise = true;

            if let Some(connections) = &connections {
                if opt_disordered {
                    self.optimize_disordered_waters(&waters, connections);
                }
            }

            // The placement order is based on the best energy around each hydrogen anchor point
            // But if it returns [], it means the most favorable spots for placing water molecules
            // are definitively not that favorable, likely they are all outside the box.

            let water_orders = self.optimize_placement_order_grid(&waters, Some(1.));

            // And now we sample the position of all the water individually. The
            // acceptance of the new position/orientation is based on the metropolis
            // acceptance-rejection criteria. But per default, everything is accepted
            // with a favorable energy (< 0 kcal/mol).
            for i in &water_orders {
                let mut water = waters[*i].clone();
                let energy_position = self.optimize_position_grid(&mut water, add_noise, Some(1.));

                // The first great filter
                if utils::boltzmann_acceptance_rejection(
                    &vec![energy_position],
                    &vec![self.energy_cutoff],
                    self.temperature,
                )[0]
                {
                    // Build the explicit water
                    water.build_explicit_water(&self.water_model);

                    // Optimize the orientation
                    let energy_orientation = self.optimize_orientation_grid(&mut water);

                    // The last great energy filter
                    if utils::boltzmann_acceptance_rejection(
                        &energy_orientation,
                        &vec![self.energy_cutoff],
                        self.temperature,
                    )[0]
                    {
                        shells.push(Shell {
                            shell_id: shell_id + 1,
                            energy_position: Some(energy_position),
                            energy_orientation,
                        });
                        // data.append((shell_id + 1, energy_position, energy_orientation))
                        self.update_maps(&water);
                    } else {
                        unfavorable_water_indices.push(*i);
                    }
                } else {
                    unfavorable_water_indices.push(*i);
                };
            }
            // Keep only the good waters
            waters = water_orders
                .iter()
                .filter(|&i| !unfavorable_water_indices.contains(i))
                .map(|&i| waters[i].clone())
                .collect();

            // Keep connections of the good waters
            if let Some(connections) = connections {
                let mut connections_vec = connections.to_vec();
                connections_vec.retain(|connection| {
                    !unfavorable_water_indices.contains(&connection.molecule_j)
                });
                for (i, connection) in connections.iter_mut().enumerate() {
                    connection.molecule_j = i;
                }
                df.connections = connections_vec;
            }

            // Add water shell informations
            df.shells = shells;
        }

        (waters, df)
    }

    fn neighbor_points_grid(
        &self,
        water: &Water,
        from_edges: Option<f32>,
    ) -> (Vec<Vec3<f32>>, Vec<f32>) {
        // ) -> Vec<(Vec3<f32>, f32)> { // TODO port: review: this might be more intuitive
        let oxygen_type = water.atom_types(Some(&vec![0]));
        assert_eq!(1, oxygen_type.len());
        let oxygen_type = oxygen_type.first().unwrap();

        // This is how we select the allowed positions:
        // 1. Get all the point coordinates on the grid around the anchor (sphere). If the anchor type
        // is donor, we have to reduce the radius by 1 angstrom. Because the hydrogen atom is closer
        // to the water molecule than the heavy atom.
        // 2. Compute angles between all the coordinates and the anchor
        // 3. Select coordinates with an angle superior or equal to the choosen angle
        // 4. Get their energy

        let mut coord_sphere = if water.hb_type == "donor" {
            self.ad_map.neighbor_points(
                &water.hb_anchor,
                self.max_distance - 1.,
                self.min_distance - 1.,
            )
        } else {
            self.ad_map
                .neighbor_points(&water.hb_anchor, self.max_distance, self.min_distance)
        };

        if let Some(from_edges) = from_edges {
            let is_close = self.ad_map.is_close_to_edge(&coord_sphere, from_edges);
            coord_sphere = coord_sphere
                .into_iter()
                .zip(is_close.into_iter())
                .filter_map(|(coord, close)| if !close { Some(coord) } else { None })
                .collect();
        }

        // It is mandatory to normalize the hb_vector, otherwise you don't get the angle right
        let hb_vector = water.hb_anchor + (water.hb_vector + water.hb_anchor).normalized();

        coord_sphere = coord_sphere
            .into_iter()
            .filter(|coord| {
                let angle = utils::get_angle(coord, &water.hb_anchor, &hb_vector, true);
                angle >= self.angle
            })
            .collect();

        let energy_sphere = self
            .ad_map
            .energy_coordinates(&coord_sphere, &oxygen_type, "linear");

        (coord_sphere, energy_sphere)
    }

    fn optimize_placement_order_grid(
        &mut self,
        waters: &[Water],
        from_edges: Option<f32>,
    ) -> Vec<usize> {
        let mut energies = vec![];

        for water in waters {
            let (_, mut energy_sphere) = self.neighbor_points_grid(water, from_edges);

            if !energy_sphere.is_empty() {
                for es in energy_sphere.iter_mut() {
                    if *es == 0. {
                        *es = f32::INFINITY;
                    }
                }
                if let Some(min_value) = energy_sphere
                    .into_iter()
                    .reduce(|a, b| if a < b { a } else { b })
                {
                    energies.push(min_value)
                }
            } else {
                energies.push(f32::INFINITY)
            }
        }

        // Pick order based on Boltzmann choices
        let order: Vec<usize> =
            utils::boltzmann_choices(&energies, self.temperature, Some(energies.len()));

        let new_energies = order.iter().map(|&i| energies[i]).collect::<Vec<f32>>();

        if !order.is_empty() {
            let decisions = utils::boltzmann_acceptance_rejection(
                &new_energies,
                &vec![self.energy_cutoff],
                self.temperature,
            );

            return order
                .iter()
                .zip(decisions.iter()) // Pair each index with its corresponding boolean
                .filter_map(|(&idx, &keep)| if keep { Some(idx) } else { None }) // Keep only where `decisions` is true
                .collect();
        } else {
            return vec![];
        }
    }

    /// Optimize the position of the spherical water molecule.
    /// The movement of the water is contrained by the distance and
    /// the angle with the anchor.
    fn optimize_position_grid(
        &self,
        water: &mut Water,
        add_noise: bool,
        from_edges: Option<f32>,
    ) -> f32 {
        let oxygen_type = water.atom_types(Some(&vec![0]));
        assert_eq!(1, oxygen_type.len());
        let oxygen_type = oxygen_type.first().unwrap();

        let (coord_sphere, mut energy_sphere) = self.neighbor_points_grid(water, from_edges);

        if !energy_sphere.is_empty() {
            // Pick position based on Boltzmann choices

            let idx = utils::boltzmann_choices(&energy_sphere, self.temperature, None);

            if !idx.is_empty() {
                let mut new_coord = coord_sphere[idx[0]];

                if add_noise {
                    let limit = self.ad_map.spacing / 2.;
                    let mut rng = rand::thread_rng();
                    let noise = Vec3::new(
                        rng.gen_range(-limit..=limit),
                        rng.gen_range(-limit..=limit),
                        rng.gen_range(-limit..=limit),
                    );
                    new_coord += noise;
                }

                // Update the coordinates
                water.translate(&utils::vector(&water.coordinates[1], &new_coord));
                return energy_sphere[idx[0]];
            }
        }

        // If we do not find anything, at least we return the energy
        // of the current water molecule.
        // TODO port: energy_coordinates first arg likely incorrect
        self.ad_map
            .energy_coordinates(&vec![water.coordinates[1]], &oxygen_type, "linear")[0]
    }

    fn optimize_disordered_waters(&mut self, waters: &[Water], connections: &[Bond]) -> Vec<f32> {
        let mut disordered_energies: Vec<f32> = vec![];
        self.rotation = 10.0;

        // Number of rotation necessary to do a full spin
        let n_rotation = (360.0 / self.rotation).floor() as i32 - 1;
        let rotation = self.rotation.to_radians();

        if let Some(receptor) = &self.receptor {
            for row in &receptor.rotatable_bonds {
                let mut energies = vec![];
                let mut angles = vec![];
                let mut rot_waters = vec![];

                // Get index of all the waters attached
                // to a disordered group by looking at the connections
                let molecule_j = filter_and_extract_molecule_j(connections, &row);
                rot_waters.extend(molecule_j.iter().map(|&j| waters[j].clone()));

                if !rot_waters.is_empty() {
                    // Get energy of the favorable disordered waters
                    // TODO port: pass additional parameters
                    // energy_waters = np.array([ad_map.energy(w.atom_informations(), ignore_electrostatic=True, ignore_desolvation=True) for w in rot_waters])
                    // TODO port other parameters for energy()
                    let energy_waters_nested = rot_waters
                        .iter()
                        .map(|w| {
                            self.ad_map
                                .energy(&w.atoms(), vec![], false, true, true, "linear", true)
                                .flatten()
                        })
                        .collect::<Vec<Vec<f32>>>();
                    let energy_waters: Vec<f32> = energy_waters_nested
                        .into_iter()
                        .flat_map(|v| v.into_iter())
                        .collect();
                    for mut w in &energy_waters {
                        w = &w.min(0.);
                    }

                    energies.push(energy_waters.into_iter().sum());
                    // Current angle of the disordered group
                    let mut current_angle = utils::dihedral(row.dihedral_pars(), false);
                    angles.push(current_angle);

                    // Find all the atoms that depends on these atoms. This
                    // will be useful when we will want to rotate a whole sidechain.
                    // Atom children has to be initialized before
                    // molecule._OBMol.FindChildren(atom_children, match[2], match[3])
                    // print np.array(atom_children)

                    // Atoms 3 and 2 define the rotation axis
                    let p1 = row.atom_k_xyz;
                    let p2 = row.atom_j_xyz;

                    // Scan all the angles
                    for i in 0..n_rotation {
                        // TODO: Performance wise, we should not update water
                        // coordinates everytime. Coordinates should be extracted
                        // before doing the optimization and only at the end
                        // we update the coordinates of the water molecules.
                        for rot_water in rot_waters.iter_mut() {
                            let p0 = rot_water.coordinates(Some(&[1]))[0];
                            let p_new = utils::rotate_point(p0, p1, p2, rotation);
                            rot_water.update_coordinates(&p_new, 1);
                        }

                        // Get energy and update the current angle (increment rotation)
                        // TODO port: same block of code as bit above, maybe refactor
                        let energy_waters_nested = rot_waters
                            .iter()
                            .map(|w| {
                                self.ad_map
                                    .energy(&w.atoms(), vec![], false, true, true, "linear", true)
                                    .flatten()
                            })
                            .collect::<Vec<Vec<f32>>>();
                        let mut energy_waters: Vec<f32> = energy_waters_nested
                            .into_iter()
                            .flat_map(|v| v.into_iter())
                            .collect();
                        for mut w in &energy_waters {
                            w = &w.min(0.);
                        }
                        energies.push(energy_waters.into_iter().sum());

                        current_angle += rotation;
                        angles.push(current_angle);
                    }

                    // Pick state based on Boltzmann choices
                    let i = utils::boltzmann_choices(&energies, self.temperature, None)[0];

                    disordered_energies.push(energies[i]);

                    // Calculate the best angle, based on how much we rotated

                    let best_angle =
                        ((360. - current_angle.to_degrees()) + angles[i].to_degrees()).to_radians();
                    // Update coordinates to the choosen state
                    for rot_water in rot_waters.iter_mut() {
                        // TODO port: what's up with indices here
                        // p0 = rot_water.coordinates(1)[0]
                        let p0 = rot_water.coordinates(Some(&[1]))[0];
                        let p_new = utils::rotate_point(p0, p1, p2, best_angle);
                        rot_water.update_coordinates(&p_new, 1);
                        // Update also the anchor point
                        rot_water.hb_anchor =
                            utils::rotate_point(rot_water.hb_anchor, p1, p2, best_angle);
                        rot_water.hb_vector =
                            utils::rotate_point(rot_water.hb_vector, p1, p2, best_angle);
                    }
                }
            }
        }

        disordered_energies
    }

    /// Optimize the orientation of the TIP5P water molecule using the grid.
    /// TODO port: is return type correct
    fn optimize_orientation_grid(&self, water: &mut Water) -> Vec<f32> {
        let oxygen_xyz = &water.coordinates;
        let water_info = water.atom_informations(None);
        // oxygen_xyz = water.coordinates(1)
        // water_info = water.atom_informations()
        // energies = np.zeros(self._water_orientations.shape[0])
        let num_orientations = self.water_orientations.len();
        let energies: Vec<f32> = vec![0.0; num_orientations];

        // Translate the coordinates
        // let water_orientations = self.water_orientations + oxygen_xyz;
        // TODO port: double check, for sure this is verbose
        let water_orientations: Vec<Vec3<f32>> = self
            .water_orientations
            .iter()
            .zip(oxygen_xyz.iter())
            .map(|(orientation, offset)| {
                orientation
                    .iter()
                    .zip([offset.x as f32, offset.y as f32, offset.z as f32].iter())
                    .map(|(value, offset_value)| value + offset_value)
                    .collect()
            })
            .collect();

        // Get the energies for each atom
        // Oxygen first
        // energies += ad_map.energy_coordinates(oxygen_xyz, water_info["t"][0])
        // TODO port: we shouldn't assume a first element
        let atom_type = &water_info[0].t;
        self.ad_map
            .energy_coordinates(oxygen_xyz, &atom_type, "linear");
        //... and then hydrogens/lone-pairs
        // for i, atom_type in enumerate(water_info["t"][1:]):
        for (i, atom) in water_info.iter().enumerate() {
            // TODO port: first parameter is most likely incorrect (not using i) but it's what compiles currently
            // energies += ad_map.energy_coordinates(water_orientations[:,i], atom_type)
            self.ad_map
                .energy_coordinates(&water_orientations, &atom.t, "linear");
        }

        // Pick orientation based on Boltzmann choices
        // idx = utils.boltzmann_choices(energies, self._temperature)

        let idx = utils::boltzmann_choices(&energies, self.temperature, None);

        if idx.len() > 0 {
            // Update the coordinates with the selected orientation, except oxygen (i + 2)
            // TODO port: again we've to omit an index to compile, review
            // let new_orientation = water_orientations[idx[0]];
            let new_orientation = water_orientations;
            for (i, xyz) in new_orientation.iter().enumerate() {
                water.update_coordinates(&xyz, i + 2);
            }

            // TODO port: and again omitting an index..
            // return energies[idx[0]];
            return energies;
        }

        // If we do not find anything, at least we return the energy
        // of the current water molecule.
        let current_energy =
            self.ad_map
                .energy(&water.atoms(), vec![], false, true, true, "linear", true);

        current_energy.flatten()
    }

    /// If we choose the closest point in the grid and not the coordinates of the
    /// oxygen as the center of the grid, it is because we want to avoid any edge effect
    /// when we will combine the small box to the bigger box, and also the energy is
    /// not interpolated but it is coming from the grid directly.
    fn update_maps(&mut self, water: &Water) {
        let box_size = Vec3 {
            x: 8.,
            y: 8.,
            z: 8.,
        };

        let npts = ((box_size / self.ad_map.spacing).round() / 2.).floor() * 2. + 1.;
        let water_xyz = &water.coordinates;

        // The center is the closest grid point from the water molecule
        // TODO port: review &water.coordinates[0], should not assume 1 element
        let center_xyz =
            self.ad_map
                .neighbor_points(&water.coordinates[0], self.ad_map.spacing, 0.);
        // TODO port: added a [0] here just to make it compile, most likely more adjustments needed
        let center_index = self.ad_map.cartesian_to_index(&center_xyz[0]);
        let center_index_ft = center_index.map(|i| i as f32);

        // Get grid indices in the receptor box
        // Necessary, in order to add the water map to the receptor map
        let size_half_box = (npts - Vec3::new(1., 1., 1.)) / 2.;

        let i_min = (center_index_ft - size_half_box).map2(npts, |value, max| value.clamp(0., max));
        let i_max = (center_index_ft + size_half_box + Vec3::new(1., 1., 1.))
            .map2(npts, |value, max| value.clamp(0., max));

        let x_range = i_min.x..i_max.x;
        let y_range = i_min.y..i_max.y;
        let z_range = i_min.z..i_max.z;

        // Return the ranges as a tuple
        let indices = (x_range, y_range, z_range);

        // Convert the indices in xyz coordinates, and translate it in order
        // that the water molecule is at the origin
        // TODO port: confirm, is there no more elegant way to do this
        let step = 1.0;
        let x = generate_float_range(i_min.x, i_max.x, step);
        let y = generate_float_range(i_min.y, i_max.y, step);
        let z = generate_float_range(i_min.z, i_max.z, step);

        let (X, Y, Z) = meshgrid(&x, &y, &z);

        let x_flat = X.iter().cloned().collect::<Vec<f32>>();
        let y_flat = Y.iter().cloned().collect::<Vec<f32>>();
        let z_flat = Z.iter().cloned().collect::<Vec<f32>>();

        let num_points = x_flat.len();
        let mut grid_index_array2 = Array2::<f32>::zeros((num_points, 3));

        for (i, ((&x, &y), &z)) in x_flat.iter().zip(&y_flat).zip(&z_flat).enumerate() {
            grid_index_array2[[i, 0]] = x;
            grid_index_array2[[i, 1]] = y;
            grid_index_array2[[i, 2]] = z;
        }

        let grid_index: Vec<Vec3<i32>> = convert_array2_to_vec3(grid_index_array2);

        let mut grid_xyz = self.ad_map.index_to_cartesian(grid_index);
        grid_xyz.iter_mut().for_each(|point| *point -= water_xyz[0]);

        // Get the quaternion between the current water and the water reference
        // The water is translating to the origin before
        let origin = water_xyz[0];
        let translated_water = water_xyz
            .iter()
            .map(|&point| point - origin)
            .collect::<Vec<Vec3<f32>>>();

        let q = utils::quaternion_rotate(&translated_water, &self.water_ref.coordinates(None));

        // Define the vector to rotate
        let vec = Vector3::new(0.0, 1.0, 0.0);
        // Define the quaternion for rotation (90 degrees around X-axis)
        let rotation = UnitQuaternion::from_axis_angle(&Vector3::x_axis(), 90.0_f32.to_radians());
        // Rotate the vector using the quaternion
        let rotated_vec = rotation * vec;
        let vek_vec = Vec3::new(rotated_vec.x, rotated_vec.y, rotated_vec.z);

        let quat = nalgebra::Quaternion::new(q.w, q.x, q.y, q.z);
        let unit_quat = UnitQuaternion::from_quaternion(quat);

        // Rotate the grid!
        let rotated_grid_xyz: Vec<Vec3<f32>> = grid_xyz
            .iter()
            .map(|point| {
                // TODO port: avoid these transformations back and forth nalgebra and vek, use only nalgebra maybe
                let vec = Vector3::new(point.x, point.y, point.z);
                let rotated = rotation * vec;
                Vec3::new(rotated.x, rotated.y, rotated.z)
            })
            .collect();

        let mut atom_types = vec!["SW"];
        if self.water_model == "tip3p" {
            atom_types.extend(vec!["OW", "HW"]);
        } else if self.water_model == "tip5p" {
            atom_types.extend(vec!["OT", "HT", "LP"]);
        }

        // Interpolation baby!
        for atom_type in atom_types {
            let mut energy =
                self.water_map
                    .energy_coordinates(&rotated_grid_xyz, atom_type, "linear");
            let reshaped_energy = swap_and_reshape(energy, &x, &y, &z);
            // TODO port
            // v += energy;

            // TODO port what does this do?
            // self.ad_map.maps_interpn[atom_type] = self
            //     .ad_map
            //     .generate_affinity_map_interpn(&self.ad_map.maps[atom_type]);
        }
    }
}

// TODO port: return type unlikely to be correct, also, implement
// this ports this line:
// energy = np.swapaxes(energy.reshape((y.shape[0], x.shape[0], z.shape[0])), 0, 1)
fn swap_and_reshape(energy: Vec<f32>, x: &Vec<f32>, y: &Vec<f32>, z: &Vec<f32>) -> Vec<f32> {
    todo!()
}

fn convert_array2_to_vec3(array: Array2<f32>) -> Vec<Vec3<i32>> {
    array
        .rows()
        .into_iter()
        .map(|row| Vec3::new(row[0] as i32, row[1] as i32, row[2] as i32))
        .collect()
}

// TODO port: review
pub fn meshgrid(x: &[f32], y: &[f32], z: &[f32]) -> (Array3<f32>, Array3<f32>, Array3<f32>) {
    let nx = x.len();
    let ny = y.len();
    let nz = z.len();

    let mut x_grid = Array3::<f32>::zeros((nx, ny, nz));
    let mut y_grid = Array3::<f32>::zeros((nx, ny, nz));
    let mut z_grid = Array3::<f32>::zeros((nx, ny, nz));

    for (i, &xi) in x.iter().enumerate() {
        x_grid.index_axis_mut(Axis(0), i).fill(xi);
    }

    for (j, &yj) in y.iter().enumerate() {
        y_grid.index_axis_mut(Axis(1), j).fill(yj);
    }

    for (k, &zk) in z.iter().enumerate() {
        z_grid.index_axis_mut(Axis(2), k).fill(zk);
    }

    (x_grid, y_grid, z_grid)
}

fn flatten_and_stack(x_grid: Array3<f32>, y_grid: Array3<f32>, z_grid: Array3<f32>) -> Array2<f32> {
    // Flatten the 3D grids into 1D arrays
    let x_flat = x_grid.iter().cloned().collect::<Vec<f32>>();
    let y_flat = y_grid.iter().cloned().collect::<Vec<f32>>();
    let z_flat = z_grid.iter().cloned().collect::<Vec<f32>>();

    // Stack the flattened arrays into a single 2D array
    let num_points = x_flat.len();
    let mut grid_index = Array2::<f32>::zeros((num_points, 3));

    for (i, ((&x, &y), &z)) in x_flat.iter().zip(&y_flat).zip(&z_flat).enumerate() {
        grid_index[[i, 0]] = x;
        grid_index[[i, 1]] = y;
        grid_index[[i, 2]] = z;
    }

    grid_index
}

fn generate_float_range(start: f32, end: f32, step: f32) -> Vec<f32> {
    let mut values = Vec::new();
    let mut current = start;

    while current < end {
        values.push(current);
        current += step;
    }

    values
}

fn filter_and_extract_molecule_j<'a>(connections: &'a [Bond], bond: &Bond) -> Vec<usize> {
    connections
        .iter()
        .filter(|&conn| conn.atom_i == bond.atom_i || Some(conn.atom_i) == bond.atom_j)
        .map(|b| b.molecule_j) // Extract molecule_j
        .collect()
}
