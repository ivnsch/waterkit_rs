use std::collections::{HashMap, HashSet};

use crate::{
    atom::Atom,
    fld_parser::parse_fld,
    map_parser::{parse_map, parse_map_infos, GridInformation},
};
use anyhow::{anyhow, Result};
use kd_tree::KdTree;
use ndarray::Array3;
use vek::Vec3;

#[derive(Debug, Clone)]
pub struct Map {
    pub kd_tree: KdTree<[f32; 3]>,
    pub spacing: f32,
    pub xmin: f32,
    pub xmax: f32,
    pub ymin: f32,
    pub ymax: f32,
    pub zmin: f32,
    pub zmax: f32,
    pub maps: HashMap<String, Array3<f32>>,
    pub maps_interpn: HashMap<String, Interpolator>,
    pub npts: Vec3<f32>,
    pub map_files: Vec<String>,
    pub labels: Vec<String>,
    pub center: Vec3<f32>,
    pub edges: (Vec<f32>, Vec<f32>, Vec<f32>),
}

#[derive(Debug, Clone)]
pub struct Interpolator {
    affinity_map: Array3<f32>,
}

impl Interpolator {
    // Note assuming always trilinear. Are other methods actually used in original library?
    fn interpolate(&self, points: &[Vec3<f32>], method: &str) -> Vec<f32> {
        let mut results = Vec::with_capacity(points.len());

        // println!("affinity map: {:?}", self.affinity_map);

        for point in points {
            results.push(interpolate_point(&self.affinity_map, point));
        }

        results
    }
}

fn interpolate_point(affinity_map: &Array3<f32>, point: &Vec3<f32>) -> f32 {
    let (x_len, y_len, z_len) = affinity_map.dim();

    let x = point.x;
    let y = point.y;
    let z = point.z;

    // Compute the integer coordinates of the cube surrounding (x, y, z)
    let x0 = x.floor() as isize;
    let y0 = y.floor() as isize;
    let z0 = z.floor() as isize;

    let x1 = (x0 + 1).min(x_len as isize - 1);
    let y1 = (y0 + 1).min(y_len as isize - 1);
    let z1 = (z0 + 1).min(z_len as isize - 1);

    // Extrapolation: Clamp values within the grid boundaries
    // TODO port: better extrapolation method? what does python grid interpolator do here?
    let x0 = x0.clamp(0, x_len as isize - 1) as usize;
    let y0 = y0.clamp(0, y_len as isize - 1) as usize;
    let z0 = z0.clamp(0, z_len as isize - 1) as usize;
    let x1 = x1.clamp(0, x_len as isize - 1) as usize;
    let y1 = y1.clamp(0, y_len as isize - 1) as usize;
    let z1 = z1.clamp(0, z_len as isize - 1) as usize;

    // Perform trilinear interpolation
    // https://en.wikipedia.org/wiki/Trilinear_interpolation

    // Retrieve values from the 8 corners of the cube
    let c000 = affinity_map[[x0, y0, z0]];
    let c001 = affinity_map[[x0, y0, z1]];
    let c010 = affinity_map[[x0, y1, z0]];
    let c011 = affinity_map[[x0, y1, z1]];
    let c100 = affinity_map[[x1, y0, z0]];
    let c101 = affinity_map[[x1, y0, z1]];
    let c110 = affinity_map[[x1, y1, z0]];
    let c111 = affinity_map[[x1, y1, z1]];

    // Compute interpolation weights
    // NOTE leaving out denominators as our grid spacing is 1, which makes them 1
    // NOTE also, we clamp because we clamped the cube's vertices to grid boundaries
    // TODO port: last note points to issue - we shouldn't really clamp here,
    // it's ok if far away point generates large values
    // but using this later to compute probability weights causes errors so leaving it like this for now
    let xd = (x - x0 as f32).clamp(0.0, 1.0);
    let yd = (y - y0 as f32).clamp(0.0, 1.0);
    let zd = (z - z0 as f32).clamp(0.0, 1.0);

    let c00 = c000 * (1.0 - xd) + c100 * xd;
    let c01 = c001 * (1.0 - xd) + c101 * xd;
    let c10 = c010 * (1.0 - xd) + c110 * xd;
    let c11 = c011 * (1.0 - xd) + c111 * xd;

    let c0 = c00 * (1.0 - yd) + c01 * yd;
    let c1 = c10 * (1.0 - yd) + c11 * yd;

    let interpolated_value = c0 * (1.0 - zd) + c1 * zd;
    interpolated_value
}

impl Map {
    pub async fn new(map_files: Vec<String>, labels: Vec<String>) -> Result<Map> {
        if labels.len() != map_files.len() {
            return Err(anyhow!(
                "Map files and labels must have the same length, labels: {}, map_files: {}",
                labels.len(),
                map_files.len()
            ));
        }

        let mut prev_grid_information = None;

        let mut maps: HashMap<String, String> = HashMap::new();

        // Get information (center, spacing, nelements) from each grid
        // and make sure they are all identical

        for (label, map_file) in labels.iter().zip(map_files.iter()) {
            let grid_information = grid_information_from_map(&map_file).await?;

            if let Some(prev_grid_information) = prev_grid_information {
                if prev_grid_information != grid_information {
                    return Err(anyhow!(
                        "grid {} is different from the previous one.",
                        label
                    ));
                }
            }
            prev_grid_information = Some(grid_information);
            maps.insert(label.to_string(), map_file.to_string());
        }

        // we require that all grid informations are the same, so we can use any.
        let grid_information =
            prev_grid_information.ok_or_else(|| anyhow!("No grid informations, exit"))?;

        // Compute min and max coordinates
        // Half of each side
        let l = ((grid_information.nelements - 1.) * grid_information.spacing) / 2.;
        // Minimum and maximum coordinates
        let min = grid_information.center - l;
        let max = grid_information.center + l;

        let mut affinity_maps: HashMap<String, Array3<f32>> = HashMap::new();
        let mut maps_interpn: HashMap<String, Interpolator> = HashMap::new();

        // Read all the affinity maps
        for (label, map_file) in maps {
            let affinity_map = read_affinity_map(&map_file).await?;
            maps_interpn.insert(label.clone(), generate_affinity_map_interpn(&affinity_map));
            affinity_maps.insert(label, affinity_map);
        }

        let edges_x = linspace(min.x, max.x, grid_information.nelements.x as usize);
        let edges_y = linspace(min.y, max.y, grid_information.nelements.y as usize);
        let edges_z = linspace(min.z, max.z, grid_information.nelements.z as usize);

        Ok(Map {
            kd_tree: build_kdtree_from_grid(&edges_x, &edges_y, &edges_z),
            spacing: grid_information.spacing,
            npts: grid_information.nelements,
            center: grid_information.center,
            map_files,
            labels,
            xmin: min.x,
            ymin: min.y,
            zmin: min.z,
            xmax: max.x,
            ymax: max.y,
            zmax: max.z,
            maps: affinity_maps,
            maps_interpn: maps_interpn,
            edges: (edges_x, edges_y, edges_z),
        })
    }

    /// Check if coordinates are in the map.
    /// Args:
    ///     xyz (array_like): 2d array of the 3d coordinates
    /// Returns:
    ///     ndarray: 1d Numpy array of boolean
    /// TODO port: assuming that xyz is just a Vec3<f32> if it's really 2d needs update
    /// TODO port: if it's correct like this, update docs
    pub fn is_in_map(&self, xyz: Vec3<f32>) -> bool {
        let x = xyz.x;
        let y = xyz.y;
        let z = xyz.z;

        let x_in = self.xmin <= x && x <= self.xmax;
        let y_in = self.ymin <= y && y <= self.ymax;
        let z_in = self.zmin <= z && z <= self.zmax;

        let all_in = x_in && y_in && z_in;

        all_in
    }

    /// Get energy interaction of a molecule based of the grid.
    /// Args:
    ///     nd (ndarray): Structure numpy array with columns (i, xyz, q, t)
    ///     ignore_atom_types (list): list of atom types/terms to ignore (default: None)
    ///     ignore_electrostatic (bool): to ignore electrostatic term (default: False)
    ///     ignore_desolvation (bool): to ignore desolvation term (default: False)
    ///     method (str): Interpolate method (default: linear)
    /// Returns:
    ///     float: Grid energy interaction
    pub fn energy(
        &self,
        atoms: &[Atom],
        ignore_atom_types: Vec<String>,
        ignore_vdw_hb: bool,
        ignore_electrostatics: bool,
        ignore_desolvation: bool,
        method: &str,
    ) -> f32 {
        let mut energies = vec![];
        let mut vdw_hb = vec![];
        let elec = 0.;
        let desolv = 0.;

        let ignore_atom_types: HashSet<_> = ignore_atom_types.into_iter().collect();

        let mut atom_types: HashSet<_> = atoms.iter().map(|a| a.t.clone()).collect();
        atom_types = atom_types.difference(&ignore_atom_types).cloned().collect();

        let mut index: Vec<usize> = vec![];

        for atom_type in atom_types {
            let xyz = if atoms.len() > 1 {
                Some(
                    atoms
                        .into_iter()
                        .enumerate()
                        .filter_map(|(i, a)| {
                            if a.t == atom_type {
                                index.push(i);
                                Some(a.coords)
                            } else {
                                None
                            }
                        })
                        .collect(),
                )
            } else {
                atoms.first().map(|a| vec![a.coords])
            };

            if let Some(xyz) = xyz {
                if !ignore_vdw_hb {
                    // vdw_hb = self._maps_interpn[atom_type](xyz, method=method)
                    vdw_hb = self.maps_interpn[&atom_type].interpolate(&xyz, method);
                    energies.push(vdw_hb);
                }

                if !ignore_electrostatics {
                    let q: Vec<f32> = if atoms.len() > 1 {
                        index
                            .iter()
                            .map(|&i| atoms[i].clone())
                            .map(|a| a.q)
                            .collect()
                    } else {
                        atoms.iter().map(|a| a.q).collect()
                    };

                    let interp = self.maps_interpn["Electrostatics"].interpolate(&xyz, method);
                    let elec: Vec<f32> = q
                        .iter()
                        .zip(interp.iter())
                        .map(|(q_val, interp_val)| q_val * interp_val)
                        .collect();
                    energies.push(elec);
                }

                if !ignore_desolvation {
                    let desolv = self.maps_interpn["Desolvation"].interpolate(&xyz, method);
                    energies.push(desolv);
                }
            }
        }

        let flat = energies.into_iter().flatten();
        flat.sum()
    }

    /// Grid coordinates around a point at a certain distance.
    /// Args:
    ///     xyz (array_like): 3d coordinate of a point
    ///     radius (float): max radius
    ///     min_radius (float): min radius (default: 0)
    /// Returns:
    ///     ndarray: 2d Numpy array of coordinates
    pub fn neighbor_points(&self, xyz: &Vec3<f32>, radius: f32, min_radius: f32) -> Vec<Vec3<f32>> {
        // TODO port: HACK! had to increase radius because otherwise never finds any points
        // data etc. seems to be same as python so don't know what's issue? need to investigate more
        let radius = radius * 3.;
        let arr = self.kd_tree.within_radius(&xyz.into_array(), radius);
        // println!(
        //     "kd_tree within radius: {}, xyz: {:?}, arr len: {:?}, flat: {:?} ",
        //     radius,
        //     xyz,
        //     arr.len(),
        //     self.kd_tree.len()
        // );
        let mut coordinates: Vec<Vec3<f32>> = arr.into_iter().map(|r| Vec3::from(*r)).collect();

        if min_radius > 0. {
            coordinates.retain(|coord| (*xyz - *coord).magnitude() >= min_radius);
        }

        coordinates
    }

    /// Check if the points xyz is at a certain distance of the edge of the box.
    /// Args:
    ///     xyz (array_like): 2d array of the 3d coordinates
    ///     distance (float): distance
    /// Returns:
    ///     ndarray: 1d Numpy array of boolean
    pub fn is_close_to_edge(&self, xyz: &[Vec3<f32>], distance: f32) -> Vec<bool> {
        xyz.into_iter()
            .map(|xyz| self.is_close_to_edge_point(xyz, distance))
            .collect()
    }

    pub fn is_close_to_edge_point(&self, xyz: &Vec3<f32>, distance: f32) -> bool {
        let x = xyz.x;
        let y = xyz.y;
        let z = xyz.z;

        let x_close = (self.xmin - x).abs() <= distance || (self.xmax - x).abs() <= distance;
        let y_close = (self.ymin - y).abs() <= distance || (self.ymax - y).abs() <= distance;
        let z_close = (self.zmin - z).abs() <= distance || (self.zmax - z).abs() <= distance;
        let close_to = x_close || y_close || z_close;

        close_to
    }

    /// Grid energy of each coordinates xyz.
    /// Args:
    ///     xyz (array_like): 2d array of 3d coordinates
    ///     atom_type (str): name of the atom type
    ///     method (str): Interpolate method (default: linear)
    /// Returns:
    ///     ndarray: 1d Numpy array of the energy values
    pub fn energy_coordinates(&self, xyz: &[Vec3<f32>], atom_type: &str, method: &str) -> Vec<f32> {
        self.maps_interpn[atom_type].interpolate(xyz, method)
    }

    // Return the closest grid index of the cartesian grid coordinates
    pub fn cartesian_to_index(&self, xyz: &Vec3<f32>) -> Vec3<i32> {
        let (mins, _) = calculate_bounds(&self.kd_tree);
        let idx: Vec3<i32> = ((xyz - Vec3::from(mins)) / self.spacing).map(|v| v.round() as i32);

        // All the index values outside the grid are clipped (limited) to the nearest index
        idx.map2(self.npts, |i, n| i.clamp(0, n as i32))
    }

    /// Return the cartesian grid coordinates associated to the grid index
    pub fn index_to_cartesian(&self, idx: Vec<Vec3<i32>>) -> Vec<Vec3<f32>> {
        idx.into_iter()
            .map(|i| self.index_to_cartesian_internal(i))
            .collect()
    }

    fn index_to_cartesian_internal(&self, idx: Vec3<i32>) -> Vec3<f32> {
        let (mins, _) = calculate_bounds(&self.kd_tree);
        let mins = Vec3::from(mins);
        mins + idx.map(|i| i as f32) * self.spacing
    }

    /// Return a interpolate function from the grid and the affinity map.
    /// This helps to interpolate the energy of coordinates off the grid.
    // pub fn generate_affinity_map_interpn(&self, affinity_map: &[Vec<f32>]) -> Interpolator {
    pub fn generate_affinity_map_interpn(&self, affinity_map: &Array3<f32>) -> Interpolator {
        // TODO regular grid interpolator
        Interpolator {
            affinity_map: affinity_map.clone(),
        }
    }

    /// Read a fld file.
    /// The AutoDock map files are read using the information contained
    /// into the fld file. The fld file is created by AutoGrid.
    /// Args:
    ///     fld_file (str): pathname of the AutoGrid fld file.
    /// Returns:
    ///     Map: Instance of Map object.
    pub async fn from_fld(fld_file: &str) -> Result<Map> {
        let parsed_fld = parse_fld(fld_file).await?;

        Map::new(parsed_fld.map_files, parsed_fld.labels).await
    }
}

fn calculate_bounds(tree: &KdTree<[f32; 3]>) -> ([f32; 3], [f32; 3]) {
    let mut mins = [f32::MAX; 3];
    let mut maxs = [f32::MIN; 3];

    for &point in tree.iter() {
        for i in 0..3 {
            mins[i] = mins[i].min(point[i]);
            maxs[i] = maxs[i].max(point[i]);
        }
    }

    (mins, maxs)
}

/// Read grid information in the map file
async fn grid_information_from_map(path: &str) -> Result<GridInformation> {
    parse_map_infos(path).await
}

/// Take a grid file and extract gridcenter, spacing and npts info
async fn read_affinity_map(path: &str) -> Result<Array3<f32>> {
    let map = parse_map(path).await?;
    // Some sorceries happen here --> swap x and z axes
    let affinities = map.affinities;
    let npts = map.infos.nelements + 1.;
    let shape = (npts[2] as usize, npts[1] as usize, npts[0] as usize);
    let mut affinity_3d =
        Array3::from_shape_vec(shape, affinities).expect("Failed to reshape into 3D array");

    affinity_3d.swap_axes(0, 2);

    Ok(affinity_3d)
}

fn generate_affinity_map_interpn(affinity_map: &Array3<f32>) -> Interpolator {
    Interpolator {
        affinity_map: affinity_map.clone(),
    }
}

fn linspace(start: f32, stop: f32, num: usize) -> Vec<f32> {
    if num < 2 {
        return vec![start];
    }
    let step = (stop - start) / (num as f32 - 1.0);
    (0..num).map(|i| start + step * i as f32).collect()
}

fn build_kdtree_from_grid(edges_x: &[f32], edges_y: &[f32], edges_z: &[f32]) -> KdTree<[f32; 3]> {
    // empirically equivalent to python:
    // X, Y, Z = np.meshgrid(x, y, z)
    // xyz = np.stack((X.ravel(), Y.ravel(), Z.ravel()), axis=-1)
    let mut kd_tree_points: Vec<[f32; 3]> = Vec::new();
    for &x in edges_x {
        for &y in edges_y {
            for &z in edges_z {
                kd_tree_points.push([x, y, z]);
            }
        }
    }
    KdTree::build_by_ordered_float(kd_tree_points)
}
