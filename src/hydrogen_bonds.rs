use crate::force_fields_parser::HbAtom;
use std::collections::HashMap;
use vek::Vec3;

#[derive(Debug, Clone)]
pub struct HydrogenBond {
    pub atom_i: usize,
    pub vector_xyz: Vec3<f32>,
    pub anchor_type: String,
    pub atom_types: HashMap<String, HbAtom>,
}
