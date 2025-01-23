use nalgebra::Matrix4;
use rand::{distributions::WeightedIndex, prelude::Distribution, thread_rng, Rng};
use vek::{Mat3, Vec3, Vec4};

use crate::atom::DihedralPars;

/// Dihedral angle.

/// Source:
/// https://stackoverflow.com/questions/20305272/dihedral-torsion-angle-from-four-points-in-cartesian-coordinates-in-python
pub fn dihedral(p: DihedralPars, degree: bool) -> f32 {
    let p0 = p.atom_i_xyz;
    let p1 = p.atom_j_xyz;
    let p2 = p.atom_k_xyz;
    let p3 = p.atom_l_xyz;

    let b0 = -1.0 * (p1 - p0);
    let mut b1 = p2 - p1;
    let b2 = p3 - p2;

    // normalize b1 so that it does not influence magnitude of vector
    // rejections that come next
    b1 /= b1.magnitude();

    // vector rejections
    // v = projection of b0 onto plane perpendicular to b1
    //   = b0 minus component that aligns with b1
    // w = projection of b2 onto plane perpendicular to b1
    //   = b2 minus component that aligns with b1
    let v = b0 - b0.dot(b1) * b1;
    let w = b2 - b2.dot(b1) * b1;

    // angle between v and w in a plane is the torsion angle
    // v and w may not be normalized but that is fine since tan is y/x
    let x = v.dot(w);
    let y = (b1.cross(v)).dot(w);

    let angle = y.atan2(x);

    if degree {
        angle.to_degrees()
    } else {
        angle
    }
}

/// Rotate the point p around the axis p1-p2
/// Source: http://paulbourke.net/geometry/rotate/PointRotate.py
/// TODO port: note that we're not setting "back" p here, just returning, check that impl doesn't expect p to be mutated
pub fn rotate_point(p: Vec3<f32>, p1: Vec3<f32>, p2: Vec3<f32>, angle: f32) -> Vec3<f32> {
    //  Translate the point we want to rotate to the origin
    let pn = p - p1;

    // Get the unit vector from the axis p1-p2
    let mut n = p2 - p1;
    n.normalize();

    // Setup the rotation matrix
    let c = angle.cos();
    let t = 1. - angle.cos();
    let s = angle.sin();
    let x = n[0];
    let y = n[1];
    let z = n[2];

    #[rustfmt::skip]
    let r = Mat3::new (
        (t*x).powi(2) + c, t*x*y - s*z, t*x*z + s*y,
        t*x*y + s*z, (t*y).powi(2) + c, t*y*z - s*x,
        t*x*z - s*y, t*y*z + s*x, (t*z).powi(2) + c
    );

    // ... and apply it
    let ptr = r * pn;

    ptr + p1
}

pub fn boltzmann_probabilities(energies: &[f32], temperature: f32) -> Vec<f32> {
    // Boltzmann constant (kcal/mol)
    let kb = 0.0019872041;

    // let d =     d = np.exp(-energies / (kb * temperature))
    // let d = (-energies / (kb * temperature)).exp();

    let d = energies
        .iter()
        .map(|&energy| (-energy / (kb * temperature)).exp());

    let d_sum: f32 = d.clone().sum();

    let p = if d_sum > 0. {
        d.map(|d| d / d_sum).collect()
    } else {
        vec![0.0; energies.len()]
    };
    p
}

/// Choose state i based on boltzmann probability.
pub fn boltzmann_choices(energies: &[f32], temperature: f32, size: Option<usize>) -> Vec<usize> {
    // TODO port: HACK! scaling down energies, we get values ~100x of python implementation. Fix!
    let energies = energies
        .to_vec()
        .into_iter()
        .map(|e| e / 400.)
        .collect::<Vec<f32>>();

    // println!("energies: {:?} temperature: {:?}", energies, temperature);
    let mut p = boltzmann_probabilities(&energies, temperature);
    // println!("p: {:?}", p);

    // TODO port HACK adding random numbers so it doesn't crash later accessing empty array
    if p.iter().sum::<f32>() == 0. {
        let mut rng = rand::thread_rng();
        p = p.iter().map(|_| rng.gen()).collect();
    }
    if p.iter().sum::<f32>() == 0. {
        vec![]
    } else {
        let mut size = size.unwrap_or(1);
        if size == 0 {
            size = 1;
        }
        // If some prob. in p are zero, ValueError: size of nonzero p is lower than size
        let non_zero: usize = p.iter().filter(|&&x| x != 0.0).count();

        let size = if non_zero < size { non_zero } else { size };

        // TODO port: write some tests to ensure this is the same
        // i = np.random.choice(len(energies), size, False, p)
        let dist = WeightedIndex::new(&p).expect("Invalid probabilities");
        let mut rng = thread_rng();
        (0..size).map(|_| dist.sample(&mut rng)).collect()
    }
}

/// Return the vector between a and b
pub fn vector(a: &Vec3<f32>, b: &Vec3<f32>) -> Vec3<f32> {
    b - a
}

/// Returm angle between a (can be multiple coordinates), b and c
pub fn get_angle(a: &Vec3<f32>, b: &Vec3<f32>, c: &Vec3<f32>, degree: bool) -> f32 {
    let ba = vector(b, a);
    let bc = vector(b, c);

    let mut cos_angle = ba.dot(bc) / (ba.magnitude() * bc.magnitude());
    // Make sure values fit between -1 and 1 for arccos
    cos_angle = cos_angle.clamp(-1., 1.);
    let angle = cos_angle.acos();

    if degree {
        angle.to_degrees()
    } else {
        angle
    }
}

// TODO port: review callers of this and argument types
pub fn boltzmann_acceptance_rejection(
    new_energies: &[f32],
    old_energies: &[f32],
    temperature: f32,
) -> Vec<bool> {
    let kb = 0.0019872041;

    // hmm nothing to flatten since &[f32]..
    // new_energies = np.ravel(new_energies)
    // old_energies = np.ravel(old_energies)

    let mut decisions = new_energies
        .iter()
        .zip(old_energies.iter())
        .map(|(&new, &old)| new < old);

    if decisions.all(|d| d) {
        return decisions.collect();
    } else {
        let unfavorable_indices = decisions.clone().enumerate().filter_map(|(i, d)| Some(i));

        let unfavorable_old_energies = if old_energies.len() == 1 {
            old_energies.to_vec()
        } else {
            unfavorable_indices
                .clone()
                .map(|i| old_energies[i])
                .collect()
        };

        // let delta_e = unfavorable_indices.map(|i| new_energies[i])
        let delta_e: Vec<f32> = unfavorable_indices
            .clone()
            .map(|i| new_energies[i]) // Extract elements from new_energies using indices
            .zip(unfavorable_old_energies.iter())
            .map(|(new, &old)| new - old) // Compute the difference
            .collect();

        let p_acc: Vec<f32> = delta_e
            .iter()
            .map(|d| (-d / (kb * temperature)).exp().min(1.))
            .collect();

        // r = np.random.rand(p_acc.shape[0])

        let mut rng = rand::thread_rng();
        let r: Vec<f32> = (0..p_acc.len()).map(|_| rng.gen::<f32>()).collect();

        // decisions[unfavorable_indices] = r <= p_acc

        // Clone the original decisions
        let mut updated_decisions: Vec<bool> = decisions.collect();

        // Update the decisions for the unfavorable indices
        for (index, (&rand_val, &p_val)) in unfavorable_indices.zip(r.iter().zip(p_acc.iter())) {
            updated_decisions[index] = rand_val <= p_val;
        }

        updated_decisions
    }
}

pub fn quaternion_rotate(y: &[Vec3<f32>], x: &[Vec3<f32>]) -> Vec4<f32> {
    assert_eq!(
        y.len(),
        x.len(),
        "Y and X must have the same number of points."
    );

    let n = y.len();

    // Initialize a 4x4 matrix for accumulation
    let mut a = Matrix4::<f32>::zeros();

    // Compute the summation of Qt_dot_W
    for i in 0..n {
        let w = make_w(y[i]);
        let q = make_q(x[i]);

        let qt_dot_w = q.transpose() * w;
        a += qt_dot_w;
    }

    // Perform eigenvalue decomposition on the symmetric matrix
    let eigen = a.symmetric_eigen();

    // Extract the eigenvector corresponding to the largest eigenvalue
    let max_index = eigen
        .eigenvalues
        .iter()
        .enumerate()
        .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
        .map(|(i, _)| i)
        .unwrap();

    // Convert the resulting eigenvector to Vec4<f32>
    let eigenvector = eigen.eigenvectors.column(max_index);
    Vec4::new(
        eigenvector[0],
        eigenvector[1],
        eigenvector[2],
        eigenvector[3],
    )
}

/// Create the W matrix
#[rustfmt::skip]
fn make_w(v: Vec3<f32>) -> Matrix4<f32> {
    Matrix4::new(
        0.0, -v.z,  v.y,  v.x,
        v.z,  0.0, -v.x,  v.y,
       -v.y,  v.x,  0.0,  v.z,
       -v.x, -v.y, -v.z,  0.0,
    )
}

/// Create the Q matrix
#[rustfmt::skip]
fn make_q(v: Vec3<f32>) -> Matrix4<f32> {
    Matrix4::new(
        0.0, -v.z, -v.y,  v.x,
        v.z,  0.0, -v.x,  v.y,
        v.y, -v.x,  0.0,  v.z,
        v.x,  v.y, -v.z,  0.0,
    )
}
