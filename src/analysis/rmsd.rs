//! Root-Mean-Square Deviation (RMSD) and protein superposition algorithms

use nalgebra::{Matrix3, SVD, Vector3};

/// Calculates the Root-Mean-Square Deviation (RMSD) between two sets of coordinates
pub fn calculate_rmsd_between_coordinate_sets(
    reference_coordinates: &[glam::Vec3],
    moving_coordinates: &[glam::Vec3],
) -> Result<f32, String> {
    if reference_coordinates.len() != moving_coordinates.len() {
        return Err(format!(
            "Coordinate sets must have the same length for RMSD calculation ({} vs {} coordinates).",
            reference_coordinates.len(),
            moving_coordinates.len()
        ));
    }
    
    if reference_coordinates.is_empty() {
        return Ok(0.0);
    }

    let mut sum_of_squared_distances = 0.0;
    for (current_coordinate_index, reference_point) in reference_coordinates.iter().enumerate() {
        let moving_point = moving_coordinates[current_coordinate_index];
        sum_of_squared_distances += reference_point.distance_squared(moving_point);
    }

    Ok((sum_of_squared_distances / reference_coordinates.len() as f32).sqrt())
}

/// Computes the optimal rotation matrix to superimpose 'moving' onto 'reference' using Kabsch algorithm

pub fn compute_kabsch_optimal_rotation(

    reference_centered_coordinates: &[glam::Vec3],

    moving_centered_coordinates: &[glam::Vec3],

) -> Result<glam::Mat3, String> {

    if reference_centered_coordinates.len() != moving_centered_coordinates.len() {

        return Err(format!(

            "Coordinate sets must have the same length for Kabsch algorithm ({} vs {} coordinates).",

            reference_centered_coordinates.len(),

            moving_centered_coordinates.len()

        ));

    }



    if reference_centered_coordinates.is_empty() {

        return Ok(glam::Mat3::IDENTITY);

    }



    let mut covariance_matrix = Matrix3::<f64>::zeros();

    for (current_coordinate_index, reference_point) in reference_centered_coordinates.iter().enumerate() {

        let moving_point = moving_centered_coordinates[current_coordinate_index];

        

        let reference_vector = Vector3::new(reference_point.x as f64, reference_point.y as f64, reference_point.z as f64);

        let moving_vector = Vector3::new(moving_point.x as f64, moving_point.y as f64, moving_point.z as f64);

        

        covariance_matrix += moving_vector * reference_vector.transpose();

    }



    let svd_result = SVD::new(covariance_matrix, true, true);

    let u_matrix = svd_result.u.ok_or("SVD failed to produce U matrix")?;

    let v_t_matrix = svd_result.v_t.ok_or("SVD failed to produce V_t matrix")?;



    let mut determinant_sign_correction = (v_t_matrix.transpose() * u_matrix.transpose()).determinant();

    if determinant_sign_correction < 0.0 {

        determinant_sign_correction = -1.0;

    } else {

        determinant_sign_correction = 1.0;

    }



    let mut reflection_matrix = Matrix3::<f64>::identity();

    reflection_matrix[(2, 2)] = determinant_sign_correction;



    let optimal_rotation_matrix = v_t_matrix.transpose() * reflection_matrix * u_matrix.transpose();

    

    Ok(glam::Mat3::from_cols_array(&[

        optimal_rotation_matrix[(0, 0)] as f32, optimal_rotation_matrix[(1, 0)] as f32, optimal_rotation_matrix[(2, 0)] as f32,

        optimal_rotation_matrix[(0, 1)] as f32, optimal_rotation_matrix[(1, 1)] as f32, optimal_rotation_matrix[(2, 1)] as f32,

        optimal_rotation_matrix[(0, 2)] as f32, optimal_rotation_matrix[(1, 2)] as f32, optimal_rotation_matrix[(2, 2)] as f32,

    ]))

}
