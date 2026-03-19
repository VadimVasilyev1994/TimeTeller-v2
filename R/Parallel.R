## calc_train_likelis_dev removed — superseded by calc_train_likelis_dev_test (below)

## ---------------------------------------------------------------------------
## Parallel versions of theta and flat contribution (unified train/test)
## ---------------------------------------------------------------------------

theta_calc_parallel <- function(object, mode = 'train', epsilon = NULL, eta = NULL) {
  if (mode == 'train') {
    averaged_likelis_rescaled <- t(object[['Train_Data']][['Averaged_Likelis_Post_Thresh_Train']])
  } else {
    epsilon <- object[['Train_Data']][['epsilon']]
    eta <- object[['Train_Data']][['eta']]
    averaged_likelis_rescaled <- t(object[['Test_Data']][['Averaged_Likelis_Post_Thresh_Test']])
  }

  num_samples <- dim(averaged_likelis_rescaled)[1]
  num_points <- dim(averaged_likelis_rescaled)[2]

  thetas <- foreach(i = 1:num_samples, .combine = 'c', .inorder = TRUE, .packages = c('stats','splines')) %dopar% {
    curr_sample <- exp(averaged_likelis_rescaled[i,])
    curr_lrf_curve <- curr_sample / max(curr_sample)
    peak_pos <- which.max(curr_lrf_curve)
    curr_curve <- suppressWarnings(
      eta * (1 + epsilon + cos(2 * pi * ((1:num_points) / num_points - peak_pos / num_points)))
    )
    eval_seq <- seq(1, num_points, length.out = 500L)
    lrf_curve_spline <- stats::predict(splines::periodicSpline(1:num_points, curr_lrf_curve, period = num_points), eval_seq)
    curve_spline <- stats::predict(splines::periodicSpline(1:num_points, curr_curve, period = num_points), eval_seq)
    sum(lrf_curve_spline$y > curve_spline$y) / 500L
  }

  if (mode == 'train') {
    object[['Train_Data']][['Thetas_Train']] <- thetas
    object[['Train_Data']][['epsilon']] <- epsilon
    object[['Train_Data']][['eta']] <- eta
  } else {
    object[['Test_Data']][['Thetas_Test']] <- thetas
  }
  return(object)
}

# Backward-compatible wrappers
theta_calc_train_dev <- function(object, epsilon = 0.4, eta = 0.35) {
  theta_calc_parallel(object, mode = 'train', epsilon = epsilon, eta = eta)
}
theta_calc_test_dev <- function(object) {
  theta_calc_parallel(object, mode = 'test')
}

calc_flat_theta_contrib_parallel <- function(object, mode = 'train') {
  epsilon <- object[['Train_Data']][['epsilon']]
  eta <- object[['Train_Data']][['eta']]

  if (mode == 'train') {
    averaged_likelis_rescaled <- t(object[['Train_Data']][['Averaged_Likelis_Post_Thresh_Train']])
    thetas_vec <- object[['Train_Data']][['Thetas_Train']]
  } else {
    averaged_likelis_rescaled <- t(object[['Test_Data']][['Averaged_Likelis_Post_Thresh_Test']])
    thetas_vec <- object[['Test_Data']][['Thetas_Test']]
  }

  num_samples <- dim(averaged_likelis_rescaled)[1]
  num_points <- dim(averaged_likelis_rescaled)[2]

  flat_contribution <- foreach(i = 1:num_samples, .combine = 'c', .inorder = TRUE, .export = 'shift_ts', .packages = c('stats','splines')) %dopar% {
    ind_ts <- exp(averaged_likelis_rescaled[i,])
    ind_ts <- shift_ts(ind_ts, round(num_points / 2))
    ind_lrf_curve <- ind_ts / max(ind_ts)
    peak_pos <- which.max(ind_lrf_curve)
    curr_curve <- suppressWarnings(
      eta * (1 + epsilon + cos(2 * pi * ((1:num_points) / num_points - peak_pos / num_points)))
    )
    eval_seq <- seq(1, num_points, length.out = 500L)
    lrf_curve_spline <- stats::predict(splines::periodicSpline(1:num_points, ind_lrf_curve, period = num_points), eval_seq)
    curve_spline <- stats::predict(splines::periodicSpline(1:num_points, curr_curve, period = num_points), eval_seq)
    theta_index <- which(lrf_curve_spline$y > curve_spline$y)
    flat_regions <- which(abs(diff(lrf_curve_spline$y)) < 1e-4)
    sum(theta_index %in% flat_regions) / length(theta_index) * 100
  }

  flat_contribution_df <- data.frame(flat_contributions = flat_contribution, thetas = thetas_vec)

  if (mode == 'train') {
    object[['Train_Data']][['Flat_Contrib_to_Theta_df']] <- flat_contribution_df
  } else {
    object[['Test_Data']][['Flat_Contrib_to_Theta_df']] <- flat_contribution_df
  }
  return(object)
}

# Backward-compatible wrappers
calc_flat_theta_contrib_train_dev <- function(object) {
  calc_flat_theta_contrib_parallel(object, mode = 'train')
}
calc_flat_theta_contrib_test_dev <- function(object) {
  calc_flat_theta_contrib_parallel(object, mode = 'test')
}

## theta_calc_test_dev and calc_flat_theta_contrib_test_dev are now
## handled by the unified parallel functions above with wrappers.


###
# NOTE: These parallel likelihood functions use dmvnorm(..., log = TRUE) directly,
# so the output is already on the log scale. This differs from the serial
# calc_likelis() which calls dmvnorm() without log=TRUE and wraps in log() after.
# Both produce identical log-likelihood arrays.

calc_likelis_parallel <- function(object, mode = 'train') {
  svd_data <- object[['Projections']][['SVD_Per_Time_Point']]
  fitted_mvn_data <- object[['Projections']][['Fitted_MVN_Interpolated']]
  num_PC <- object[['PC_Num']]

  if (mode == 'train') {
    exp_data <- object[['Train']][['Normalised_Train_Exp_Data']]
  } else {
    exp_data <- object[['Test_Data']][['Normalised_Test_Exp_Data']]
  }

  n_interp <- dim(fitted_mvn_data)[2]

  out_list <- foreach(i = 1:length(names(svd_data)), .combine = 'c', .inorder = TRUE, .multicombine = TRUE, .packages = c('Matrix','mvtnorm','foreach')) %dopar% {
    project_exp_mat <- svd_data[[i]] %*% exp_data
    samples_mat <- t(project_exp_mat)  # n_samples x num_PC
    mat <- matrix(NA, nrow = n_interp, ncol = dim(exp_data)[2])

    # Pre-compute PD-corrected covariance matrices (once per interp point)
    for (j in 1:n_interp) {
      curr_sigma <- matrix(fitted_mvn_data[(num_PC+1):(num_PC+num_PC^2),j,i], nrow = num_PC)
      curr_eig <- eigen(curr_sigma, symmetric = TRUE, only.values = TRUE)$values
      if (any(curr_eig < 0)) {
        sigma_used <- as.matrix(nearPD(curr_sigma, base.matrix = TRUE, ensureSymmetry = TRUE, eig.tol = 1e-05, conv.tol = 1e-06, posd.tol = 1e-07)$mat)
      } else {
        sigma_used <- curr_sigma
      }
      # Vectorised: evaluate density for all samples at once
      mat[j, ] <- mvtnorm::dmvnorm(samples_mat, mean = fitted_mvn_data[1:num_PC,j,i], sigma = sigma_used, checkSymmetry = FALSE, log = TRUE)
    }
    mat
  }

  likelihood_array <- array(out_list, dim = c(n_interp, dim(exp_data)[2], dim(fitted_mvn_data)[3]))

  if (mode == 'train') {
    object[['Train_Data']][['Train_Likelihood_Array']] <- likelihood_array
  } else {
    object[['Test_Data']][['Test_Likelihood_Array']] <- likelihood_array
  }
  return(object)
}

# Backward-compatible wrappers
calc_train_likelis_dev_test <- function(object) {
  calc_likelis_parallel(object, mode = 'train')
}
calc_test_likelis_dev_test <- function(object) {
  calc_likelis_parallel(object, mode = 'test')
}
