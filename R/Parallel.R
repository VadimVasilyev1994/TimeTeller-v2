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
    lrf_curve_spline <- stats::predict(splines::periodicSpline(1:num_points, curr_lrf_curve, period = num_points), seq(1, num_points, length.out = 1000))
    curve_spline <- stats::predict(splines::periodicSpline(1:num_points, curr_curve, period = num_points), seq(1, num_points, length.out = 1000))
    sum(lrf_curve_spline$y > curve_spline$y) / length(lrf_curve_spline$y)
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
    lrf_curve_spline <- stats::predict(splines::periodicSpline(1:num_points, ind_lrf_curve, period = num_points), seq(1, num_points, length.out = 1000))
    curve_spline <- stats::predict(splines::periodicSpline(1:num_points, curr_curve, period = num_points), seq(1, num_points, length.out = 1000))
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

calc_train_likelis_dev_test <- function(object) {
  svd_data <- object[['Projections']][['SVD_Per_Time_Point']]
  fitted_mvn_data <- object[['Projections']][['Fitted_MVN_Interpolated']]
  train_exp_data <- object[['Train']][['Normalised_Train_Exp_Data']]
  num_PC <- object[['PC_Num']]

  out_list <- foreach(i = 1:length(names(svd_data)), .combine = 'c', .inorder = TRUE, .multicombine = TRUE, .packages = c('Matrix','mvtnorm','foreach')) %dopar% {
    project_exp_mat <- svd_data[[i]] %*% train_exp_data
    mat <- matrix(NA, nrow = dim(fitted_mvn_data)[2], ncol = dim(train_exp_data)[2])
    for (ind_num in 1:dim(train_exp_data)[2]) {
      vec <- numeric(length = dim(fitted_mvn_data)[2])
      for (j in 1:dim(fitted_mvn_data)[2]) {
        curr_sigma <- matrix(fitted_mvn_data[(num_PC+1):(num_PC+num_PC^2),j,i], nrow = num_PC)
        curr_eig <- eigen(curr_sigma, symmetric = TRUE, only.values = TRUE)$values
        if (any(curr_eig < 0)) {
          sigma_used <- nearPD(curr_sigma, base.matrix = TRUE, ensureSymmetry = TRUE, eig.tol = 1e-05, conv.tol = 1e-06, posd.tol = 1e-07)$mat
        } else {
          sigma_used <- curr_sigma
        }

        vec[j] <- mvtnorm::dmvnorm(project_exp_mat[,ind_num], mean = fitted_mvn_data[1:num_PC,j,i], sigma = sigma_used, checkSymmetry = FALSE, log = TRUE)

      }
      mat[ ,ind_num] <- vec
    }
  mat
  }

  object[['Train_Data']][['Train_Likelihood_Array']] <- array(out_list, dim = c(dim(fitted_mvn_data)[2], dim(train_exp_data)[2], dim(fitted_mvn_data)[3]))
  return(object)
}


calc_test_likelis_dev_test <- function(object) {
  svd_data <- object[['Projections']][['SVD_Per_Time_Point']]
  fitted_mvn_data <- object[['Projections']][['Fitted_MVN_Interpolated']]
  test_exp_data <- object[['Test_Data']][['Normalised_Test_Exp_Data']]
  num_PC <- object[['PC_Num']]

  out_list <- foreach(i = 1:length(names(svd_data)), .combine = 'c', .inorder = TRUE, .multicombine = TRUE, .packages = c('Matrix','mvtnorm','foreach')) %dopar% {
    project_exp_mat <- svd_data[[i]] %*% test_exp_data
    mat <- matrix(NA, nrow = dim(fitted_mvn_data)[2], ncol = dim(test_exp_data)[2])
    for (ind_num in 1:dim(test_exp_data)[2]) {
      vec <- numeric(length = dim(fitted_mvn_data)[2])
      for (j in 1:dim(fitted_mvn_data)[2]) {
        curr_sigma <- matrix(fitted_mvn_data[(num_PC+1):(num_PC+num_PC^2),j,i], nrow = num_PC)
        curr_eig <- eigen(curr_sigma, symmetric = TRUE, only.values = TRUE)$values
        if (any(curr_eig < 0)) {
          sigma_used <- nearPD(curr_sigma, base.matrix = TRUE, ensureSymmetry = TRUE, eig.tol = 1e-05, conv.tol = 1e-06, posd.tol = 1e-07)$mat
        } else {
          sigma_used <- curr_sigma
        }

        vec[j] <- mvtnorm::dmvnorm(project_exp_mat[,ind_num], mean = fitted_mvn_data[1:num_PC,j,i], sigma = sigma_used, checkSymmetry = FALSE, log = TRUE)

      }
      mat[ ,ind_num] <- vec
    }
    mat
  }

  object[['Test_Data']][['Test_Likelihood_Array']] <- array(out_list, dim = c(dim(fitted_mvn_data)[2], dim(test_exp_data)[2], dim(fitted_mvn_data)[3]))
  return(object)
}
