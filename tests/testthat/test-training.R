# --- Tests for the training pipeline ---

# Train once for all tests in this file (expensive operation)
tt_trained <- suppressMessages(quiet(train_model(
  exp_matrix = bjarn_data$expr_mat,
  genes      = bjarn_data$probes_used,
  group_1    = bjarn_data$group,
  time       = bjarn_data$time,
  log_thresh = -5
)))

test_that("train_model returns a list with expected structure", {
  expect_type(tt_trained, "list")
  expect_true("Train" %in% names(tt_trained))
  expect_true("Train_Data" %in% names(tt_trained))
  expect_true("Projections" %in% names(tt_trained))
  expect_true("Metadata" %in% names(tt_trained))
})

test_that("training metadata is stored correctly", {
  meta <- tt_trained$Metadata$Train
  expect_equal(meta$Genes_Used, bjarn_data$probes_used)
  expect_equal(length(meta$Time), ncol(bjarn_data$expr_mat))
  expect_equal(length(meta$Group_1), ncol(bjarn_data$expr_mat))
})

test_that("train_model produces valid theta values", {
  thetas <- tt_trained$Train_Data$Results_df$Theta
  # Theta must be between 0 and 1
  expect_true(all(thetas >= 0 & thetas <= 1))
  # One theta per sample
  expect_equal(length(thetas), ncol(bjarn_data$expr_mat))
})

test_that("train_model produces valid time predictions", {
  pred_times <- tt_trained$Train_Data$Results_df$time_1st_peak
  # Predictions must be in [0, 24)
  expect_true(all(pred_times >= 0 & pred_times < 24, na.rm = TRUE))
  # Prediction errors should be reasonable (< 6h mean for a good training set)
  mae <- mean(abs(tt_trained$Train_Data$Results_df$Pred_Error))
  expect_true(mae < 6)
})

test_that("normalisation produces expected dimensions", {
  norm_data <- tt_trained$Train$Normalised_Train_Exp_Data
  # Rows = genes, cols = samples
  expect_equal(nrow(norm_data), length(bjarn_data$probes_used))
  expect_equal(ncol(norm_data), ncol(bjarn_data$expr_mat))
})

test_that("SVD projections have correct structure", {
  svd_data <- tt_trained$Projections$SVD_Per_Time_Point
  # One projection per time point
  n_timepoints <- length(unique(bjarn_data$time))
  expect_equal(length(svd_data), n_timepoints)
  # Each projection matrix: num_PC x num_genes
  expect_equal(nrow(svd_data[[1]]), tt_trained$PC_Num)
  expect_equal(ncol(svd_data[[1]]), length(bjarn_data$probes_used))
})

test_that("likelihood array has correct dimensions", {
  likelis <- tt_trained$Train_Data$Train_Likelihood_Array
  # 3D array: [interp_points x samples x local_projections]
  expect_equal(length(dim(likelis)), 3)
  expect_equal(dim(likelis)[2], ncol(bjarn_data$expr_mat))
  # Values should be finite or -Inf (log(0) at zero-density interpolation points)
  expect_true(all(!is.na(likelis) & !is.nan(likelis)))
})

test_that("flat contribution values are bounded", {
  flat_df <- tt_trained$Train_Data$Flat_Contrib_to_Theta_df
  expect_true(all(flat_df$flat_contributions >= 0 & flat_df$flat_contributions <= 100))
})

test_that("epsilon and eta hyperparameters are stored", {
  expect_equal(tt_trained$Train_Data$epsilon, 0.4)
  expect_equal(tt_trained$Train_Data$eta, 0.35)
})
