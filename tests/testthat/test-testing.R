# --- Tests for the testing (prediction) pipeline ---

# Train on subset, test on held-out
train_ind <- which(bjarn_data$group %in% unique(bjarn_data$group)[1:8])
test_ind  <- which(bjarn_data$group %in% unique(bjarn_data$group)[9:10])

tt_full <- suppressMessages(quiet(train_model(
  exp_matrix = bjarn_data$expr_mat[, train_ind],
  genes      = bjarn_data$probes_used,
  group_1    = bjarn_data$group[train_ind],
  time       = bjarn_data$time[train_ind],
  log_thresh = -5
)))

tt_full <- suppressMessages(quiet(test_model(
  tt_full,
  exp_matrix   = bjarn_data$expr_mat[, test_ind],
  test_group_1 = bjarn_data$group[test_ind],
  test_time    = bjarn_data$time[test_ind],
  log_thresh   = -5
)))

test_that("test_model returns expected structure", {
  expect_true("Test_Data" %in% names(tt_full))
  expect_true("Results_df" %in% names(tt_full$Test_Data))
  expect_true("Test_Likelihood_Array" %in% names(tt_full$Test_Data))
  expect_true("Test_Projections" %in% names(tt_full$Test_Data))
})

test_that("test results have correct number of samples", {
  n_test <- length(test_ind)
  expect_equal(nrow(tt_full$Test_Data$Results_df), n_test)
  expect_equal(length(tt_full$Test_Data$Thetas_Test), n_test)
})

test_that("test theta values are bounded [0, 1]", {
  thetas <- tt_full$Test_Data$Results_df$Theta
  expect_true(all(thetas >= 0 & thetas <= 1))
})

test_that("test predictions are in [0, 24)", {
  pred_times <- tt_full$Test_Data$Results_df$time_1st_peak
  expect_true(all(pred_times >= 0 & pred_times < 24, na.rm = TRUE))
})

test_that("test prediction error is reasonable", {
  # On bjarn_data with 8 training / 2 test individuals, MAE should be < 6h
  mae <- mean(abs(tt_full$Test_Data$Results_df$Pred_Error))
  expect_true(mae < 6)
})

test_that("test metadata is stored correctly", {
  meta <- tt_full$Metadata$Test
  expect_equal(length(meta$Group_1), length(test_ind))
  expect_equal(length(meta$Time), length(test_ind))
})

test_that("test likelihood array dimensions are correct", {
  likelis <- tt_full$Test_Data$Test_Likelihood_Array
  expect_equal(length(dim(likelis)), 3)
  expect_equal(dim(likelis)[2], length(test_ind))
  # Values should be finite or -Inf (log(0) at zero-density interpolation points)
  expect_true(all(!is.na(likelis) & !is.nan(likelis)))
})

test_that("normalise_test_data produces correct dimensions for intergene", {
  # intergene normalisation should preserve gene x sample dimensions
  norm_mat <- tt_full$Test_Data$Normalised_Test_Exp_Data
  expect_equal(nrow(norm_mat), length(bjarn_data$probes_used))
  expect_equal(ncol(norm_mat), length(test_ind))
})
