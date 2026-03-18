# --- Tests for unified train/test functions ---
# Verify that the unified functions with mode='train' and mode='test'
# produce results in the expected slots and are self-consistent.

# Reuse the trained + tested object from test-testing.R context
train_ind <- which(bjarn_data$group %in% unique(bjarn_data$group)[1:8])
test_ind  <- which(bjarn_data$group %in% unique(bjarn_data$group)[9:10])

tt_obj <- suppressMessages(quiet(train_model(
  exp_matrix = bjarn_data$expr_mat[, train_ind],
  genes      = bjarn_data$probes_used,
  group_1    = bjarn_data$group[train_ind],
  time       = bjarn_data$time[train_ind],
  log_thresh = -5
)))

tt_obj <- suppressMessages(quiet(test_model(
  tt_obj,
  exp_matrix   = bjarn_data$expr_mat[, test_ind],
  test_group_1 = bjarn_data$group[test_ind],
  test_time    = bjarn_data$time[test_ind],
  log_thresh   = -5
)))

test_that("get_final_likelis writes to correct slots for train", {
  expect_true(!is.null(tt_obj$Train_Data$Averaged_Likelis_Post_Thresh_Train))
  expect_true(!is.null(tt_obj$Train_Data$Max_Likelis_Train))
  expect_equal(tt_obj$Train_Data$LogThresh_Train, -5)
})

test_that("get_final_likelis writes to correct slots for test", {
  expect_true(!is.null(tt_obj$Test_Data$Averaged_Likelis_Post_Thresh_Test))
  expect_true(!is.null(tt_obj$Test_Data$Max_Likelis_Test))
  expect_equal(tt_obj$Test_Data$LogThresh_Test, -5)
})

test_that("theta_calc produces bounded values for both modes", {
  train_thetas <- tt_obj$Train_Data$Thetas_Train
  test_thetas <- tt_obj$Test_Data$Thetas_Test
  expect_true(all(train_thetas >= 0 & train_thetas <= 1))
  expect_true(all(test_thetas >= 0 & test_thetas <= 1))
})

test_that("second_peaks_fun produces Results_df for both modes", {
  train_res <- tt_obj$Train_Data$Results_df
  test_res <- tt_obj$Test_Data$Results_df

  # Both should have the core columns
  core_cols <- c("npeaks", "time_1st_peak", "Theta", "PercFlat", "FlatContrib", "Pred_Error", "Actual_Time")
  expect_true(all(core_cols %in% colnames(train_res)))
  expect_true(all(core_cols %in% colnames(test_res)))
})

test_that("backward-compatible wrappers call unified functions", {
  # These should not error — they just delegate to the unified version
  # Re-run on a fresh copy to verify wrappers work
  obj_copy <- tt_obj
  obj_copy <- get_final_likelis_train(obj_copy, log_thresh = -5)
  expect_true(!is.null(obj_copy$Train_Data$Averaged_Likelis_Post_Thresh_Train))

  obj_copy <- get_final_likelis_test(obj_copy, log_thresh = -5)
  expect_true(!is.null(obj_copy$Test_Data$Averaged_Likelis_Post_Thresh_Test))
})

test_that("choose_logthresh returns expected data frame", {
  lt_df <- suppressMessages(quiet(choose_logthresh(tt_obj, max_log = 0, min_log = -4, by_step = -2, train_or_test = 'train')))
  expect_s3_class(lt_df, "data.frame")
  expected_cols <- c("Perc_Flat", "PeakNum_Ratio", "MaxLik_Ratio", "Mean_Theta", "Median_Theta", "FlatContrib", "LogThresh")
  expect_true(all(expected_cols %in% colnames(lt_df)))
  # LogThresh column should match the requested sequence
  expect_equal(lt_df$LogThresh, c(0, -2, -4))
})
