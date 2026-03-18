# --- Tests for plotting functions ---
# These verify that plots are constructed without error and return
# the expected object types. They do NOT test visual appearance.

train_ind <- which(bjarn_data$group %in% unique(bjarn_data$group)[1:8])
test_ind  <- which(bjarn_data$group %in% unique(bjarn_data$group)[9:10])

tt_plot <- suppressMessages(quiet(train_model(
  exp_matrix = bjarn_data$expr_mat[, train_ind],
  genes      = bjarn_data$probes_used,
  group_1    = bjarn_data$group[train_ind],
  time       = bjarn_data$time[train_ind],
  log_thresh = -5
)))

tt_plot <- suppressMessages(quiet(test_model(
  tt_plot,
  exp_matrix   = bjarn_data$expr_mat[, test_ind],
  test_group_1 = bjarn_data$group[test_ind],
  test_time    = bjarn_data$time[test_ind],
  log_thresh   = -5
)))

test_that("plot_reps returns a ggplot", {
  p <- plot_reps(tt_plot, gene = bjarn_data$probes_used[1])
  expect_s3_class(p, "ggplot")
})

test_that("plot_genes returns a ggplot", {
  p <- plot_genes(tt_plot, genes = bjarn_data$probes_used[1:3])
  expect_s3_class(p, "ggplot")
})

test_that("exprs_vs_PredTime_plot returns a ggplot", {
  p <- suppressWarnings(exprs_vs_PredTime_plot(
    tt_plot, genes = bjarn_data$probes_used[1:2],
    theta_thresh = 0.3, xlim_l = 0, xlim_u = 24
  ))
  expect_s3_class(p, "ggplot")
})

test_that("plotPCs_test returns a ggplot", {
  p <- plotPCs_test(tt_plot)
  expect_s3_class(p, "ggplot")
})

test_that("choose_logthresh_plot returns a ggplot", {
  lt_df <- suppressMessages(quiet(choose_logthresh(tt_plot, max_log = 0, min_log = -4, by_step = -2, train_or_test = 'train')))
  p <- choose_logthresh_plot(lt_df)
  expect_s3_class(p, "ggplot")
})

test_that("plot_raw_likelis runs without error", {
  # Base R plot — just check it doesn't error
  expect_no_error(plot_raw_likelis(tt_plot, sample_num = 1, logthresh = -5, train_or_test = 'test'))
})

test_that("plot_ind_curve runs without error", {
  expect_no_error(plot_ind_curve(tt_plot, sample_num = 1, logthresh = -5, train_or_test = 'test'))
})

test_that("theme_tt returns a ggplot theme", {
  th <- theme_tt()
  expect_s3_class(th, "theme")
})
