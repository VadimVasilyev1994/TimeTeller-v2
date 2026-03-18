# --- Tests for utility functions in Utils.R ---

test_that("circular_difference is symmetric and bounded", {
  # Symmetric: diff(a,b) == -diff(b,a)
  expect_equal(circular_difference(6, 18), -circular_difference(18, 6))

  # Bounded within [-12, 12] hours
  diffs <- circular_difference(seq(0, 23, by = 1), rep(12, 24))
  expect_true(all(diffs >= -12 & diffs <= 12))

  # Identical times give zero
  expect_equal(circular_difference(5, 5), 0)

  # Wrap-around: 23h and 1h are 2h apart, not 22h
  expect_equal(abs(circular_difference(23, 1)), 2)
})

test_that("shift_ts produces correct circular shift", {
  ts <- c(1, 2, 5, 3, 1)  # peak at position 3
  shifted <- shift_ts(ts, target_peak = 1)
  # Peak should now be at position 1
  expect_equal(which.max(shifted), 1)
  # Length preserved
  expect_equal(length(shifted), length(ts))
  # Values preserved (just reordered)
  expect_equal(sort(shifted), sort(ts))
})

test_that("check_suggested_pkg errors on missing package", {
  expect_error(check_suggested_pkg("nonexistent_package_xyz"),
               "is required for this functionality")
})

test_that("normalize_per_sample produces zero mean and unit SD per row", {
  mat <- matrix(c(10, 20, 30, 5, 15, 25), nrow = 2, byrow = TRUE)
  colnames(mat) <- c("A", "B", "C")
  result <- normalize_per_sample(as.data.frame(mat))

  # Each row should have mean ~0 and sd ~1
  row_means <- apply(result, 1, mean)
  row_sds <- apply(result, 1, sd)
  expect_equal(row_means, c(0, 0), tolerance = 1e-10)
  expect_equal(row_sds, c(1, 1), tolerance = 1e-10)
})

test_that("nature_palette returns correct number of colours", {
  expect_length(nature_palette(), 6)
  expect_length(nature_palette(n = 10), 10)
  expect_length(nature_palette(n = 3, type = "diverging"), 3)
  # All hex colours
  expect_true(all(grepl("^#[0-9A-F]{6}$", nature_palette(), ignore.case = TRUE)))
})

test_that("get_group_indices returns correct indices", {
  # Create a minimal mock object
  mock_obj <- list(
    Full_Original_Data = matrix(1:20, nrow = 5, ncol = 4),
    Metadata = list(Train = list(
      Group_1 = c("A", "A", "B", "B"),
      Group_2 = c("X", "Y", "X", "Y"),
      Group_3 = rep(NA, 4),
      Replicate = rep(NA, 4)
    ))
  )
  # All samples when no filter
  expect_equal(get_group_indices(mock_obj), 1:4)
  # Filter by group1
  expect_equal(get_group_indices(mock_obj, group1 = "A"), c(1, 2))
  # Filter by group1 + group2 intersection
  expect_equal(get_group_indices(mock_obj, group1 = "A", group2 = "Y"), 2)
})
