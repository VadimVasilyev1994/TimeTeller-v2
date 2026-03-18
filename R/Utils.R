perpchipend <- function(h1, h2, del1, del2) {
  d <- ((2 * h1 + h2) * del1 - h1 * del2)/(h1 + h2)
  if (sign(d) != sign(del1)) {
    d <- 0
  }
  else if ((sign(del1) != sign(del2)) && (abs(d) > abs(3 * del1))) {
    d <- 3 * del1
  }
  return(d)
}

perpchip_slopes <- function(h, delta) {
  n <- length(h) + 1
  d <- numeric(length(h))
  k <- which(sign(delta[1:(n - 2)]) * sign(delta[2:(n - 1)]) > 0) + 1
  w1 <- 2 * h[k] + h[k - 1]
  w2 <- h[k] + 2 * h[k - 1]
  d[k] <- (w1 + w2)/(w1/delta[k - 1] + w2/delta[k])
  d[1] <- perpchipend(h[1], h[2], delta[1], delta[2])
  d[n] <- perpchipend(h[n - 1], h[n - 2], delta[n - 1], delta[n - 2])
  return(d)
}

perpchip <- function(xi,yi,x) {
  stopifnot(is.numeric(xi), is.numeric(yi), is.numeric(x))
  if (!pracma::is.sorted(xi))
    stop("Argument 'xi' must be a sorted vector of real numbers.")
  n <- length(xi)
  if (length(yi) != n)
    stop("Arguments 'xi', 'yi' must be vectors of equal length.")
  if (n <= 2)
    stop("At least three points needed for cubic interpolation.")
  h <- diff(xi)
  delta <- diff(yi) / h
  x_used <- xi
  y_used <- yi

  # Processing data to make it periodic in end conditons
  yi <- c(yi[length(yi)-1],yi,yi[2])
  xi <- c(xi[1]-(xi[length(xi)]-xi[length(xi)-1]),xi,xi[length(xi)]+(xi[2]-xi[1]))
  deltai <- c(delta[length(delta)], delta, delta[1])
  hi <- diff(xi)

  d <- perpchip_slopes(hi, deltai)
  d <- d[2:(length(d)-1)]

  # Fitting to the original non-extended data now
  a <- (3 * delta - 2 * d[1:(n - 1)] - d[2:n])/h
  b <- (d[1:(n - 1)] - 2 * delta + d[2:n])/h^2
  k <- rep(1, length(x))
  for (j in 2:(n - 1)) {
    k[x_used[j] <= x] <- j
  }
  s <- x - x_used[k]
  v <- y_used[k] + s * (d[k] + s * (a[k] + s * b[k]))
  return(v)
}

circular_difference <- function(time_1, time_2) {
  time_1_rad <- time_1 * 2 * pi / 24
  time_2_rad <- time_2 * 2 * pi / 24
  abs_diff_rad <- atan2(sin(time_1_rad-time_2_rad), cos(time_1_rad-time_2_rad))
  diff_hours <- abs_diff_rad * 24 / 2 / pi
  return(diff_hours)
}

# Subset sample indices by optional group filters.
# Returns the intersection of indices matching each non-NULL filter.
get_group_indices <- function(object, group1 = NULL, group2 = NULL, group3 = NULL, replicate = NULL) {
  n <- ncol(object[['Full_Original_Data']])
  meta <- object[['Metadata']][['Train']]
  idx_gr1 <- if (is.null(group1)) seq_len(n) else which(meta[['Group_1']] %in% group1)
  idx_gr2 <- if (is.null(group2)) seq_len(n) else which(meta[['Group_2']] %in% group2)
  idx_gr3 <- if (is.null(group3)) seq_len(n) else which(meta[['Group_3']] %in% group3)
  idx_rep <- if (is.null(replicate)) seq_len(n) else which(meta[['Replicate']] %in% replicate)
  Reduce(intersect, list(idx_gr1, idx_gr2, idx_gr3, idx_rep))
}

# Build a group label vector by pasting all metadata fields
build_group_vec <- function(object, index_used) {
  meta <- object[['Metadata']][['Train']]
  paste(meta[['Group_1']], meta[['Group_2']], meta[['Group_3']], meta[['Replicate']], sep = '_')[index_used]
}

# Circular shift: move the peak of ts to target_peak position
shift_ts <- function(ts, target_peak) {
  n <- length(ts)
  shift_dist <- target_peak - which.max(ts)
  indices <- ((seq_len(n) - 1 - shift_dist) %% n) + 1
  ts[indices]
}

extract_elements_func <- function(data) {
  return(c(unname(data$mu), unname(c(data$sigma))))
}

extract_elements_func_robust <- function(data) {
  return(c(unname(data$center), unname(c(data$cov))))
}

quiet <- function(x) {
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}

deviation_cv_corrected <- function(list_cv) {
  pred_errors <- purrr::map(list_cv, as_mapper(~ .x$Test_Data$Results_df$Pred_Error - median(.x$Test_Data$Results_df$Pred_Error)))
  return(unlist(pred_errors))
}

volume_func <- function(curr_sigma, alpha = 0.95, dims = 3) {
  crit_val <- qchisq(p = alpha, df = dims, lower.tail=FALSE) ^ (dims/2)
  det_sqrt <- sqrt(det(curr_sigma))
  gamma_func <- dims * gamma(dims/2)
  volume <- (2 * pi^(dims/2)) / gamma_func * crit_val^(dims/2) * det_sqrt
  return(volume)
}

give_names <- function(vec) {
  names(vec) <- paste0('PC',1:length(vec))
  return(vec)
}

## Helper: check that a suggested package is installed before use
check_suggested_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      "Package '", pkg, "' is required for this functionality but is not installed.\n",
      "Install it with: install.packages('", pkg, "')",
      call. = FALSE
    )
  }
}

prepare_raw_counts <- function(exp_matrix) {
  check_suggested_pkg('edgeR')
  obj_dge <- edgeR::DGEList(counts=as.matrix(exp_matrix))
  keep <- edgeR::filterByExpr(obj_dge)
  obj_dge <- obj_dge[keep, , keep.lib.sizes=FALSE]
  logCPM <- edgeR::cpm(obj_dge, log=TRUE, prior.count = 2)
  return(logCPM)
}

get_se_obj <- function(object) {
  check_suggested_pkg('SummarizedExperiment')
  counts_train <- object$Full_Original_Data
  counts_test <- object$Test_Data$Full_Test_Data
  common_names <- intersect(rownames(counts_train), rownames(counts_test))
  meta_train <- object$Train_Data$Results_df %>% dplyr::select(c(Group, Group_2, Group_3, Replicate, time_1st_peak, Actual_Time, Theta)) %>%
    dplyr::mutate(Dataset = 'Train')
  meta_test <- object$Test_Data$Results_df %>% dplyr::select(c(Group_1, Group_2, Group_3, Replicate, time_1st_peak, Actual_Time, Theta)) %>%
    dplyr::mutate(Group_1 = as.character(Group_1)) %>% dplyr::rename(Group = Group_1) %>% dplyr::mutate(Dataset = 'Test')
  combined_counts <- cbind(counts_train[common_names, ], counts_test[common_names, ])
  combined_meta <- rbind(meta_train, meta_test) %>% dplyr::mutate(across(where(is.character), ~na_if(., '')))
  SEobj <- SummarizedExperiment::SummarizedExperiment(assays = list(normalised_counts = combined_counts),
                                colData = combined_meta)
  return(SEobj)
}

get_tt_res <- function(object, train_or_test = 'test') {
  if (train_or_test == 'test') {
    return(object$Test_Data$Results_df)
  } else if (train_or_test == 'train') {
    return(object$Train_Data$Results_df)
  }
}

## ---------------------------------------------------------------------------
## Package-wide ggplot theme and colour palettes
## ---------------------------------------------------------------------------

# Tol Bright colourblind-safe palette (7 colours)
tt_palette <- c("#4477AA", "#EE6677", "#228833", "#CCBB44",
                "#66CCEE", "#AA3377", "#BBBBBB")

#' Nature-style colour palette
#'
#' Provides qualitative, sequential, or diverging colour palettes suitable
#' for publication figures.
#'
#' @param n number of colours to generate (NULL returns the base palette)
#' @param type one of 'qualitative', 'sequential_blue', or 'diverging'
#' @return character vector of hex colour codes
#' @export
nature_palette <- function(n = NULL, type = "qualitative") {
  palettes <- list(
    qualitative    = c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F", "#8491B4"),
    sequential_blue = c("#F7FBFF", "#DEEBF7", "#C6DBEF", "#9ECAE1", "#6BAED6", "#4292C6"),
    diverging      = c("#D7191C", "#FDAE61", "#FFFFBF", "#ABD9E9", "#2C7BB6")
  )
  pal <- palettes[[type]]
  if (!is.null(n)) pal <- colorRampPalette(pal)(n)
  pal
}

# Internal theme for consistent TimeTeller plots
#' @noRd
theme_tt <- function(base_size = 12) {
  theme_classic(base_size = base_size) %+replace%
    theme(
      plot.title        = element_text(size = base_size + 1, face = "bold", hjust = 0),
      plot.subtitle     = element_text(size = base_size - 1, color = "grey40"),
      axis.title        = element_text(size = base_size),
      axis.text         = element_text(size = base_size - 1),
      legend.title      = element_text(size = base_size - 1, face = "bold"),
      legend.text       = element_text(size = base_size - 2),
      legend.position   = "bottom",
      legend.background = element_blank(),
      strip.background  = element_rect(fill = "grey95", colour = NA),
      strip.text        = element_text(size = base_size - 1, face = "bold"),
      panel.grid.major  = element_line(colour = "grey90", linewidth = 0.3),
      plot.margin       = margin(10, 10, 10, 10)
    )
}

# Per-sample intergene normalisation helper (samples x genes matrix)
normalize_per_sample <- function(data_matrix) {
  normalized <- t(apply(data_matrix, 1, function(x) {
    (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
  }))
  colnames(normalized) <- colnames(data_matrix)
  as.data.frame(normalized)
}
