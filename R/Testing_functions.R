normalise_test_data <- function(object, exp_matrix, test_grouping_vars, test_group_1, test_group_2, test_group_3, test_replicate, test_time) {

  norm_method <- object[['Normalisation_choice']]

  if (norm_method == 'intergene') {
    normalised_mat <- apply(exp_matrix, 2, function(x) {(x - mean(x)) / sd(x)})
  }

  if (norm_method == 'clr') {
    # Centred log-ratio: log(x) minus mean of log(x) per sample
    normalised_mat <- apply(exp_matrix, 2, function(x) log(x) - mean(log(x)))
  }

  else if(norm_method == 'timecourse') {
    splitting_vec <- paste0(replace_na(test_group_1,''), replace_na(test_group_2,''), replace_na(test_group_3,''), replace_na(test_replicate,''))
    if(all(sapply(splitting_vec, function(x) x==''))) {stop("No groups specified for timecourse normalisation.")}
    split_df <- split(as.data.frame(t(exp_matrix)), splitting_vec)
    split_df_ts_norm <- lapply(split_df, function(x) as.data.frame(scale(x)))
    df_normalised <- unsplit(split_df_ts_norm, splitting_vec)
    normalised_mat <- t(df_normalised)
  }

  else if(norm_method == 'timecourse_matched') {
    group_metrics_list <- object[['Train']][['Gene_Per_Group_Info']]
    means <- group_metrics_list$Mean
    sds <- group_metrics_list$SD
    splitting_vec <- data.frame(Test_Group_1 = test_group_1, Test_Group_2 = test_group_2, Test_Group_3 = test_group_3, Test_Replicate = test_replicate) %>%
      tidyr::unite('test_timeseries_matched_name', all_of(test_grouping_vars), remove = FALSE, sep = '|') %>% dplyr::pull(test_timeseries_matched_name)
    split_df <- split(as.data.frame(t(exp_matrix[object$Metadata$Train$Genes_Used, ])), splitting_vec)
    object[['Check']] <- split_df
    df_names <- names(split_df)
    cat(paste('# Data found for the following groups: ', names(split_df)[names(split_df) %in% rownames(means)],'\n'))
    cat(paste('# Data NOT found for the following groups: ', names(split_df)[!(names(split_df) %in% rownames(means))],' ==> For these standard timecourse normalisaiton will be used\n'))
    split_df_ts_norm <- lapply(seq_along(split_df), function(i)
      as.data.frame(base::scale(split_df[[i]], center = (if (names(split_df)[i] %in% rownames(means)) as.vector(t(means[names(split_df)[i],])) else TRUE),
                                scale = (if (names(split_df)[i] %in% rownames(means)) as.vector(t(sds[names(split_df)[i],])) else TRUE))))
    df_normalised <- unsplit(split_df_ts_norm, splitting_vec)
    normalised_mat <- t(df_normalised)
  }

  else if(norm_method == 'combined') {
    group_metrics_list <- object[['Train']][['Gene_Per_Group_Info']]
    means <- group_metrics_list$Mean
    sds <- group_metrics_list$SD
    splitting_vec <- data.frame(Test_Group_1 = test_group_1, Test_Group_2 = test_group_2, Test_Group_3 = test_group_3, Test_Replicate = test_replicate) %>%
      tidyr::unite('test_timeseries_matched_name', all_of(test_grouping_vars), remove = FALSE, sep = '|') %>% dplyr::pull(test_timeseries_matched_name)
    split_df <- split(as.data.frame(t(base::scale(exp_matrix[object$Metadata$Train$Genes_Used, ]))), splitting_vec)
    object[['Check']] <- split_df
    df_names <- names(split_df)
    cat(paste('# Data found for the following groups: ', names(split_df)[names(split_df) %in% rownames(means)],'\n'))
    cat(paste('# Data NOT found for the following groups: ', names(split_df)[!(names(split_df) %in% rownames(means))],' ==> For these standard timecourse normalisaiton will be used\n'))
    split_df_ts_norm <- lapply(seq_along(split_df), function(i)
      as.data.frame(base::scale(split_df[[i]], center = (if (names(split_df)[i] %in% rownames(means)) as.vector(t(means[names(split_df)[i],])) else TRUE),
                                scale = (if (names(split_df)[i] %in% rownames(means)) as.vector(t(sds[names(split_df)[i],])) else TRUE))))
    df_normalised <- unsplit(split_df_ts_norm, splitting_vec)
    normalised_mat <- t(df_normalised)
  }

  return(normalised_mat)
}

add_test_data <- function(object, exp_matrix, test_grouping_vars, test_group_1, test_group_2, test_group_3, test_replicate, test_time, mat_normalised_test) {
  if(missing(test_group_1))   {test_group_1 <- as.character(rep(NA, dim(exp_matrix)[2]))}
  if(missing(test_group_2))   {test_group_2 <- as.character(rep(NA, dim(exp_matrix)[2]))}
  if(missing(test_group_3))   {test_group_3 <- as.character(rep(NA, dim(exp_matrix)[2]))}
  if(missing(test_replicate)) {test_replicate <- as.character(rep(NA, dim(exp_matrix)[2]))}
  if(missing(test_time))      {test_time <- as.character(rep(NA, dim(exp_matrix)[2]))}

  ll <- list(test_group_1,test_group_2,test_group_3,test_time,test_replicate)
  if (!all(sapply(ll,length)==length(ll[[1]]))) {stop("Supplied metadata is not all of equal length. Please check inputs")}
  if (!all(object[['Metadata']][['Train']][['Genes_Used']] %in% rownames(exp_matrix))) {stop("Not all training genes found in the supplied data. Please check inputs")}

  genes_used <- object[['Metadata']][['Train']][['Genes_Used']]
  if (mat_normalised_test) {
    object[['Test_Data']][['Full_Test_Data_Raw']] <- 'Not Provided'
    object[['Test_Data']][['Full_Test_Data']] <- as.matrix(exp_matrix)
  } else {
    object[['Test_Data']][['Full_Test_Data_Raw']] <- as.matrix(exp_matrix)
    object[['Test_Data']][['Full_Test_Data']] <- prepare_raw_counts(exp_matrix)
  }
  test_data <- as.matrix(exp_matrix[genes_used,])
  object[['Test_Data']][['Test_Exp_Data']] <- test_data
  object[['Test_Data']][['Normalised_Test_Exp_Data']] <- normalise_test_data(object = object, exp_matrix = test_data, test_grouping_vars = test_grouping_vars, test_group_1 = test_group_1, test_group_2 = test_group_2, test_group_3 = test_group_3, test_replicate = test_replicate, test_time = test_time)

  object[['Metadata']][['Test']][['Group_1']] <- test_group_1
  object[['Metadata']][['Test']][['Group_2']] <- test_group_2
  object[['Metadata']][['Test']][['Group_3']] <- test_group_3
  object[['Metadata']][['Test']][['Time']] <- test_time
  object[['Metadata']][['Test']][['Replicate']] <- test_replicate

  return(object)
}

## calc_test_likelis is now unified as calc_likelis in Training_functions.R
## with a backward-compatible wrapper.

## get_final_likelis_test, theta_calc_test, calc_flat_theta_contrib_test
## are now unified in Training_functions.R — wrappers preserved there for
## backward compatibility and call-site compatibility.

## second_peaks_fun_test is now unified as second_peaks_fun in
## Training_functions.R with a backward-compatible wrapper.

#' Shift predicted phases
#'
#' Shifts TimeTeller predicted times to more realistic time window IF there is a second likelihood peak located in that window
#'
#' @param object list containing TimeTeller training and testing models
#' @param minT lower bound for defined time window
#' @param maxT upper bound for defined time window
#'
#' @author Vadim Vasilyev
#'
#'
#' @return Returned is the updated object with corrected timing information
#' @export
#'

shift_outlier_times <- function(object, minT = 8, maxT = 20) {
  data <- object[['Test_Data']][['Results_df']] %>%
    dplyr::mutate(flagged_1st = if_else(time_1st_peak < minT | time_1st_peak > maxT, TRUE, FALSE), .after = time_1st_peak) %>%
    dplyr::mutate(flagged_2nd = if_else(time_2nd_peak < minT | time_2nd_peak > maxT, TRUE, FALSE), .after = time_2nd_peak) %>%
    dplyr::mutate(Corrected_Time = if_else(flagged_1st, if_else(flagged_2nd, time_1st_peak, time_2nd_peak), time_1st_peak),
                  Was_Shifted = if_else(flagged_1st, if_else(flagged_2nd, 0, 1), 0))

  object[['Test_Data']][['Results_df']] <- data

  message('Number of Time Predictions that were shifted to the defined window if applicable: ', sum(object[['Test_Data']][['Results_df']]$Was_Shifted, na.rm = TRUE))
  return(object)
}
