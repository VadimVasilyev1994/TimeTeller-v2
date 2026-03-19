make_data_object <- function(exp_matrix, genes, group_1, group_2, group_3, time, replicate, mat_normalised) {
  timeteller_list <- list()
  if (missing(group_1)) {
    message('Missing group information. Please provide high level group information or use intergene normalisation.')
    group_1 <- rep(NA, dim(exp_matrix)[2])
  }
  if (missing(group_2)) {group_2 <- rep(NA, dim(exp_matrix)[2])}
  if (missing(group_3)) {group_3 <- rep(NA, dim(exp_matrix)[2])}
  if (missing(replicate)) {replicate <- rep(NA, dim(exp_matrix)[2])}
  if (missing(time)) {stop("Model can't be trained without training times. Please provide time vector")}
  ll <- list(group_1,group_2,group_3,time,replicate)
  if (!all(sapply(ll,length)==length(ll[[1]]))) {stop("Supplied metadata is not all of equal length. Please check inputs")}
  if (!all(genes %in% rownames(exp_matrix))) {stop("Not all training genes found in the supplied data. Please check inputs")}
  if (mat_normalised) {
    timeteller_list[['Full_Original_Data_Raw']] <- 'Not Provided'
    timeteller_list[['Full_Original_Data']] <- as.matrix(exp_matrix)
  } else {
    timeteller_list[['Full_Original_Data_Raw']] <- as.matrix(exp_matrix)
    timeteller_list[['Full_Original_Data']] <- prepare_raw_counts(exp_matrix)
  }
  exp_df <- as.data.frame(t(exp_matrix[genes,]))
  colnames(exp_df) <- paste('Gene_', genes, sep='')
  full_df <- exp_df %>% dplyr::mutate(Group = as.character(group_1), Group_2 = as.character(group_2), Group_3 = as.character(group_3), Time = as.factor(time), Replicate = as.character(replicate))
  full_df <- full_df %>% dplyr::mutate_at(c('Group','Group_2','Group_3','Replicate'), ~tidyr::replace_na(.,""))
  cat('Please check the information provided:\n')
  cat(paste('# Number of genes used:',length(genes), ' (',paste(genes,collapse = ' '),')\n'))
  cat(paste('# Number of Group_1 levels:',length(unique(group_1)), ' (',paste(unique(group_1),collapse = ' '),')\n'))
  cat(paste('# Number of Group_2 levels:',length(unique(group_2)), ' (',paste(unique(group_2),collapse = ' '),')\n'))
  cat(paste('# Number of Group_3 levels:',length(unique(group_3)), ' (',paste(unique(group_3),collapse = ' '),')\n'))
  cat(paste('# Number of Replicate levels:',length(unique(replicate)), ' (',paste(unique(replicate),collapse = ' '),')\n'))
  cat(paste('# Number of unique time points Used for training:',length(unique(time)), ' (',paste(unique(time),collapse = ' '),')\n'))
  timeteller_list[['Train']][['Data']] <- full_df
  timeteller_list[['Train']][['Train_Exp_Matrix']] <- exp_matrix[genes,]
  timeteller_list[['Metadata']][['Train']][['Genes_Used']] <- genes
  timeteller_list[['Metadata']][['Train']][['Group_1']] <- group_1
  timeteller_list[['Metadata']][['Train']][['Group_2']] <- group_2
  timeteller_list[['Metadata']][['Train']][['Group_3']] <- group_3
  timeteller_list[['Metadata']][['Train']][['Time']] <- time
  timeteller_list[['Metadata']][['Train']][['Replicate']] <- replicate
  return(timeteller_list)
}

deal_with_replicates <- function(object, treat_independently = TRUE) {
  data <- object[['Train']][['Data']]
  if (treat_independently) {
    cat('Replicates (if there are any) will be treated as independent observations for the training model\n')
    object[['Train']][['Data']] <- data
    return(object)
  }

  cat('Replicates for the same groups and time will be averaged before training the model\n')
  gene_cols <- grep('^Gene_', colnames(data), value = TRUE)

  # Count replicates per group-time combination
  counts_df <- as.data.frame(table(data$Group, data$Group_2, data$Group_3, data$Time)) %>%
    filter(Freq > 1)

  if (nrow(counts_df) > 0) {
    # Average gene expression across replicates within each group-time combination
    averaged <- data %>%
      group_by(Group, Group_2, Group_3, Time) %>%
      summarise(across(all_of(gene_cols), mean), .groups = 'drop') %>%
      mutate(Replicate = '')
    object[['Train']][['Data']] <- averaged
    cat('There were', nrow(counts_df), 'groups that were averaged\n')
  } else {
    cat('No replicates were found so no averaging was done. Please check inputs if replicates expected\n')
    object[['Train']][['Data']] <- data
  }
  return(object)
}

allocate_timeseries_ind <- function(object, combine_for_norm = FALSE) {
  data <- object[['Train']][['Data']]
  data$Time_mod24 <- as.numeric(as.character(data$Time)) %% 24

  if (combine_for_norm == TRUE) {
    data <- data %>% dplyr::group_by(Group, Group_2, Group_3) %>% dplyr::mutate(timeseries_name = factor(paste0(as.character(Group), as.character(Group_2), as.character(Group_3)))) %>%
      dplyr::mutate(timeseries_ind = as.numeric(timeseries_name))
  } else {
    data <- data %>% dplyr::group_by(Group, Group_2, Group_3, Replicate) %>% dplyr::mutate(timeseries_name = factor(paste0(as.character(Group), as.character(Group_2), as.character(Group_3), as.character(Replicate)))) %>%
      dplyr::mutate(timeseries_ind = as.numeric(timeseries_name))
  }

  object[['Train']][['Data']] <- data
  return(object)
}

normalised_df <- function(object, method = 'intergene', grouping_vars = c('Group')) {
  data <- object[['Train']][['Data']]
  genes <- colnames(data %>% ungroup() %>% dplyr::select(starts_with('Gene_')))
  object[['Normalisation_choice']] <- method

  # Helper: shared postprocessing for all normalisation methods
  # Adds time_group factor, arranges by time, and stores normalised data + matrix
  store_normalised <- function(object, norm_df) {
    norm_df <- norm_df %>%
      mutate(time_group = factor(
        paste('Time_', Time_mod24, sep = ''),
        levels = paste('Time_', sort(unique(norm_df$Time_mod24)), sep = '')
      )) %>%
      arrange(Time_mod24)
    object[['Train']][['Normalised_Data']] <- norm_df
    object[['Train']][['Normalised_Train_Exp_Data']] <- norm_df %>%
      ungroup() %>% dplyr::select(matches("normalised")) %>% as.matrix() %>% t()
    return(object)
  }

  if(method == 'timecourse') {
    cat('Using timecourse normalisation\n')
    norm_df <- data %>% group_by(timeseries_ind) %>%
      mutate(across(all_of(genes), .fns = list(normalised = ~(.x - mean(.x))/sd(.x)))) %>%
      rename_at(vars(contains("_normalised")), list(~paste("normalised", gsub("_normalised", "", .), sep = "_"))) %>%
      ungroup()
    return(store_normalised(object, norm_df))
  }

  else if(method == 'intergene'){
    cat('Using intergene normalisation\n')
    norm_df <- data %>% ungroup() %>% rowwise() %>%
      mutate(mean = mean(c_across(all_of(genes))), sd = sd(c_across(all_of(genes)))) %>%
      ungroup() %>%
      mutate(across(all_of(genes), .fns = list(normalised = ~(.x - mean)/sd))) %>%
      rename_at(vars(contains("_normalised")), list(~paste("normalised", gsub("_normalised", "", .), sep = "_")))
    return(store_normalised(object, norm_df))
  }

  else if(method == 'clr'){
    cat('Using clr normalisation\n')
    norm_df <- data %>% ungroup() %>% rowwise() %>%
      mutate(geom_mean = exp(mean(log(c_across(all_of(genes)))))) %>%
      ungroup() %>%
      mutate(across(all_of(genes), .fns = list(normalised = ~ (log(.x) - log(geom_mean))))) %>%
      rename_at(vars(contains("_normalised")), list(~paste("normalised", gsub("_normalised", "", .), sep = "_")))
    return(store_normalised(object, norm_df))
  }

  else if(method == 'timecourse_matched'){
    cat('Using timecourse-matched normalisation. Pay attention to the group names used when testing\n')
    norm_df <- data %>% group_by(timeseries_ind) %>%
      mutate(across(all_of(genes), .fns = list(normalised = ~(.x - mean(.x))/sd(.x)))) %>%
      rename_at(vars(contains("_normalised")), list(~paste("normalised", gsub("_normalised", "", .), sep = "_")))
    object <- store_normalised(object, norm_df)

    # Store per-group mean/SD for matched normalisation of test data
    group_statistics <- data %>% group_by(across(all_of(grouping_vars))) %>% dplyr::select(all_of(genes)) %>%
      summarise_all(list(Mean = mean, SD = sd)) %>%
      rename_at(vars(contains("Gene_")), list(~paste("", gsub("Gene_", "", .), sep = ""))) %>%
      tidyr::unite('timeseries_matched_name', all_of(grouping_vars), remove = FALSE, sep = '|')
    group_statistics <- as.data.frame(group_statistics)
    rownames(group_statistics) <- group_statistics$timeseries_matched_name
    group_statistics$timeseries_matched_name <- NULL
    group_metrics_list <- sapply(c("Mean", "SD"), function(x) group_statistics[endsWith(names(group_statistics), x)], simplify = FALSE)
    object[['Train']][['Gene_Per_Group_Info']] <- group_metrics_list
    return(object)
  }

  else if(method == 'combined'){
    cat('Using intergene then timecourse-matched normalisation. Pay attention to the group names used when testing\n')
    # Step 1: intergene normalisation
    intergene_norm_df <- data %>% ungroup() %>% rowwise() %>%
      mutate(mean = mean(c_across(all_of(genes))), sd = sd(c_across(all_of(genes)))) %>%
      ungroup() %>%
      mutate(across(all_of(genes), .fns = list(normalised = ~(.x - mean)/sd))) %>%
      rename_at(vars(contains("_normalised")), list(~paste("Intergene", gsub("_normalised", "", .), sep = "_")))

    new_genes <- colnames(intergene_norm_df %>% ungroup() %>% dplyr::select(starts_with('Intergene_')))

    # Step 2: timecourse normalisation on intergene-normalised values
    norm_df <- intergene_norm_df %>% group_by(timeseries_ind) %>%
      mutate(across(all_of(new_genes), .fns = list(normalised = ~(.x - mean(.x))/sd(.x)))) %>%
      rename_at(vars(contains("_normalised")), list(~paste("normalised", gsub("_normalised", "", .), sep = "_")))
    object <- store_normalised(object, norm_df)

    # Store per-group mean/SD for matched normalisation of test data
    group_statistics <- intergene_norm_df %>% group_by(across(all_of(grouping_vars))) %>% dplyr::select(all_of(new_genes)) %>%
      summarise_all(list(Mean = mean, SD = sd)) %>%
      rename_at(vars(contains("Intergene_Gene_")), list(~paste("", gsub("Intergene_Gene_", "", .), sep = ""))) %>%
      tidyr::unite('timeseries_matched_name', all_of(grouping_vars), remove = FALSE, sep = '|')
    group_statistics <- as.data.frame(group_statistics)
    rownames(group_statistics) <- group_statistics$timeseries_matched_name
    group_statistics$timeseries_matched_name <- NULL
    group_metrics_list <- sapply(c("Mean", "SD"), function(x) group_statistics[endsWith(names(group_statistics), x)], simplify = FALSE)
    object[['Train']][['Gene_Per_Group_Info']] <- group_metrics_list
    return(object)
  }
}

get_projections <- function(object, num_PC) {
  data <- object[['Train']][['Normalised_Data']]
  object[['PC_Num']] <- num_PC
  group_list <- data %>% dplyr::ungroup() %>% dplyr::select(starts_with('normalised_')) %>% split.data.frame(data$time_group)
  group_list <- lapply(group_list, as.matrix)
  group_list <- lapply(group_list, base::t)
  object[['Projections']][['Training_Exp_Per_Time_Point']] <- group_list
  pr_comp_list <- purrr::map(group_list, as_mapper(~ base::t(svd(.x)$u[,1:num_PC])))
  # pr_comp_list <- purrr::map(group_list, as_mapper(~ base::t(svd(sweep(.x, 1, apply(.x, 1, mean), '-'))$u[,1:num_PC])))
  var_explained_list <- purrr::map(group_list, as_mapper(~ cumsum(svd(.x)$d^2 / sum(svd(.x)$d^2) * 100)))
  # var_explained_list <- purrr::map(group_list, as_mapper(~ cumsum(svd(sweep(.x, 1, apply(.x, 1, mean), '-'))$d^2 / sum(svd(.x)$d^2) * 100)))
  object[['Projections']][['SVD_Per_Time_Point']] <- pr_comp_list
  object[['Projections']][['SVD_Per_Time_Point_Var_Explained']] <- var_explained_list
  all_projections <- list()
  for (j in seq_along(names(pr_comp_list))) {
    local_proj_list <- list()
    for (i in seq_along(names(group_list))) {
      local_proj_list[[names(group_list)[i]]] <- pr_comp_list[[j]] %*% group_list[[i]]
    }
    pc_vec <- paste('PC', seq(1:num_PC), sep = '')
    df <-  as.data.frame(t(do.call(cbind, local_proj_list)))
    colnames(df) <- pc_vec
    df <- df %>% dplyr::mutate(sample_times = data$Time_mod24, Group = data$Group, Group_2 = data$Group_2, Group_3 = data$Group_3)
    all_projections[[names(pr_comp_list)[j]]] <- df
  }
  object[['Projections']][['All_Projections']] <- all_projections
  return(object)
}

get_mvn_original_data <- function(object, cov_estimate = 'normal', alpha_par = 0.7) {
  if (cov_estimate == 'normal') {cat('Using standard Covariance and Location estimation\n')}
  if (cov_estimate == 'robust') {cat('Using Robust Covariance and Location Estimation (Fast MCD)\n')}
  object[['Projections']][['Cov_Method']] <- cov_estimate
  num_PC <- object[['PC_Num']]
  projection_list <- object[['Projections']][['All_Projections']]
  likeli_array_original <- base::array(data = NA, dim = c(num_PC + num_PC^2, length(names(projection_list)), length(names(projection_list))))
  dimnames(likeli_array_original) <- list(c(paste0('mean_PC',1:num_PC),paste0('sigma',1:num_PC^2)), names(projection_list), names(projection_list))
  for (i in 1:length(names(projection_list))) {
    curr_projection <- projection_list[[i]]
    data.split <- split(curr_projection[, 1:num_PC], curr_projection$sample_times)
    if (cov_estimate == 'normal') {
      data.mvnorm <- lapply(data.split, function(x) Rfast::mvnorm.mle(as.matrix(x)))
      mat <- purrr::map(data.mvnorm, extract_elements_func)
      mat_final <- base::matrix(unlist(mat), nrow = num_PC + num_PC^2)
      likeli_array_original[,,i] <- mat_final
    }
    if (cov_estimate == 'robust') {
      check_suggested_pkg("rrcov")
      data.mvnorm <- lapply(data.split, function(x) rrcov::CovMcd(as.matrix(x), alpha = alpha_par, nsamp = 'deterministic'))
      mat <- purrr::map(data.mvnorm, extract_elements_func_robust)
      mat_final <- base::matrix(unlist(mat), nrow = num_PC + num_PC^2)
      likeli_array_original[,,i] <- mat_final
    }
  }

  object[['Projections']][['Fitted_MVN_Original']] <- likeli_array_original
  return(object)
}

get_mvn_interpolated <- function(object, num_interp_points = 144, interp_method = 'perpchip', cov_path = 'spline') {
  likelihood_array <- object[['Projections']][['Fitted_MVN_Original']]
  numPC <- object[['PC_Num']]
  dims_array <- dim(likelihood_array)
  interpolated_array <- base::array(data = NA, dim = c(dims_array[1], num_interp_points, dims_array[3]))
  times <- as.numeric(sub(".*Time_", "", dimnames(likelihood_array)[[2]]))
  times_periodic <- c(times, times[1]+24)
  object[['Metadata']][['Train']][['min_T_mod24']] <- min(times_periodic)
  object[['Metadata']][['Train']][['max_T_mod24']] <- max(times_periodic)
  if (cov_path == 'spline') {
    for (j in 1:dims_array[3]) {
      for (i in 1:dims_array[1]) {
        ts_to_interpolate <- likelihood_array[i,,j]
        ts_to_interpolate_periodic <- c(ts_to_interpolate, ts_to_interpolate[1])

        if (interp_method == 'perpchip') {
          interpolated_array[i,,j] <- perpchip(times_periodic, ts_to_interpolate_periodic, seq(min(times_periodic), max(times_periodic), length.out = num_interp_points))
        }
        if (interp_method == 'standard') {
          interpolated_array[i,,j] <- stats::spline(times_periodic, ts_to_interpolate_periodic, n = num_interp_points, method = "periodic")$y
        }

      }
    }
    dimnames(interpolated_array) <- list(dimnames(likelihood_array)[[1]], paste0('Time_',round(seq(min(times_periodic), max(times_periodic), length.out = num_interp_points),2)), dimnames(likelihood_array)[[3]])
    object[['Projections']][['Fitted_MVN_Interpolated']] <- interpolated_array
    return(object)
  }


  if (cov_path == 'fisherrao') {
    for (j in 1:dims_array[3]) {

      num_points <- num_interp_points + dims_array[3]
      times_num <- as.numeric(substring(colnames(likelihood_array), first = 6))
      interpolation_points <- c()
      for (i in 1:length(times_num)) {
        interpolation_points[i] <- length(seq(0, 1, length.out = floor(((abs(circular_difference(times_num[i],times_num[ifelse(i+1>length(times_num),1,i+1)])))*60) / (24*60/num_points))))
      }
      interp_res <- list()
      for (i in 1:length(times_num)) {
        cov_mat_1 <- matrix(likelihood_array[(numPC+1):(numPC+numPC^2),i,j], nrow = 3)
        cov_mat_2 <- matrix(likelihood_array[(numPC+1):(numPC+numPC^2),ifelse(i+1>length(times_num),1,i+1),j], nrow = 3)
        interp_points <- seq(0, 1, length.out = interpolation_points[i])

        mat_list <- matrix(NA, nrow = numPC^2, ncol = length(interp_points))
        for (k in 1:length(interp_points)) {
          first_term <- pracma::sqrtm(cov_mat_1)$B
          second_term <- solve(pracma::sqrtm(cov_mat_1)$B) %*% cov_mat_2 %*% solve(pracma::sqrtm(cov_mat_1)$B)
          second_term <- expm::expm(expm::logm(second_term) * interp_points[k])
          third_term <- pracma::sqrtm(cov_mat_1)$B

          interp_mat <- first_term %*% second_term %*% third_term
          mat_list[,k] <- as.vector(interp_mat)
        }
        mat_names <- paste0('Time_', round(seq(times_num[i], ifelse(i+1>length(times_num),times_num[1]+24,times_num[i+1]), length.out = interpolation_points[i]),2))
        colnames(mat_list) <- mat_names
        if (i>1) {mat_list <- mat_list[,-1]}
        interp_res[[i]] <- mat_list
      }

      check <- do.call(cbind,interp_res)
      check <- check[, !duplicated(colnames(check))]
      interpolated_array <- interpolated_array[ , 1:dim(check)[2], ]
      interpolated_array[(numPC+1):(numPC+numPC^2), ,j] <- check


      for (i in 1:numPC) {
        ts_to_interpolate <- likelihood_array[i,,j]
        ts_to_interpolate_periodic <- c(ts_to_interpolate, ts_to_interpolate[1])

        if (interp_method == 'perpchip') {
          interpolated_array[i,,j] <- perpchip(times_periodic, ts_to_interpolate_periodic, seq(min(times_periodic), max(times_periodic), length.out = dim(check)[2]))
        }
        if (interp_method == 'standard') {
          interpolated_array[i,,j] <- stats::spline(times_periodic, ts_to_interpolate_periodic, n = dim(check)[2], method = "periodic")$y
        }

      }

      dimnames(interpolated_array)[[2]] <- colnames(check)

    }
    dimnames(interpolated_array) <- list(dimnames(likelihood_array)[[1]], dimnames(interpolated_array)[[2]], dimnames(likelihood_array)[[3]])


    object[['Projections']][['Fitted_MVN_Interpolated']] <- interpolated_array
    return(object)
  }

  if (cov_path == 'wasserstein') {
    # Wasserstein (Bures) geodesic for covariance interpolation.
    # For adjacent knots K0, K1, the geodesic at parameter t is:
    #   K(t) = (1-t)^2*K(0) + t^2*K(1) + t(1-t)[(K(0)*K(1))^{1/2} + (K(1)*K(0))^{1/2}]
    # This guarantees PD at every interpolation point (Mallasto et al. 2020, Eq. 9).
    # Means are interpolated via splines as in the other paths.

    for (j in 1:dims_array[3]) {
      num_points <- num_interp_points + dims_array[3]
      times_num <- as.numeric(substring(colnames(likelihood_array), first = 6))

      # Compute number of interpolation points per segment (proportional to time gap)
      interpolation_points <- c()
      for (i in 1:length(times_num)) {
        next_i <- ifelse(i + 1 > length(times_num), 1, i + 1)
        gap <- abs(circular_difference(times_num[i], times_num[next_i]))
        interpolation_points[i] <- length(seq(0, 1, length.out = floor((gap * 60) / (24 * 60 / num_points))))
      }

      interp_res <- list()
      for (i in 1:length(times_num)) {
        next_i <- ifelse(i + 1 > length(times_num), 1, i + 1)
        K0 <- matrix(likelihood_array[(numPC + 1):(numPC + numPC^2), i, j], nrow = numPC)
        K1 <- matrix(likelihood_array[(numPC + 1):(numPC + numPC^2), next_i, j], nrow = numPC)
        interp_t <- seq(0, 1, length.out = interpolation_points[i])

        # Pre-compute the cross terms (K0 K1)^{1/2} and (K1 K0)^{1/2}
        # Use: (AB)^{1/2} = A^{1/2} (A^{1/2} B A^{1/2})^{1/2} A^{-1/2}
        K0_sqrt <- expm::sqrtm(K0)
        K0_inv_sqrt <- solve(K0_sqrt)
        inner <- K0_sqrt %*% K1 %*% K0_sqrt
        inner_sqrt <- expm::sqrtm(inner)
        K0K1_sqrt <- K0_sqrt %*% inner_sqrt %*% K0_inv_sqrt  # (K0 K1)^{1/2}
        K1K0_sqrt <- t(K0K1_sqrt)  # (K1 K0)^{1/2} = [(K0 K1)^{1/2}]^T for symmetric K0, K1

        mat_list <- matrix(NA, nrow = numPC^2, ncol = length(interp_t))
        for (k in seq_along(interp_t)) {
          tt <- interp_t[k]
          Kt <- (1 - tt)^2 * K0 + tt^2 * K1 + tt * (1 - tt) * (K0K1_sqrt + K1K0_sqrt)
          mat_list[, k] <- as.vector(Kt)
        }

        t_start <- times_num[i]
        t_end <- ifelse(next_i == 1, times_num[1] + 24, times_num[next_i])
        mat_names <- paste0('Time_', round(seq(t_start, t_end, length.out = interpolation_points[i]), 2))
        colnames(mat_list) <- mat_names
        if (i > 1) mat_list <- mat_list[, -1]
        interp_res[[i]] <- mat_list
      }

      check <- do.call(cbind, interp_res)
      check <- check[, !duplicated(colnames(check))]
      interpolated_array <- interpolated_array[, 1:dim(check)[2], ]
      interpolated_array[(numPC + 1):(numPC + numPC^2), , j] <- check

      # Interpolate means via splines (same as other paths)
      for (i in 1:numPC) {
        ts_to_interpolate <- likelihood_array[i, , j]
        ts_to_interpolate_periodic <- c(ts_to_interpolate, ts_to_interpolate[1])

        if (interp_method == 'perpchip') {
          interpolated_array[i, , j] <- perpchip(times_periodic, ts_to_interpolate_periodic, seq(min(times_periodic), max(times_periodic), length.out = dim(check)[2]))
        }
        if (interp_method == 'standard') {
          interpolated_array[i, , j] <- stats::spline(times_periodic, ts_to_interpolate_periodic, n = dim(check)[2], method = "periodic")$y
        }
      }

      dimnames(interpolated_array)[[2]] <- colnames(check)
    }

    dimnames(interpolated_array) <- list(dimnames(likelihood_array)[[1]], dimnames(interpolated_array)[[2]], dimnames(likelihood_array)[[3]])
    object[['Projections']][['Fitted_MVN_Interpolated']] <- interpolated_array
    return(object)
  }
}

#' Compute log-likelihood array across all local projections
#'
#' For each local SVD projection (one per training time point), project the
#' expression data and evaluate the multivariate normal density at every
#' interpolation point. The result is a 3D array:
#'   (interpolation_points x samples x local_projections)
#'
#' @param object TimeTeller list object
#' @param mode either 'train' or 'test'
#' @return Updated object with log-likelihood array (and test projections if mode='test')
#' @keywords internal
calc_likelis <- function(object, mode = 'train', diagnose_pd = FALSE) {
  svd_data <- object[['Projections']][['SVD_Per_Time_Point']]
  fitted_mvn_data <- object[['Projections']][['Fitted_MVN_Interpolated']]
  num_PC <- object[['PC_Num']]

  # Select expression data based on mode
  if (mode == 'train') {
    exp_data <- object[['Train']][['Normalised_Train_Exp_Data']]
  } else {
    exp_data <- object[['Test_Data']][['Normalised_Test_Exp_Data']]
  }

  n_interp <- dim(fitted_mvn_data)[2]
  n_samples <- dim(exp_data)[2]
  n_projections <- dim(fitted_mvn_data)[3]
  likelihood_array <- base::array(data = NA, dim = c(n_interp, n_samples, n_projections))

  # Store test projections for downstream visualisation
  if (mode == 'test') projections_list <- list()

  # Track nearPD corrections if diagnostics requested
  if (diagnose_pd) {
    pd_total <- 0L
    pd_corrected <- 0L
  }

  for (i in seq_along(names(svd_data))) {
    project_exp_mat <- svd_data[[i]] %*% exp_data
    if (mode == 'test') projections_list[[i]] <- project_exp_mat

    # --- Pre-compute PD-corrected covariance matrices for this projection ---
    sigma_list <- vector("list", n_interp)
    mean_list <- vector("list", n_interp)
    for (j in seq_len(n_interp)) {
      curr_sigma <- matrix(fitted_mvn_data[(num_PC + 1):(num_PC + num_PC^2), j, i], nrow = num_PC)
      curr_eig <- eigen(curr_sigma, symmetric = TRUE, only.values = TRUE)$values
      if (diagnose_pd) pd_total <- pd_total + 1L
      if (any(curr_eig < 0)) {
        if (diagnose_pd) pd_corrected <- pd_corrected + 1L
        sigma_list[[j]] <- as.matrix(nearPD(curr_sigma, base.matrix = TRUE, ensureSymmetry = TRUE,
                                            eig.tol = 1e-05, conv.tol = 1e-06, posd.tol = 1e-07)$mat)
      } else {
        sigma_list[[j]] <- curr_sigma
      }
      mean_list[[j]] <- fitted_mvn_data[1:num_PC, j, i]
    }

    # --- Vectorised density evaluation: all samples at each interp point ---
    samples_mat <- t(project_exp_mat)
    for (j in seq_len(n_interp)) {
      likelihood_array[j, , i] <- mvtnorm::dmvnorm(samples_mat,
                                                     mean = mean_list[[j]],
                                                     sigma = sigma_list[[j]],
                                                     checkSymmetry = FALSE)
    }

    cat("\rFinished", i, "of", length(names(svd_data)), "\n")
  }

  # Report nearPD diagnostics
  if (diagnose_pd) {
    pct <- round(100 * pd_corrected / pd_total, 2)
    message("nearPD diagnostic: ", pd_corrected, " / ", pd_total,
            " covariance matrices required correction (", pct, "%)")
  }

  # Store results in appropriate slots
  if (mode == 'train') {
    object[['Train_Data']][['Train_Likelihood_Array']] <- log(likelihood_array)
  } else {
    names(projections_list) <- names(svd_data)
    object[['Test_Data']][['Test_Projections']] <- projections_list
    object[['Test_Data']][['Test_Likelihood_Array']] <- log(likelihood_array)
  }
  return(object)
}

# Backward-compatible wrappers
calc_train_likelis <- function(object, diagnose_pd = FALSE) {
  calc_likelis(object, mode = 'train', diagnose_pd = diagnose_pd)
}
calc_test_likelis <- function(object, diagnose_pd = FALSE) {
  calc_likelis(object, mode = 'test', diagnose_pd = diagnose_pd)
}

## ---------------------------------------------------------------------------
## Unified internal functions (train/test share the same computation)
## ---------------------------------------------------------------------------

#' Apply log threshold and compute averaged likelihoods
#' @param object TimeTeller list object
#' @param log_thresh log threshold value
#' @param mode either 'train' or 'test'
#' @return Updated object with thresholded and averaged likelihoods
#' @keywords internal
get_final_likelis <- function(object, log_thresh, mode = 'train') {
  # Resolve slot names based on mode
  if (mode == 'train') {
    data_slot <- 'Train_Data'
    likelis_key <- 'Train_Likelihood_Array'
    thresh_key <- 'LogThresh_Train'
    post_thresh_key <- 'Likelis_Post_Thresh_Train'
    avg_key <- 'Averaged_Likelis_Post_Thresh_Train'
    max_key <- 'Max_Likelis_Train'
  } else {
    data_slot <- 'Test_Data'
    likelis_key <- 'Test_Likelihood_Array'
    thresh_key <- 'LogThresh_Test'
    post_thresh_key <- 'Likelis_Post_Thresh_Test'
    avg_key <- 'Averaged_Likelis_Post_Thresh_Test'
    max_key <- 'Max_Likelis_Test'
  }

  likelis_array <- object[[data_slot]][[likelis_key]]
  object[[data_slot]][[thresh_key]] <- log_thresh

  # Apply floor threshold to each local projection slice
  for (i in seq_len(dim(likelis_array)[3])) {
    likelis_array[,,i] <- pmax(likelis_array[,,i], log_thresh)
  }

  object[[data_slot]][[post_thresh_key]] <- likelis_array
  # Average across local projections
  final_averaged_likelis <- apply(likelis_array, c(1,2), mean)
  object[[data_slot]][[avg_key]] <- final_averaged_likelis
  object[[data_slot]][[max_key]] <- apply(final_averaged_likelis, 2, max)
  return(object)
}

# Backward-compatible wrappers
get_final_likelis_train <- function(object, log_thresh) {
  get_final_likelis(object, log_thresh, mode = 'train')
}
get_final_likelis_test <- function(object, log_thresh) {
  get_final_likelis(object, log_thresh, mode = 'test')
}

#' Compute Theta (clock dysfunction metric) for each sample
#'
#' Theta measures the proportion of the likelihood curve that exceeds a
#' cosine envelope, capturing how peaked vs diffuse the likelihood is.
#' A well-functioning clock produces a sharply peaked likelihood (low Theta);
#' a disrupted clock produces a broad/flat likelihood (high Theta).
#'
#' @param object TimeTeller list object
#' @param mode either 'train' or 'test'
#' @param epsilon envelope amplitude parameter (default 0.4)
#' @param eta envelope scaling parameter (default 0.35)
#' @return Updated object with Theta values
#' @keywords internal
theta_calc <- function(object, mode = 'train', epsilon = NULL, eta = NULL) {
  if (mode == 'train') {
    averaged_likelis_rescaled <- t(object[['Train_Data']][['Averaged_Likelis_Post_Thresh_Train']])
  } else {
    epsilon <- object[['Train_Data']][['epsilon']]
    eta <- object[['Train_Data']][['eta']]
    averaged_likelis_rescaled <- t(object[['Test_Data']][['Averaged_Likelis_Post_Thresh_Test']])
  }

  num_samples <- dim(averaged_likelis_rescaled)[1]
  num_points <- dim(averaged_likelis_rescaled)[2]
  n_eval <- 500L  # spline evaluation points (sufficient for theta precision)
  eval_seq <- seq(1, num_points, length.out = n_eval)
  thetas <- numeric(num_samples)

  for (i in seq_len(num_samples)) {
    curr_sample <- exp(averaged_likelis_rescaled[i,])
    curr_lrf_curve <- curr_sample / max(curr_sample)
    peak_pos <- which.max(curr_lrf_curve)
    curr_curve <- suppressWarnings(
      eta * (1 + epsilon + cos(2 * pi * ((1:num_points) / num_points - peak_pos / num_points)))
    )

    lrf_curve_spline <- stats::predict(
      periodicSpline(1:num_points, curr_lrf_curve, period = num_points), eval_seq
    )
    curve_spline <- stats::predict(
      periodicSpline(1:num_points, curr_curve, period = num_points), eval_seq
    )

    thetas[i] <- sum(lrf_curve_spline$y > curve_spline$y) / n_eval
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
theta_calc_train <- function(object, epsilon = 0.4, eta = 0.35) {
  theta_calc(object, mode = 'train', epsilon = epsilon, eta = eta)
}
theta_calc_test <- function(object) {
  theta_calc(object, mode = 'test')
}

#' Compute flat region contribution to Theta
#'
#' Quantifies what proportion of Theta is attributable to flat (saturated)
#' regions in the likelihood curve — i.e. regions clamped at the log threshold.
#' High flat contribution suggests Theta is inflated by the threshold rather
#' than genuine clock dysfunction.
#'
#' @param object TimeTeller list object
#' @param mode either 'train' or 'test'
#' @return Updated object with flat contribution data frame
#' @keywords internal
calc_flat_theta_contrib <- function(object, mode = 'train') {
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
  n_eval <- 500L
  eval_seq <- seq(1, num_points, length.out = n_eval)
  flat_contribution <- numeric(num_samples)

  for (i in seq_len(num_samples)) {
    ind_ts <- exp(averaged_likelis_rescaled[i,])
    ind_ts <- shift_ts(ind_ts, round(num_points / 2))
    ind_lrf_curve <- ind_ts / max(ind_ts)

    peak_pos <- which.max(ind_lrf_curve)
    curr_curve <- suppressWarnings(
      eta * (1 + epsilon + cos(2 * pi * ((1:num_points) / num_points - peak_pos / num_points)))
    )
    lrf_curve_spline <- stats::predict(
      periodicSpline(1:num_points, ind_lrf_curve, period = num_points), eval_seq
    )
    curve_spline <- stats::predict(
      periodicSpline(1:num_points, curr_curve, period = num_points), eval_seq
    )

    theta_index <- which(lrf_curve_spline$y > curve_spline$y)
    flat_regions <- which(abs(diff(lrf_curve_spline$y)) < 1e-4)
    flat_contribution[i] <- sum(theta_index %in% flat_regions) / length(theta_index) * 100
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
calc_flat_theta_contrib_train <- function(object) {
  calc_flat_theta_contrib(object, mode = 'train')
}
calc_flat_theta_contrib_test <- function(object) {
  calc_flat_theta_contrib(object, mode = 'test')
}

#' Extract peak information and build results data frame
#'
#' Finds peaks in the averaged likelihood curve, computes prediction statistics
#' (weighted mean time, circular mean/SD across local projections), and assembles
#' the full results data frame with metadata.
#'
#' @param object TimeTeller list object
#' @param mode either 'train' or 'test'
#' @param minpeakheight,minpeakdistance,nups,ndowns,threshold,npeaks passed to pracma::findpeaks
#' @return Updated object with Results_df in the appropriate data slot
#' @keywords internal
second_peaks_fun <- function(object, mode = 'train', minpeakheight = -Inf, minpeakdistance = 1,
                             nups = 1, ndowns = 0, threshold = 0, npeaks = 2) {

  # Resolve slot names
  if (mode == 'train') {
    data_slot <- 'Train_Data'
    likelis_key <- 'Train_Likelihood_Array'
    thresh_key <- 'LogThresh_Train'
    avg_key <- 'Averaged_Likelis_Post_Thresh_Train'
  } else {
    data_slot <- 'Test_Data'
    likelis_key <- 'Test_Likelihood_Array'
    thresh_key <- 'LogThresh_Test'
    avg_key <- 'Averaged_Likelis_Post_Thresh_Test'
  }

  likelis_array <- object[[data_slot]][[likelis_key]]
  logthresh <- object[[data_slot]][[thresh_key]]
  flat_contributions_df <- object[[data_slot]][['Flat_Contrib_to_Theta_df']]
  averaged_likelis <- object[[data_slot]][[avg_key]]

  npoints <- dim(averaged_likelis)[1]

  # --- Find peaks in averaged likelihood curves ---
  peaks_list <- apply(averaged_likelis, 2, pracma::findpeaks, nups = nups, ndowns = ndowns,
                      sortstr = TRUE, minpeakheight = minpeakheight, minpeakdistance = minpeakdistance,
                      threshold = threshold, npeaks = npeaks, simplify = FALSE)

  npeaks_list <- map(peaks_list, ~nrow(.x))
  npeaks_list[sapply(npeaks_list, is.null)] <- NA
  npeaks_vec <- unlist(npeaks_list)

  log_values <- map(peaks_list, ~.x[,1])
  log_values[sapply(log_values, is.null)] <- NA

  times <- map(peaks_list, ~.x[,2])
  times[sapply(times, is.null)] <- NA
  times_hours <- map(times, ~.x / npoints * 24)

  # Convert peak indices to hours, offset by the minimum training time
  min_t <- object[['Metadata']][['Train']][['min_T_mod24']]
  times_1st_peak <- (map_dbl(times_hours, ~.x[1]) + min_t) %% 24
  times_2nd_peak <- (map_dbl(times_hours, ~.x[2]) + min_t) %% 24
  log_likeli_1st_peak <- map_dbl(log_values, ~.x[1])
  log_likeli_2nd_peak <- map_dbl(log_values, ~.x[2])

  peaks_df <- data.frame(
    npeaks = npeaks_vec,
    time_1st_peak = round(times_1st_peak, 2),
    time_2nd_peak = round(times_2nd_peak, 2),
    max_1st_peak = round(log_likeli_1st_peak, 3),
    max_2nd_peak = round(log_likeli_2nd_peak, 3)
  )
  peaks_df <- peaks_df %>%
    mutate(times_diff = time_1st_peak - time_2nd_peak,
           max_diff = max_1st_peak - max_2nd_peak)

  # --- Flat likelihood metrics ---
  percent_flat <- function(likelis, thresh) {
    round(sum(likelis == thresh) / npoints * 100, 1)
  }
  peaks_df$PercFlat <- apply(averaged_likelis, 2, percent_flat, logthresh)
  peaks_df$FlatContrib <- flat_contributions_df$flat_contributions
  peaks_df$Theta <- flat_contributions_df$thetas

  # --- Per-projection peak statistics (before thresholding) ---
  local_max_vals <- apply(likelis_array, c(2,3), max)
  peaks_df <- peaks_df %>%
    mutate(local_smallest_peak = round(apply(local_max_vals, 1, min), 3),
           local_biggest_peak = round(apply(local_max_vals, 1, max), 3),
           local_mean = round(apply(local_max_vals, 1, mean), 3),
           local_sd = round(apply(local_max_vals, 1, sd), 3))

  local_times <- (apply(likelis_array, c(2,3), which.max) / npoints * 24 + min_t) %% 24
  time_weights <- t(apply(local_max_vals, 1, function(x) (x - min(x)) / (max(x) - min(x))))

  # Weighted average time prediction using max likelihood values
  weighted_times_vec <- apply(sweep(local_times * time_weights, 1, apply(time_weights, 1, sum), '/'), 1, sum)

  # Circular mean and deviation of local time predictions
  local_times_rad <- local_times / 24 * 2 * pi
  mean_local_times <- (24 + suppressWarnings(apply(local_times_rad, 1, circular::mean.circular)) * 24 / 2 / pi) %% 24
  sd_local_times <- suppressWarnings(apply(local_times_rad, 1, circular::meandeviation)) * 24 / 2 / pi

  peaks_df <- peaks_df %>%
    mutate(mean_time_pred = round(mean_local_times, 2),
           sd_time_pred = round(sd_local_times, 2),
           min_time_pred = round(apply(local_times, 1, min), 2),
           max_time_pred = round(apply(local_times, 1, max), 2),
           weighted_mean_time_pred = round(weighted_times_vec, 2)) %>%
    mutate(time_pred_diff = time_1st_peak - weighted_mean_time_pred)

  # --- Attach metadata and compute prediction error ---
  if (mode == 'test') {
    # Fallback: if no peaks found above threshold, use weighted mean of local predictions
    peaks_df$time_1st_peak <- ifelse(is.na(peaks_df$time_1st_peak),
                                     peaks_df$weighted_mean_time_pred,
                                     peaks_df$time_1st_peak)
    peaks_df$Actual_Time <- as.numeric(object[['Metadata']][['Test']][['Time']]) %% 24
    peaks_df$Pred_Error <- circular_difference(peaks_df$time_1st_peak, peaks_df$Actual_Time)
    peaks_df <- peaks_df %>%
      mutate(Group_1 = object[['Metadata']][['Test']][['Group_1']],
             Group_2 = object[['Metadata']][['Test']][['Group_2']],
             Group_3 = object[['Metadata']][['Test']][['Group_3']],
             Replicate = object[['Metadata']][['Test']][['Replicate']])
  } else {
    peaks_df$Actual_Time <- object$Train$Normalised_Data$Time_mod24
    peaks_df$Pred_Error <- circular_difference(peaks_df$time_1st_peak, peaks_df$Actual_Time)
    peaks_df$Group <- object$Train$Normalised_Data$Group
    peaks_df$Group_2 <- object$Train$Normalised_Data$Group_2
    peaks_df$Group_3 <- object$Train$Normalised_Data$Group_3
    peaks_df$Replicate <- object$Train$Normalised_Data$Replicate
  }

  object[[data_slot]][['Results_df']] <- peaks_df
  return(object)
}

# Backward-compatible wrappers
second_peaks_fun_train <- function(object, minpeakheight = -Inf, minpeakdistance = 1,
                                   nups = 1, ndowns = 0, threshold = 0, npeaks = 2) {
  second_peaks_fun(object, mode = 'train', minpeakheight = minpeakheight,
                   minpeakdistance = minpeakdistance, nups = nups, ndowns = ndowns,
                   threshold = threshold, npeaks = npeaks)
}
second_peaks_fun_test <- function(object, minpeakheight = -Inf, minpeakdistance = 1,
                                  nups = 1, ndowns = 0, threshold = 0, npeaks = 2) {
  second_peaks_fun(object, mode = 'test', minpeakheight = minpeakheight,
                   minpeakdistance = minpeakdistance, nups = nups, ndowns = ndowns,
                   threshold = threshold, npeaks = npeaks)
}
