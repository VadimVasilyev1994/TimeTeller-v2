# Train the TimeTeller model

Trains the model using provided expression matrix and metadata.

## Usage

``` r
train_model(
  exp_matrix,
  genes,
  group_1,
  group_2,
  group_3,
  time,
  replicate,
  mat_normalised = TRUE,
  treat_independently = TRUE,
  combine_for_norm = FALSE,
  parallel_comp = FALSE,
  cores = 4,
  method = "intergene",
  grouping_vars = c("Group"),
  num_PC = 3,
  log_thresh,
  epsilon = 0.4,
  eta = 0.35,
  cov_estimate = "normal",
  alpha_par = 0.75,
  num_interp_points = 144,
  interp_method = "perpchip",
  cov_path = "spline",
  minpeakheight = -Inf,
  minpeakdistance = 1,
  nups = 1,
  ndowns = 0,
  threshold = 0,
  npeaks = 2,
  diagnose_pd = FALSE
)
```

## Arguments

- exp_matrix:

  Matrix or data frame with features in rows and samples in columns

- genes:

  genes used to train the TimeTeller model

- group_1:

  vector containing metadata for each sample (eg individual, organ)

- group_2:

  vector containing metadata for each sample (eg individual, organ).
  Leave empty if unavailable

- group_3:

  vector containing metadata for each sample (eg individual, organ).
  Leave empty if unavailable

- time:

  vector containing timing of each sample. Please appreciate the
  difference between CT and ZT

- replicate:

  vector containing replicate information if available

- mat_normalised:

  was the expression matrix supplied already normalised (eg RMA for
  microarray, CPM/TPM for RNA, etc) or not. If not,
  [`edgeR::cpm`](https://rdrr.io/pkg/edgeR/man/cpm.html) function will
  be used (pseudocount of 2) to obtain log normalised values for
  downstream analysis. Default is TRUE

- treat_independently:

  should replicates be averaged (might be useful for technical
  replicates, not recommended for biological replicates). Default is
  TRUE

- combine_for_norm:

  should replicates be combined for normalisation. If TRUE, replicates
  will be combined for normalisation purposes (eg in `timecourse`
  normalisation, replicates will be combined and treated as one
  timeseries). Default is FALSE

- parallel_comp:

  if TRUE, parallel computation using `foreach` and `doParallel`
  packages will be used. Default is FALSE

- cores:

  if using parallel computation, how many cores should be used. Default
  is 4

- method:

  method used for normalisation. Must be one of the `'intergene'`,
  `'timecourse'`, `'timecourse_matched'`, `'clr'`

- grouping_vars:

  group variables (eg `'group_1'`, `'group_2'` or their combination
  `c('group_1','group_2')`) used for `timecourse_matched` normalisation

- num_PC:

  how many principal components should be used for the analysis. Default
  is 3

- log_thresh:

  log threshold selected. This is an important parameter and should be
  chosen carefully. Please read
  [`help(choose_logthresh_plot)`](https://vadimvasilyev1994.github.io/TimeTeller/reference/choose_logthresh_plot.md)
  for further information

- epsilon:

  hyperparameter value used in theta calculation. No reason to change
  this, more information at
  <https://www.biorxiv.org/content/10.1101/2023.03.14.532177v1>. Default
  is 0.40

- eta:

  hyperparameter value used in theta calculation. No reason to change
  this, more information at
  <https://www.biorxiv.org/content/10.1101/2023.03.14.532177v1>. Default
  is 0.35

- cov_estimate:

  type of covariance estimate used. If `'normal'` is selected,
  [`Rfast::mvnorm.mle`](https://rdrr.io/pkg/Rfast/man/mvnorm.mle.html)
  will be used. If `'robust'` is selected,
  [`rrcov::CovMcd`](https://rdrr.io/pkg/rrcov/man/CovMcd.html) will be
  used instead. This can be useful when there is a lot of training data,
  some of which is suspected to have high technical noise. Since we
  don't want to ignore biological noise (eg chronotype) and it's often
  hard to tell the two apart, default is `'normal'`

- alpha_par:

  this is following `rrcov` package (in case `'robust'` was chosen
  above) and is the parameter controlling the size of the subsets over
  which the determinant will be minimized. Allowed values are between
  0.5 and 1

- num_interp_points:

  number of interpolation points to use for the mean spline and
  covariance matrices. Default is 144 and this corresponds to 24 \* 60 /
  144 = 10 minutes

- interp_method:

  method used for spline interpolation. Must be one of the `'perpchip'`
  (which uses an adapted periodic version of `pchip` from `pracma`
  package) or `'standard'` (which uses periodic version of standard
  `spline` function from `splines` package).

- cov_path:

  method for covariance matrix interpolation. Must be one of the
  `'spline'` (in which case covariance matrix will be interpolated
  element-wise using `interp_method` selection) or `'fisherrao'` (in
  which case the geodesics based on the Fisher–Rao metric for Gaussian
  distributions will be used).

- minpeakheight:

  this is used internally with `pracma` package to locate peaks. No
  reason to change, however included for future testing / development.
  Please consult
  [`?pracma::findpeaks`](https://rdrr.io/pkg/pracma/man/findpeaks.html),
  if needed, for further information

- minpeakdistance:

  this is used internally with `pracma` package to locate peaks. No
  reason to change, however included for future testing / development

- nups:

  this is used internally with `pracma` package to locate peaks. No
  reason to change, however included for future testing / development

- ndowns:

  this is used internally with `pracma` package to locate peaks. No
  reason to change, however included for future testing / development

- threshold:

  this is used internally with `pracma` package to locate peaks. No
  reason to change, however included for future testing / development

- npeaks:

  this is used internally with `pracma` package to locate peaks. No
  reason to change, however included for future testing / development

- diagnose_pd:

  if TRUE, prints the proportion of interpolated covariance matrices
  that required nearPD correction. Useful for assessing interpolation
  quality. Default is FALSE

## Value

Returned is the rich object of class `list` containing the TimeTeller
model for further analysis

## References

Vlachou, D., Veretennikova, M., Usselmann, L., Vasilyev, V., Ott, S.,
Bjarnason, G.A., Dallmann, R., Levi, F. and Rand, D.A., 2023.
TimeTeller: a tool to probe the circadian clock as a multigene dynamical
system. bioRxiv, pp.2023-03.

Amari, S.I. and Nagaoka, H., 2000. Methods of information geometry (Vol.
191). American Mathematical Soc..

## Author

Vadim Vasilyev

## Examples

``` r
library("TimeTeller")
tt_model <- train_model(exp_matrix = bjarn_data$expr_mat, genes = bjarn_data$probes_used,
                        group_1 = bjarn_data$group, time = bjarn_data$time, log_thresh = -5)
#> Please check the information provided:
#> # Number of genes used: 9  ( 202861_at 205251_at 209750_at 209782_s_at 209824_s_at 210971_s_at 213462_at 221045_s_at 39549_at )
#> # Number of Group_1 levels: 10  ( Ind_1 Ind_2 Ind_3 Ind_4 Ind_5 Ind_6 Ind_7 Ind_8 Ind_9 Ind_10 )
#> # Number of Group_2 levels: 1  ( NA )
#> # Number of Group_3 levels: 1  ( NA )
#> # Number of Replicate levels: 1  ( NA )
#> # Number of unique time points Used for training: 6  ( 8 12 16 20 24 28 )
#> Replicates (if there are any) will be treated as independent observations for the training model
#> Using intergene normalisation
#> Data loaded and successfully normalised
#> Calculating projections and interpolating densities...
#> Using standard Covariance and Location estimation
#> Calculating likelihoods...
#> Finished 1 of 6 
#> Finished 2 of 6 
#> Finished 3 of 6 
#> Finished 4 of 6 
#> Finished 5 of 6 
#> Finished 6 of 6 
#> Calculating additional info...
#> Finished!
```
