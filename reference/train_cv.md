# Cross-validate the TimeTeller model

Cross-validation of the TimeTeller model using provided expression
matrix and metadata. Single or multiple groups can be selected for CV

## Usage

``` r
train_cv(
  group_to_leave_out = "group_1",
  genes,
  exp_matrix,
  test_grouping_vars,
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
  npeaks = 2
)
```

## Arguments

- group_to_leave_out:

  This is the group that cross-validation will be done on. For example,
  if `group_1` contains sample individual and `group_2` contains sample
  organ information, selecting `'group_2'` will result in CV done on
  organs as opposed to individuals (ie leave one organ out for testing,
  train on the rest). Must be one of the `'group_1'`, `'group_2'`,
  `'group_3'` or `'replicate'`

- genes:

  genes used to train the TimeTeller model

- exp_matrix:

  matrix or data frame with features in rows and samples in columns

- test_grouping_vars:

  group variables (eg `'group_1'`, `'group_2'` or their combination
  `c('group_1','group_2')`) used for `timecourse_matched` normalisation.
  For cross-validation this will be the same as `grouping_vars`

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

## Value

Returned is the object of class `list` containing the results of
TimeTeller cross-validation analysis

## Author

Vadim Vasilyev
