# Project test data on the training model

Projects the test data onto the model obtained after running
`train_model` function. Among the outputs are the estimated time for
each sample and the clock dysfunction metric (Theta)

## Usage

``` r
test_model(
  object,
  exp_matrix,
  test_grouping_vars,
  test_group_1,
  test_group_2,
  test_group_3,
  test_replicate,
  test_time,
  mat_normalised_test = TRUE,
  log_thresh,
  parallel_comp = FALSE,
  cores = 4,
  minpeakheight = -Inf,
  minpeakdistance = 1,
  nups = 1,
  ndowns = 0,
  threshold = 0,
  npeaks = 2
)
```

## Arguments

- object:

  object of class list containing timeteller training model (obtained
  after running `train_model` function)

- exp_matrix:

  matrix or data frame containing test data with features in rows and
  samples in columns

- test_grouping_vars:

  groups below that will be used for `timecourse_matched` normalisation.
  See the `vignette` for examples

- test_group_1:

  vector containing metadata for each sample (eg individual, organ).
  These are not required for `intergene` normalisation, however may be
  useful for visualisation. For `timecourse` and `timecourse_matched`
  normalisations, however, these are important and should be supplied in
  the same format as for the training data to ensure consistent results

- test_group_2:

  vector containing metadata for each sample (eg individual, organ)

- test_group_3:

  vector containing metadata for each sample (eg individual, organ)

- test_replicate:

  vector containing replicate information if available

- test_time:

  vector containing replicate information if available. This does not
  affect the results but is useful to calculate the prediction error and
  hence the algorithm accuracy

- mat_normalised_test:

  was the expression matrix supplied already normalised (eg RMA for
  microarray, CPM/TPM for RNA, etc) or not. If not,
  [`edgeR::cpm`](https://rdrr.io/pkg/edgeR/man/cpm.html) function will
  be used (pseudocount of 2) to obtain log normalised values for
  downstream analysis. Default is TRUE

- log_thresh:

  log threshold selected. This is an important parameter and should be
  chosen carefully.Please read
  [`help(choose_logthresh_plot)`](https://vadimvasilyev1994.github.io/TimeTeller/reference/choose_logthresh_plot.md)
  for further information

- parallel_comp:

  if TRUE, parallel computation using `foreach` and `doParallel`
  packages will be used. Default is FALSE

- cores:

  if using parallel computation, how many cores should be used. Default
  is 4

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

Returned is the rich object of class `list` containing the results for
both train and test models.

## Author

Vadim Vasilyev
