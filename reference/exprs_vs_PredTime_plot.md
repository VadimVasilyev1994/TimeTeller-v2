# Visualise the estimated temporal expression of test data

Plots both known temporal expression of the training data and estimated
temporal expression of the test data. Times used for the test data are
the ones predicted by TimeTeller. Supports intergene normalisation,
TimeTeller-predicted times for training data, and custom dataset labels
and gene name mappings.

## Usage

``` r
exprs_vs_PredTime_plot(
  object,
  genes,
  new_names,
  theta_thresh,
  xlim_l = 10,
  xlim_u = 20,
  dataset_manual = FALSE,
  dataset_info = "Test",
  tt_time = FALSE,
  intergene = FALSE
)
```

## Arguments

- object:

  list containing TimeTeller training and test models following
  `train_model` and `test_model` respectively

- genes:

  genes / features of interest used for plotting

- new_names:

  optional named vector mapping gene IDs to display names for facet
  labels

- theta_thresh:

  threshold used for theta classification into 'Good' and 'Bad' clocks

- xlim_l:

  lower limit for x axis. Default is 10

- xlim_u:

  upper limit for x axis. Default is 20

- dataset_manual:

  if TRUE, uses `dataset_info` as the test dataset label instead of
  'Test'

- dataset_info:

  character string to label the test dataset when
  `dataset_manual = TRUE`

- tt_time:

  if TRUE, uses TimeTeller-predicted times instead of actual times for
  training data

- intergene:

  if TRUE, applies per-sample intergene normalisation before plotting

## Value

Returned is the `ggplot` object

## Author

Vadim Vasilyev
