# Log Threshold selection

This will iterate though multiple values of logthresh and save the
relevant metrics for each. Can be used to aid in selecting an
appropriate threshold. See
[`?choose_logthresh_plot`](https://vadimvasilyev1994.github.io/TimeTeller/reference/choose_logthresh_plot.md)
for further visualisation

## Usage

``` r
choose_logthresh(
  object,
  max_log = 0,
  min_log = -12,
  by_step = -1,
  train_or_test = "test"
)
```

## Arguments

- object:

  list containing the output of `train_model` and `test_model` functions

- max_log:

  max value of logthresh to be considered. Generally recommended to
  check `[0,-12]` range

- min_log:

  min value of logthresh to be considered

- by_step:

  step size for logthresh optimisation. We have found -0.5 to be
  sufficient, however users may benefit from more granular data
  depending on the application

- train_or_test:

  whether the simulation should be run on `'train'` or `'test'` data

## Value

returns `data.frame` containing the results of log threshold selection

## Author

Vadim Vasilyev
