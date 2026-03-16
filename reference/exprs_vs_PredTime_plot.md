# Visualise the the estimated temporal expression of test data

Plots both known temporal expression of the training data and estimated
temporal expression of the test data. This means that times used for the
test data are the ones predicted by `TimeTeller`

## Usage

``` r
exprs_vs_PredTime_plot(object, genes, theta_thresh, xlim_l = 10, xlim_u = 20)
```

## Arguments

- object:

  list containing TimeTeller training and test models following
  `train_model` and `test_model` respectively

- genes:

  genes / features of interest used for plotting

- theta_thresh:

  threshold used for theta classification into 'Good' and 'Bad' clocks.
  This can be subjective and is used for visualisation and sanity check

- xlim_l:

  lower limit for x axis. Helps with visualisation when samples cluster
  in a particular part of the day (eg biopsy samples taken mostly during
  day time)

- xlim_u:

  upper limit for x axis

## Value

Returned is the `ggplot` object

## Author

Vadim Vasilyev
