# Theta calculation of individual samples

This displays `theta` and flat likelihood percentage calculation for
individual samples. Can be useful for diagnostics and troubleshooting
the model

## Usage

``` r
plot_ind_curve(object, sample_num, logthresh, train_or_test = "test")
```

## Arguments

- object:

  list containing TimeTeller training and test models following
  `train_model` and `test_model` respectively

- sample_num:

  number of sample to be displayed

- logthresh:

  which log threshold should be used for visualisation

- train_or_test:

  is the sample coming from `train` or `test` data

## Value

Returns plots of theta calculation curves produced using standard
[`graphics::plot`](https://rdrr.io/r/graphics/plot.default.html)

## Author

Vadim Vasilyev
