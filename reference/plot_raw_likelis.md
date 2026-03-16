# Likelihoods of individual samples

This displays likelihoods (for all the local projections) for individual
samples. Can be useful for diagnostics and troubleshooting the model

## Usage

``` r
plot_raw_likelis(object, sample_num, logthresh, train_or_test = "test")
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

Returns a plot of likelihood curves produced using
[`graphics::matplot`](https://rdrr.io/r/graphics/matplot.html)

## Author

Vadim Vasilyev
