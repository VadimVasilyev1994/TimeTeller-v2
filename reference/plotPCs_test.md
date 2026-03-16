# Assess precision of time predictions

Useful when assessing the structure of the test data and the precision
of `TimeTeller` phase estimates

## Usage

``` r
plotPCs_test(object, ymin = 0, ymax = 24)
```

## Arguments

- object:

  list containing TimeTeller training and test models following
  `train_model` and `test_model` respectively

- ymin:

  lower limit for y axis. Helps with visualisation when samples cluster
  in a particular part of the day (eg biopsy samples taken mostly during
  day time)

- ymax:

  upper limit for y axis

## Value

Returns test data projection onto the first 4 PCs vs estimated time

## Author

Vadim Vasilyev
