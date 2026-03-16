# Covariance matrix ellipsoid volume

This displays the volume (area) of the hyper-ellipse at a given
confidence level (<https://datavis.ca/papers/ellipses.pdf>). Can be
useful for diagnostics and further analysis of the covariance matrix
interpolation

## Usage

``` r
check_cov_interp(object, alpha = 0.95)
```

## Arguments

- object:

  list containing TimeTeller training and test models following
  `train_model`

- alpha:

  confidence level used for volume calculation

## Value

Returns the plot of volume of hyper-ellipse for all the local
projections using
[`graphics::matplot`](https://rdrr.io/r/graphics/matplot.html)

## Author

Vadim Vasilyev
