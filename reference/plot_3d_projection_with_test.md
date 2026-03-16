# Visualise the test data

Test data projected onto the training model.

## Usage

``` r
plot_3d_projection_with_test(
  object,
  selected_local_projection,
  density = FALSE,
  opacity = 0.05,
  sig_level = 0.9
)
```

## Arguments

- object:

  list containing TimeTeller training and test models following
  `train_model` and `test_model` respectively

- selected_local_projection:

  training time selected for the local projection

- density:

  should covariance matrix cloud be displayed. Default is FALSE

- opacity:

  opacity level if `density = TRUE`

- sig_level:

  confidence level when `density = TRUE` and is used to control the size
  of the ellipsoid. This uses
  [`rgl::ellipse3d`](https://rdrr.io/pkg/rgl/man/ellipse3d.html)

## Value

Returned is the `plot_ly` object

## Author

Vadim Vasilyev
