# Rhythmicity analysis of training data

Calculates rhythmicity results for the training data

## Usage

``` r
choose_genes_tt(
  object,
  group1,
  group2,
  group3,
  replicate,
  method = "population",
  parallel = TRUE,
  cores = 4
)
```

## Arguments

- object:

  list containing TimeTeller training and test models following
  `train_model` and `test_model` respectively

- group1:

  if only a subset of Group_1 (Metadata provided) should be used

- group2:

  if only a subset of Group_2 (Metadata provided) should be used

- group3:

  if only a subset of Group_3 (Metadata provided) should be used

- replicate:

  if only a subset of Replicate (Metadata provided) should be used

- method:

  method used for rhythmicity analysis. Default is `'population'` as in
  <https://tbiomed.biomedcentral.com/articles/10.1186/1742-4682-11-16>

- parallel:

  if TRUE, parallel computation using `foreach` and `doParallel`
  packages will be used. Default is TRUE

- cores:

  if using parallel computation, how many cores should be used. Default
  is 4

## Value

Returns plots of theta calculation curves produced using standard
[`graphics::plot`](https://rdrr.io/r/graphics/plot.default.html)

## References

Cornelissen, G., 2014. Cosinor-based rhythmometry. Theoretical Biology
and Medical Modelling, 11(1), pp.1-24.

## Author

Vadim Vasilyev
