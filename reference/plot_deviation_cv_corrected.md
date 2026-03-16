# Cross-validation results on the training data

Visualisation of results following `train_cv` function. Plots prediction
error against actual time. This function can be useful when systematic
deviations (eg chronotype) are suspected, in which case error will be
adjusted by the median. Care should be taken otherwise, as this could
also mask the algorithm's systematic error, therefore
`plot_deviation_cv_original` and `plot_cv_res` are the recommended first
steps

## Usage

``` r
plot_deviation_cv_corrected(list_cv)
```

## Arguments

- list_cv:

  list containing the output of `train_cv` function

## Value

Returns `TimeTeller` prediction error plotted against the actual sample
time

## Author

Vadim Vasilyev
