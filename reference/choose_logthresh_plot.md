# Log Threshold selection plot

Plots the results of `choose_logthresh` function

## Usage

``` r
choose_logthresh_plot(choose_logthresh_df, cap = 10, perc_flat = 0.02)
```

## Arguments

- choose_logthresh_df:

  `data.frame` containing the output of `choose_logthresh` function

- cap:

  parameter used for clarity of visualisation. Ratio of LogLik for 1st
  and 2nd peaks

- perc_flat:

  option to target a particular proportion of samples with flat
  likelihoods (ie their peak is below the currently selected LogThresh)

## Value

returns `ggplot` object

## Author

Vadim Vasilyev
