# Compute flat region contribution to Theta

Quantifies what proportion of Theta is attributable to flat (saturated)
regions in the likelihood curve — i.e. regions clamped at the log
threshold. High flat contribution suggests Theta is inflated by the
threshold rather than genuine clock dysfunction.

## Usage

``` r
calc_flat_theta_contrib(object, mode = "train")
```

## Arguments

- object:

  TimeTeller list object

- mode:

  either 'train' or 'test'

## Value

Updated object with flat contribution data frame
