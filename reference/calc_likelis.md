# Compute log-likelihood array across all local projections

For each local SVD projection (one per training time point), project the
expression data and evaluate the multivariate normal density at every
interpolation point. The result is a 3D array: (interpolation_points x
samples x local_projections)

## Usage

``` r
calc_likelis(object, mode = "train", diagnose_pd = FALSE)
```

## Arguments

- object:

  TimeTeller list object

- mode:

  either 'train' or 'test'

## Value

Updated object with log-likelihood array (and test projections if
mode='test')
