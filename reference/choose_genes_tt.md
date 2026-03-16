# Rhythmicity analysis of training data

Calculates population cosinor rhythmicity results for the training data,
including phase consistency (MRL) and amplitude reliability (CV of
rAMP). Genes are ranked by a weighted composite of p-value, R-squared,
MRL, and amplitude CV.

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
  cores = 6,
  rank_weights = c(pval = 1, rsquared = 1, mrl = 1, amp_cv = 1)
)
```

## Arguments

- object:

  list containing TimeTeller training model following `train_model`

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
  is 6

- rank_weights:

  named numeric vector of weights for composite ranking. Names must be
  `pval`, `rsquared`, `mrl`, `amp_cv`. Default is equal weighting.

## Value

Returns the updated object with `Rhythmicity_Results` data frame

## References

Cornelissen, G., 2014. Cosinor-based rhythmometry. Theoretical Biology
and Medical Modelling, 11(1), pp.1-24.

## Author

Vadim Vasilyev
