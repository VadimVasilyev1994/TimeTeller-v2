# Rhythmicity analysis of selected geneset

Convenience function for single cosinor and visualisation of the
selected geneset

## Usage

``` r
geneset_rhythm_info(
  object,
  geneset,
  labels,
  group1,
  group2,
  group3,
  replicate,
  method = "population"
)
```

## Arguments

- object:

  list containing TimeTeller rhythmicity results following
  `choose_genes_tt`

- geneset:

  geneset of interest

- labels:

  genes to highlight on the resulting polar plot

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

## Value

Returns an array containing single cosinor results (MESOR, Amp and
Rhythmicity test Pval) and the polar plot with summary info

## References

Cornelissen, G., 2014. Cosinor-based rhythmometry. Theoretical Biology
and Medical Modelling, 11(1), pp.1-24.

## Author

Vadim Vasilyev
