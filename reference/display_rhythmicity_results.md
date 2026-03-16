# Display the results of cosinor rhythmicity analysis

Visualisation of results following `choose_genes_tt` function

## Usage

``` r
display_rhythmicity_results(
  object,
  probes_of_interest,
  info_used = "ranks",
  organism
)
```

## Arguments

- object:

  list containing the output of `train_model` and `choose_genes_tt`
  functions

- probes_of_interest:

  genes / features of interest to highlight

- info_used:

  plot style to output. Must be one of the `'ranks'` or `'original'`.
  Choosing `'ranks'` allows for more interactivity using `ggplotly`

- organism:

  this is added for convenience. Will attempt to use
  [`gprofiler2::gconvert`](https://rdrr.io/pkg/gprofiler2/man/gconvert.html)
  to get Entrez gene symbols for more user-friendly visualisation.
  Therefore, the argument format has to be consistent with that.

## Value

Returns `plot_ly` or `ggplot` object depending on the `info_used`
selection

## Author

Vadim Vasilyev
