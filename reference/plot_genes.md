# Visualise multiple gene expression

Quick exploratory plot of temporal gene expression across groups.
Recommended use when there are no replicates for clarity of
visualisation

## Usage

``` r
plot_genes(object, genes, group1, group2)
```

## Arguments

- object:

  list containing TimeTeller training model following `train_model`

- genes:

  genes of interest for visualisation

- group1:

  groups of interest for plotting. If not supplied, all levels from
  `Group_1` will be plotted

- group2:

  groups of interest for plotting. If not supplied, all levels from
  `Group_2` will be plotted

## Value

Returned is the object of class `ggplot`

## Author

Vadim Vasilyev
