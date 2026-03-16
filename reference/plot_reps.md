# Visualise gene expression across groups and replicates

Quick exploratory plot of temporal gene expression across groups and
replicates for the training data.

## Usage

``` r
plot_reps(object, gene, group1, group2)
```

## Arguments

- object:

  list containing TimeTeller training model following `train_model`

- gene:

  gene of interest for visualisation

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
