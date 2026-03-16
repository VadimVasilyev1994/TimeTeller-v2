# Basic_TimeTeller_Analysis

``` r
library(TimeTeller)
#> Warning: replacing previous import 'Matrix::image' by 'graphics::image' when
#> loading 'TimeTeller'
#> Warning: replacing previous import 'graphics::layout' by 'plotly::layout' when
#> loading 'TimeTeller'
#> Warning: replacing previous import 'ggplot2::last_plot' by 'plotly::last_plot'
#> when loading 'TimeTeller'
#> Warning: replacing previous import 'ggplot2::%+%' by 'psych::%+%' when loading
#> 'TimeTeller'
#> Warning: replacing previous import 'ggplot2::alpha' by 'psych::alpha' when
#> loading 'TimeTeller'
#> Warning: replacing previous import 'foreach::when' by 'purrr::when' when
#> loading 'TimeTeller'
#> Warning: replacing previous import 'foreach::accumulate' by 'purrr::accumulate'
#> when loading 'TimeTeller'
#> Warning in rgl.init(initValue, onlyNULL): RGL: unable to open X11 display
#> Warning: 'rgl.init' failed, will use the null device.
#> See '?rgl.useNULL' for ways to avoid this warning.
#> Warning: replacing previous import 'Matrix::cov2cor' by 'stats::cov2cor' when
#> loading 'TimeTeller'
#> Warning: replacing previous import 'dplyr::lag' by 'stats::lag' when loading
#> 'TimeTeller'
#> Warning: replacing previous import 'plotly::filter' by 'stats::filter' when
#> loading 'TimeTeller'
#> Warning: replacing previous import 'Matrix::toeplitz' by 'stats::toeplitz' when
#> loading 'TimeTeller'
#> Warning: replacing previous import 'circular::var' by 'stats::var' when loading
#> 'TimeTeller'
#> Warning: replacing previous import 'circular::sd' by 'stats::sd' when loading
#> 'TimeTeller'
#> Warning: replacing previous import 'Matrix::update' by 'stats::update' when
#> loading 'TimeTeller'
#> Warning: replacing previous import 'Matrix::tail' by 'utils::tail' when loading
#> 'TimeTeller'
#> Warning: replacing previous import 'Matrix::head' by 'utils::head' when loading
#> 'TimeTeller'
#> Warning: replacing previous import 'purrr::transpose' by
#> 'data.table::transpose' when loading 'TimeTeller'
```

## Brief Introduction

This package builds on the work in [Original
publication](https://www.biorxiv.org/content/10.1101/622050v2) and
[Methodology](https://www.biorxiv.org/content/10.1101/2023.03.14.532177v1.full.pdf),
which showed the potential of TimeTeller to act as an independent cancer
biomarker in microarray data, extended the algorithm to RNA-seq data and
performed extensive validation on some of the available mouse, baboon
and human datasets.

TimeTeller is a supervised machine learning tool that analyses the local
circadian clock as a system. It aims to estimate the circadian clock
phase and the level of dysfunction from a single sample by modelling the
multi-dimensional state of the clock.

*TimeTeller* package implements the algorithm including methodology
improvements and tools for visualisation.

## Standard Workflow

In this vignette we will look at two datasets: human oral mucosa and

``` r
str(bjarn_data)
#> List of 4
#>  $ expr_mat   : num [1:1009, 1:60] 7.46 8.32 8.78 7.02 6.43 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:1009] "202861_at" "205251_at" "209750_at" "209782_s_at" ...
#>   .. ..$ : chr [1:60] "data1" "data2" "data3" "data4" ...
#>  $ group      : chr [1:60] "Ind_1" "Ind_1" "Ind_1" "Ind_1" ...
#>  $ time       : num [1:60] 8 12 16 20 24 28 8 12 16 20 ...
#>  $ probes_used: chr [1:9] "202861_at" "205251_at" "209750_at" "209782_s_at" ...
```
