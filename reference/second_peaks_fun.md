# Extract peak information and build results data frame

Finds peaks in the averaged likelihood curve, computes prediction
statistics (weighted mean time, circular mean/SD across local
projections), and assembles the full results data frame with metadata.

## Usage

``` r
second_peaks_fun(
  object,
  mode = "train",
  minpeakheight = -Inf,
  minpeakdistance = 1,
  nups = 1,
  ndowns = 0,
  threshold = 0,
  npeaks = 2
)
```

## Arguments

- object:

  TimeTeller list object

- mode:

  either 'train' or 'test'

- minpeakheight, minpeakdistance, nups, ndowns, threshold, npeaks:

  passed to pracma::findpeaks

## Value

Updated object with Results_df in the appropriate data slot
