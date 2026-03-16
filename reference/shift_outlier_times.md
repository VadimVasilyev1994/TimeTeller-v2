# Shift predicted phases

Shifts TimeTeller predicted times to more realistic time window IF there
is a second likelihood peak located in that window

## Usage

``` r
shift_outlier_times(object, minT = 8, maxT = 20)
```

## Arguments

- object:

  list containing TimeTeller training and testing models

- minT:

  lower bound for defined time window

- maxT:

  upper bound for defined time window

## Value

Returned is the updated object with corrected timing information

## Author

Vadim Vasilyev
