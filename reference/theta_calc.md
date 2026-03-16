# Compute Theta (clock dysfunction metric) for each sample

Theta measures the proportion of the likelihood curve that exceeds a
cosine envelope, capturing how peaked vs diffuse the likelihood is. A
well-functioning clock produces a sharply peaked likelihood (low Theta);
a disrupted clock produces a broad/flat likelihood (high Theta).

## Usage

``` r
theta_calc(object, mode = "train", epsilon = NULL, eta = NULL)
```

## Arguments

- object:

  TimeTeller list object

- mode:

  either 'train' or 'test'

- epsilon:

  envelope amplitude parameter (default 0.4)

- eta:

  envelope scaling parameter (default 0.35)

## Value

Updated object with Theta values
