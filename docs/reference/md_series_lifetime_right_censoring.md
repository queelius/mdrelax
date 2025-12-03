# Decorate masked data with series system lifetime and right-censoring

Takes component failure times and adds system failure time and
right-censoring indicators.

## Usage

``` r
md_series_lifetime_right_censoring(
  md,
  tau,
  compvar = "t",
  sysvar = "t",
  deltavar = "delta",
  failvar = "k"
)
```

## Arguments

- md:

  Data frame with component failure times (columns t1, t2, ..., tm).

- tau:

  Right-censoring times. Can be scalar or vector of length n.

- compvar:

  Column prefix for component lifetimes. Default is "t".

- sysvar:

  Column name for system lifetime. Default is "t".

- deltavar:

  Column name for censoring indicator. Default is "delta".

- failvar:

  Column name for failed component index. Default is "k".

## Value

Modified data frame with additional columns:

- System lifetime (sysvar)

- Censoring indicator (deltavar): TRUE if censored

- Failed component index (failvar)
