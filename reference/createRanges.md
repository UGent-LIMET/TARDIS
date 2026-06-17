# Create m/z and retention time ranges for target compounds

Creates ranges around given m/z and retention time based on given data
and allowed deviance.

## Usage

``` r
createRanges(msData, dbData, ppm, rtdev)
```

## Arguments

- msData:

  `MsExperiment()` object

- dbData:

  target database, output of
  [`createTargetList()`](https://pablovgd.github.io/TARDIS/reference/createTargetList.md)

- ppm:

  allowed deviance in ppm around given m/z value

- rtdev:

  allowed deviance in seconds of retention time. Defines the search
  window in the time dimension.

## Value

A list containing the m/z and retention time ranges for all given target
compounds

## Author

Pablo Vangeenderhuysen
