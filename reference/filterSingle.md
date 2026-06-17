# filter Spectra to single peak in single sample

filter Spectra to single peak in single sample

## Usage

``` r
filterSingle(spectra, dataOrigin, rt_range, mz_range)
```

## Arguments

- spectra:

  a `Spectra` object

- dataOrigin:

  file that contains the spectra

- rt_range:

  `numeric(2)` rt range of compound

- mz_range:

  `numeric(2)` mz range of compound

## Value

returns filtered `Spectra` object
