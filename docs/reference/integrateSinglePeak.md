# Integrate single chromatographic peak in single file

Integrate single chromatographic peak in single file

## Usage

``` r
integrateSinglePeak(file_path, rt, mz, ppm, rtdev, smoothing)
```

## Arguments

- file_path:

  `character(1)` path to .mzML or .mzXML file

- rt:

  `numeric(1)` search retention time for target in minutes

- mz:

  `numeric(1)` search m/z of target compound

- ppm:

  `numeric(1)` allowed ppm deviance

- rtdev:

  `numeric(1)` RT window in seconds

- smoothing:

  `logical(1)` `TRUE` or `FALSE`

## Value

results of peak integration
