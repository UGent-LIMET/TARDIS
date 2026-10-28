# Find peak points

Basic peak picking algorithm that finds a targeted chromatographic peak
given the data and expected retention time of the target.

## Usage

``` r
find_peak_points_orig(rtvector, vector, searchrt)
```

## Arguments

- rtvector:

  Vector with retention times.

- vector:

  Vector with intensities.

- searchrt:

  Expected retention time of peak.
