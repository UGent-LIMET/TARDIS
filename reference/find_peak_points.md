# Simple peak detection algorithm on chromatographic data

Basic peak picking algorithm that finds a targeted chromatographic peak
given based on retention time and intensities and expected retention
time of the target signal. The chromatographic signal is expected to
contain only signal from a single ion. See details below for information
on the algorithm.

## Usage

``` r
find_peak_points(
  rtime = numeric(),
  intensity = numeric(),
  targetRtime = numeric(),
  .check = TRUE
)
```

## Arguments

- rtime:

  `numeric` with retention times.

- intensity:

  `numeric` (same length as `rtime`) with the signal intensities.

- targetRtime:

  `numeric(1)` with the expected retention time of the target peak.

- .check:

  `logical(1)` whether input parameters should be checked. Use
  `.check = FALSE` only if validity of the input parameters was checked
  in an upstream function.

## Details

The method:

- identifies first all loacal maxima in `intensity`, then

- removes local maxima with an intensity lower than 50% of the maximum
  intensity in `intensity` and finally

- reports the above defined tentative peak apex with an retenton time
  closest to the provided `targetRtime`.

The closest local minima that are at least 2 data points from the apex
position are reported as the left and right margin of the peak.

## Author

Pablo Vangeenderhuysen
