# Changelog

## TARDIS 1.1

### Changes in 1.1.2

- Fix bug with change commas in input files to points, introduced in
  1.1.1

### Changes in 1.1.1

- Change commas in input files to points
- NA handling logic changed: impute NAs inside a peak, convert to zeros
  for NAs outside of peak

### Changes in 1.1.0

- Implement polarity filter

## TARDIS 1.0

### Changes in 1.0.1

- Quick fix for integrateSinglePeak

### Changes in 1.0.0

- Fix issue \#32
- Implement output of input parameters, see issue \#33
- Update documentation
- Small changes in GUI

## TARDIS 0.1

### Changes in 0.1.7

- Fix issue \#28

### Changes in 0.1.6

- Add case study vignette

### Changes in 0.1.5

- Fix regarding issue \#24

### Changes in 0.1.4

- Hotfix custom mass range.

### Changes in 0.1.3

- Added functionality to tardisPeaks to allow MsExperiment object as
  input.

### Changes in 0.1.2

- Various small fixes and typo corrections.
- Added quick start vignette.
- Fixed bug when intensities are all zero and/or constant.

### Changes in 0.1.1

- Refactor
  [`find_peak_points()`](https://pablovgd.github.io/TARDIS/reference/find_peak_points.md)
  and add unit tests and documentation.
- General improvement of code readability and documentation.
- Get correct rt of ISTD for RT alignment
- New function
  [`checkScans()`](https://pablovgd.github.io/TARDIS/reference/checkScans.md)
  to check faulty input files that miss scans.
- Setting intensity filter to zero disables to filter to retain `NA`.
