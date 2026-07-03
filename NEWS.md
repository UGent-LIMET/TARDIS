# TARDIS 2.0

## Changes in 2.8.0
- changed Spectra backend from MsBackendMzR to MsBackendOfflineSql, since the v2.7.1 code didn't run successfully in lab computers (used too much memory)

## Changes in 2.7.1
- fixed bug: account for different column names in createTargetList input

## Changes in 2.7.0
- generate std plots in screening phase
- works also when user inputs a qc_pattern that is not "QC"
- added smoothing parameter "smoothing_order" for the smoothing function

## Changes in 2.6.1
- Fixed bugs (smoothing function null handling)

## Changes in 2.6.0
- Fixed plotting issues (null values handled wrongly by previous imputation step)
- Added pval_cutoff parameter to GUI and functions

## Changes in 2.5.0
- Fixed null handling by linear imputation of missing values
- No null plots!

## Changes in 2.4.0.2
- Fixed logic in rtAlignment() function (retention time alignment)

## Changes in 2.4.0b
- Fixed issue: detecting overlapping peaks
- Code detects unimodal/multimodal distribution and detects peak accordingly

## Changes in 2.3.0b
- Fixed issue: retention time alignment
- Retention time is aligned, then shifted to the correct retention time.

## Changes in 2.2.0b
- Fixed issue: detecting only the highest peak
- Ignore local maxima that is lower than 10% of the maximum peak intensity, not 50% (previous)

## Changes in 2.1.0b
- Improved algorithm efficiency and parallel computing code
- Achieved lower runtime in comparison to the original code by Pablo

## Changes in 2.0.0b
- Fixed 'missing value where TRUE/FALSE needed' crash in `tardisPeaks`.
- Added robust handling for empty GUI inputs.

## Changes in 2.0.0
- Edited peakdet.R and peaks_with_tardis.R
- Added functions for repetitive code
- Fixed issue: detecting only one peak in two overlapping peaks
- Implemented parallel processing in tardisPeaks function (the number of cores need to be specified in the user interface)
- Fixed issue (still needs improvement): some peaks get cut off the sides of the plot window
- Implemented baseline correction to smoothing method
- Created newSmoothing.R: a new smoothing function (not implemented in main function)

# TARDIS 1.0

## Changes in 1.1.0
- Implement polarity filter

## Changes in 1.0.1
- Quick fix for integrateSinglePeak

## Changes in 1.0.0
- Fix issue #32
- Implement output of input parameters, see issue #33
- Update documentation
- Small changes in GUI

# TARDIS 0.1

## Changes in 0.1.7
- Fix issue #28

## Changes in 0.1.6
- Add case study vignette

## Changes in 0.1.5
- Fix regarding issue #24

## Changes in 0.1.4
- Hotfix custom mass range.

## Changes in 0.1.3
- Added functionality to tardisPeaks to allow MsExperiment object as input.

## Changes in 0.1.2
- Various small fixes and typo corrections.
- Added quick start vignette.
- Fixed bug when intensities are all zero and/or constant.

## Changes in 0.1.1
- Refactor `find_peak_points()` and add unit tests and documentation.
- General improvement of code readability and documentation.
- Get correct rt of ISTD for RT alignment
- New function `checkScans()` to check faulty input files that miss scans.
- Setting intensity filter to zero disables to filter to retain `NA`.
