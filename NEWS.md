# Stepwise code changes and validation (from v1.1.0)

## Change1
- Added filterSingle_extractEIC function in extractPeak.R and changed relevant parts in peaks_with_tardis.R to use this new function

## Change2
- Made the most basic smoothingSG function in peaks_with_tardis.R, with no changes to algorithm, only replacing repetitive smoothing code
- Added basic NA handling of return values in smoothingSG (not imputing points within the peak) - just to ensure code doesn't crash

## Change3
- Made the most basic checkValidPeak function in peaks_with_tardis.R, with no changes to algorithm, only replacing repetitive smoothing code
- Fixed minor bugs

## Change4
- Added dataHandling function (deals with different QC_pattern other than only "QC", use MsBackendOfflineSql backend instead of MsBackendMzR), change repetitive parts in tardisPeaks() to use dataHandling(), add package import MsBackendSql
- Edited createRanges.R to account for new backend (has some algorithm changes; was arbitrary from the beginning)

## Change5
- Made 5 lapply code blocks, added safe_bind() and standardize_results() (necessary for lapply), change checkvalidpeak input/output type (added null handling)
- Changed max_int_filter null value handling, to fix error "Error in if: missing value where TRUE/FALSE needed"
- edited metrics tables generation code (to preserve clean layout for the output csv tables)

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
