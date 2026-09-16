#' @title Create m/z and retention time ranges for target compounds
#' @description
#' Creates ranges around given m/z and retention time based on given data and
#' allowed deviance.
#' @param msData [MsExperiment()] object
#' @param dbData target database, output of [createTargetList()]
#' @param ppm allowed deviance in ppm around given m/z value
#' @param rtdev allowed deviance in seconds of retention time. Defines the
#'     search window in the time dimension.
#' @return A list containing the m/z and retention time ranges for all given
#'     target compounds
#' @export
#' @author Pablo Vangeenderhuysen
createRanges <- function(msData, dbData, ppm, rtdev) {
  
  spectra <- msData@spectra
  
  # Extract once (CRITICAL)
  basePeakMZ <- spectra$basePeakMZ
  rt_vals <- rtime(spectra)
  
  mz_min <- min(basePeakMZ, na.rm = TRUE)
  mz_max <- max(basePeakMZ, na.rm = TRUE)
  
  rt_min <- min(rt_vals, na.rm = TRUE)
  rt_max <- max(rt_vals, na.rm = TRUE)
  
  # ---- MZ RANGES ----
  mzmed <- as.numeric(dbData$`m/z`)
  delta_mz <- mzmed * ppm * 1e-6  # vectorized
  
  mz_low  <- mzmed - delta_mz
  mz_high <- mzmed + delta_mz
  
  # clamp to global bounds (vectorized)
  mz_low  <- pmax(mz_low,  mz_min)
  mz_high <- pmin(mz_high, mz_max)
  
  # handle inverted ranges (rare edge case)
  invalid <- mz_high < mz_low
  if (any(invalid)) {
    mz_low[invalid]  <- mz_min
    mz_high[invalid] <- mz_min
  }
  
  mzRanges <- cbind(mz_low, mz_high)
  
  # ---- RT RANGES ----
  rtmed <- as.numeric(dbData$tr)
  half_window <- rtdev / 2
  
  rt_low  <- rtmed - half_window
  rt_high <- rtmed + half_window
  
  # special edge handling (vectorized)
  left_out  <- (rt_high - 2) < rt_min
  right_out <- (rt_low + 2) > rt_max
  
  if (any(left_out)) {
    rt_low[left_out]  <- rt_min
    rt_high[left_out] <- rt_min + 10
  }
  
  if (any(right_out)) {
    rt_low[right_out]  <- rt_max - 10
    rt_high[right_out] <- rt_max
  }
  
  # clamp overlaps
  rt_low  <- pmax(rt_low,  rt_min)
  rt_high <- pmin(rt_high, rt_max)
  
  rtRanges <- cbind(rt_low, rt_high)
  
  return(list(mzRanges, rtRanges))
}