#' filter Spectra to single peak in single sample
#'
#' @param spectra a `Spectra` object
#' @param dataOrigin file that contains the spectra
#' @param rt_range `numeric(2)` rt range of compound
#' @param mz_range `numeric(2)` mz range of compound
#' @returns returns filtered `Spectra` object
filterSingle <- function(spectra, dataOrigin, rt_range, mz_range) {
    spectra <- spectra |>
        filterDataOrigin(dataOrigin) |>
        filterRt(rt_range) |>
        filterMzRange(mz_range)
    spectra
}

#' Sum Intensities of Spectra
#' @noRd
.sum_intensities <- function(x, ...) {
    if (nrow(x)) {
        cbind(
            mz = NA_real_,
            intensity = sum(x[, "intensity"], na.rm = TRUE)
        )
    } else {
        cbind(mz = NA_real_, intensity = NA_real_)
    }
}

#' @title Function to extract EIC from Spectra object
#'
#' @param spectra a `Spectra` object.
#'
#' @return two column matrix with rt and int
#'
#' @author Pablo Vangeenderhuysen
#'
#'
#' @export
extract_eic <- function(spectra) {
    sfs_agg <-
        addProcessing(spectra, .sum_intensities)
    eic <-
        cbind(
            rtime(sfs_agg),
            unlist(intensity(sfs_agg), use.names = FALSE)
        )
    rownames(eic) <- NULL
    colnames(eic) <- c("rt", "int")
    eic
}

# vectorize the loop to make it faster
filterSingle_extractEIC <- function(spectra, dataOrigin, rt_range, mz_range) {
  spectra_filtered <- spectra |>
    filterDataOrigin(dataOrigin) |>
    filterRt(rt_range) |>
    filterMzRange(mz_range)
  
  ints <- intensity(spectra_filtered)
  rts  <- rtime(spectra_filtered)
  
  # if there is NA, impute via Linear Interpolation
  if (anyNA(ints)) {
    ints <- imputeLinInterpol(ints)
  }
  # if there is NA, impute via Linear Interpolation
  if (anyNA(rts)) {
    rts <- imputeLinInterpol(rts)
  }
  
  out <- cbind(
    rt  = rts,
    # na.rm = TRUE ignores NA (treats it as 0)
    int = vapply(ints, function(x) if (length(x)) sum(x, na.rm = TRUE) else NA_real_, numeric(1))
  )
  out
}