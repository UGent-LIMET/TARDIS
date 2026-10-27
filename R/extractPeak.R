#' @title Extract multiple EICs
#'
#' @param object An `MsExperiment` object.
#' @param mzRanges A two column matrix with minimum and maximum m/z values.
#' @param rtRanges A two column matrix with minimum and maximum retention
#'     times.
#' @param massRange Optional mass window used to remove scans without peaks in
#'     that window. One additional chromatogram is extracted as a shared scan
#'     mask, avoiding a separate pass over all mzML peak data.
#'
#' @return A list containing the extracted `MChromatograms`, an optional
#'     mass-range scan mask for each sample, and the number of samples.
#' @noRd
.extractChromatograms <- function(object, mzRanges, rtRanges,
                                  massRange = NULL) {
    if (is.null(dim(mzRanges)) && length(mzRanges) == 2) {
        mzRanges <- matrix(mzRanges, nrow = 1)
    } else {
        mzRanges <- as.matrix(mzRanges)
    }
    if (is.null(dim(rtRanges)) && length(rtRanges) == 2) {
        rtRanges <- matrix(rtRanges, nrow = 1)
    } else {
        rtRanges <- as.matrix(rtRanges)
    }
    if (ncol(mzRanges) != 2 || ncol(rtRanges) != 2 ||
            nrow(mzRanges) != nrow(rtRanges)) {
        stop("'mzRanges' and 'rtRanges' must be equally sized two-column matrices.")
    }

    targetCount <- nrow(mzRanges)
    if (!is.null(massRange)) {
        mzRanges <- rbind(mzRanges, massRange)
        rtRanges <- rbind(rtRanges, range(rtRanges))
    }

    chromatograms <- xcms::chromatogram(
        object,
        mz = mzRanges,
        rt = rtRanges,
        return.type = "MChromatograms",
        BPPARAM = BiocParallel::SerialParam()
    )
    sampleCount <- ncol(chromatograms)
    massRangeEics <- NULL
    if (!is.null(massRange)) {
        massRangeEics <- lapply(
            seq_len(sampleCount),
            function(sample) {
                .mChromatogramEic(
                    chromatograms,
                    target = targetCount + 1,
                    sample = sample
                )
            }
        )
    }
    list(
        data = chromatograms,
        massRangeEics = massRangeEics,
        sampleCount = sampleCount
    )
}

#' @title Convert one MChromatogram entry to a TARDIS EIC matrix
#'
#' @param chromatograms An `MChromatograms` object.
#' @param target Integer target index.
#' @param sample Integer sample index.
#'
#' @return A two-column matrix with columns `rt` and `int`.
#' @noRd
.mChromatogramEic <- function(chromatograms, target, sample) {
    currentChromatogram <- methods::slot(chromatograms, ".Data")[[
        target,
        sample
    ]]
    eic <- cbind(
        rt = ProtGenerics::rtime(currentChromatogram),
        int = ProtGenerics::intensity(currentChromatogram)
    )
    rownames(eic) <- NULL
    eic
}

#' @title Retrieve one masked EIC from the MChromatograms result
#'
#' @param chromatograms The result of `.extractChromatograms()`.
#' @param target Integer target index.
#' @param sample Integer sample index.
#'
#' @return A two-column matrix with columns `rt` and `int`.
#' @noRd
.chromatogramEic <- function(chromatograms, target, sample) {
    eic <- .mChromatogramEic(
        chromatograms$data,
        target = target,
        sample = sample
    )
    if (!is.null(chromatograms$massRangeEics)) {
        massRangeEic <- chromatograms$massRangeEics[[sample]]
        retainedRtime <- massRangeEic[, "rt"][
            !is.na(massRangeEic[, "int"])
        ]
        eic <- eic[eic[, "rt"] %in% retainedRtime, , drop = FALSE]
    }
    eic
}
