#' Sum the intensities within one spectrum for EIC extraction
#'
#' @noRd
.sumIntensities <- function(x, ...) {
    if (nrow(x)) {
        cbind(
            mz = NA_real_,
            intensity = sum(x[, "intensity"], na.rm = TRUE)
        )
    } else {
        cbind(mz = NA_real_, intensity = NA_real_)
    }
}

#' Integrate single chromatographic peak in single file
#'
#' @param file_path `character(1)` path to .mzML or .mzXML file
#' @param rt `numeric(1)` search retention time for target in minutes
#' @param mz `numeric(1)` search m/z of target compound
#' @param ppm `numeric(1)` allowed ppm deviance
#' @param rtdev `numeric(1)` RT window in seconds
#' @param smoothing `logical(1)` `TRUE` or `FALSE`
#'
#' @import MsExperiment
#' @importFrom Spectra MsBackendMzR
#' @importFrom Spectra filterMzRange
#' @importFrom Spectra filterDataOrigin
#' @importFrom Spectra filterRt
#' @importFrom Spectra dataOrigin
#' @importFrom Spectra addProcessing
#' @importFrom signal sgolayfilt
#' @importFrom BiocParallel SnowParam
#' @importFrom xcms intensity
#' @importFrom pracma trapz
#'
#' @returns results of peak integration
#'
#' @export

integrateSinglePeak <- function(file_path, rt, mz, ppm, rtdev, smoothing) {
    data <- MsExperiment::readMsExperiment(
        spectraFiles = file_path,
        backend = MsBackendMzR(),
        BPPARAM = SnowParam(workers = 1)
    )
    dbData <- data.frame(rt * 60, mz)
    colnames(dbData) <- c("tr", "m/z")
    range <- createRanges(data, dbData, ppm, rtdev)
    mzRange <- range[[1]]
    rtRange <- range[[2]]
    spectra <- data@spectra
    spectra <-
        filterRt(spectra, rtRange)
    spectra <-
        filterMzRange(spectra, mzRange)
    sfs_agg <-
        addProcessing(spectra, .sumIntensities)
    eic <-
        cbind(
            rtime(sfs_agg),
            unlist(intensity(sfs_agg), use.names = FALSE)
        )
    rt <- eic[, 1]
    int <- eic[, 2]
    int[which(is.na(int))] <- 0
    smoothed <- sgolayfilt(int, p = 3, n = 7)
    if (smoothing == TRUE) {
        int <- smoothed
        int[int < 0] <- 0
    }
    border <- find_peak_points(rt, smoothed, dbData$tr)
    x <- rt[border["left"]:border["right"]]
    y <- int[border["left"]:border["right"]]
    rt_list <- list(rt)
    int_list <- list(int)
    x_list <- list(x)
    y_list <- list(y)
    auc <- trapz(x, y)
    pop <- length(x)
    qscore <- qscoreCalculator(x, y)
    found_rt <- rt[border["peak_index"]]
    max_int <- int[border["peak_index"]]
    plot(
        NULL,
        xlim = range(unlist(rt_list) / 60, unlist(x_list) / 60),
        ylim = c(0, max(unlist(int_list), unlist(y_list))),
        type = "n",
        main = "",
        sub = "",
        xlab = "rt",
        ylab = "int"
    )
    mapply(
        function(rt,
    int,
    x,
    y,
    rt_int_color)
        # x_y_color)
        {
            lines(rt / 60, int, col = "red") # Line plot for rt and int
            points(rt / 60, int, col = "red") # Points for rt and int
            a <- x[1] # left integration border
            b <- tail(x, 1) # right integration border
            index_a <- which.min(abs(rt - a))
            index_b <- which.min(abs(rt - b))
            polygon(
                c(rt[index_a] / 60, rt[index_a:index_b] / 60, rt[index_b] / 60),
                c(0, int[index_a:index_b], 0),
                col = adjustcolor(rt_int_color, alpha.f = 0.3),
                border = NA
            )
        },
        rt_list,
        int_list,
        x_list,
        y_list,
        rt_int_color = "red"
    )
    results <- data.frame(auc, pop, qscore, max_int, found_rt)
    results
}
