#' @title TARDIS Peak Detection
#'
#' @description
#' Main function of the TARDIS package that is called in the Shiny app.
#' Given data files and a list of targeted compounds it returns the area of
#' those peaks, optional diagnostic plots and several other parameters.
#' See vignette for a detailed tutorial.
#'
#' @param file_path `character(1)` Path to the .mzML or .mzXML files
#'     containing LC-MS data.
#' @param lcmsData  `MsExperiment` MsExperiment containing the data to be
#'     preprocessed. Sampledata should at least include run type and should
#'     match the later provided "QC" or "sample" pattern.
#' @param dbData Output of [createTargetList()]
#' @param ppm `numeric(1)` Allowed deviance from given m/z of targets in ppm.
#' @param rtdev `numeric(1)` Allowed deviance from given retention time of
#'     compound, defines search window for the peak picking algorithm.
#' @param mass_range `numeric(2)` If the user uses data with overlapping mass
#'     windows, only one mass window at the time can be analyzed, specify this
#'     here.
#' @param polarity `character(1)` Ionisation mode to be considered, can be
#'     either "positive" or "negative"
#' @param output_directory `character(1)` Provide directory to store output
#' @param plots_samples `logical(1)` Create plots for all samples
#' @param plots_QC `logical(1)` Create plots for all QCs
#' @param diagnostic_plots `logical(1)` Create diagnostic plots of 5 QCs
#'     spread across the runs
#' @param batch_positions `list` Indicate start and end file of each batch,
#'     e.g. `list(c(1,20),c(21,40))`
#' @param QC_pattern `character(1)`  Pattern of QC files
#' @param sample_pattern `character(1)` Pattern of sample files
#' @param rt_alignment `logical(1)` Align retention time based on internal
#'     standard compounds in the QC samples.
#' @param int_std_id `character` Provide ID's of internal standard compounds for
#'     retention time alignment
#' @param screening_mode `logical(1)` Run the algorithm over 5 QCs to quickly
#'     check retention time shifts
#' @param smoothing `logical(1)` Smooth the peaks with [sgolayfilt()]
#' @param max_int_filter `numeric(1)` Disregard peaks with a max. int. lower
#'     than this value
#'
#' @import MsExperiment
#' @importFrom Spectra MsBackendMzR
#' @importFrom Spectra MsBackendSql
#' @importFrom Spectra filterMzRange
#' @importFrom Spectra filterEmptySpectra
#' @importFrom Spectra filterDataOrigin
#' @importFrom Spectra filterRt
#' @importFrom Spectra filterPolarity
#' @importFrom Spectra dataOrigin
#' @importFrom Spectra addProcessing
#' @importFrom Spectra Spectra
#' @importFrom signal sgolayfilt
#' @importFrom xcms PeakGroupsParam
#' @importFrom xcms adjustRtime
#' @importFrom xcms applyAdjustedRtime
#' @importFrom xcms rtime
#' @importFrom xcms intensity
#' @importFrom pracma trapz
#' @importFrom BiocParallel SnowParam
#' @importFrom tidyr spread
#' @importFrom writexl write_xlsx
#' @importFrom dplyr summarise
#' @importFrom dplyr summarise_at
#' @importFrom dplyr group_by
#' @importFrom dplyr select
#' @importFrom dplyr first
#' @importFrom S4Vectors DataFrame
#'
#' @return returns `list` with auc table and feature table with summarized
#'     stats per compound. Outputs plots and other tables to output folder.
#'
#' @export
#'

smoothingSG <- function(tr,
                        all_files_i,
                        spectra_QC,
                        rt_input,
                        mz_input,
                        smoothing) {
  
  eic <- filterSingle_extractEIC(spectra_QC, all_files_i, rt_input, mz_input)
  
  rt  <- eic[, 1L]
  int <- eic[, 2L]
  
  int[which(is.na(int))] <- 0
  
  if (length(int) < 7) {
    if (length(int) %% 2 == 0) {
      fl <- length(int) - 1
    } else {
      fl <- length(int)
    }
  } else {
    fl <- 7
  }
  smoothed <- sgolayfilt(int, p = 3, n = fl)
  if (smoothing == TRUE) {
    int <- smoothed
    int[int < 0] <- 0
  }
  # Border detection
  border <-
    find_peak_points(rt, smoothed, tr,
                     .check = FALSE
    )  
  return(list(rt, int, border))
}

# a function to check if a peak is valid
# checkValidPeak <- function(x, y, rt, int, border, sample_name, compound_info, results) {
# 
#   if (length(unique(y)) > 1) {
#     auc <- trapz(x, y)
#     pop <- length(x)
#     qscore <- qscoreCalculator(x, y)
#     #compound_info <- dbData[j, ]
#     results <- rbind(
#       results,
#       data.frame(
#         Component = compound_info$ID,
#         Sample = sample_name,
#         AUC = auc,
#         MaxInt = int[border[3L]],
#         SNR = qscore[1],
#         peak_cor = qscore[2],
#         foundRT = rt[border[3L]],
#         pop = pop,
#         compound_info
#       )
#     )
#   } else {
#     #compound_info <- dbData[j, ]
#     
#     results <- rbind(
#       results,
#       data.frame(
#         Component = compound_info$ID,
#         Sample = sample_name,
#         AUC = NA,
#         MaxInt = NA,
#         SNR = NA,
#         peak_cor = NA,
#         foundRT = NA,
#         pop = NA,
#         compound_info
#       )
#     )
#   }
#   return(results)
# }

# a function to check if a peak is valid
checkValidPeak <- function(x, y, d, sample_name, int, rt, border) {
  # if there is NA, impute via Linear Interpolation
  if (anyNA(x)) {
    x <- imputeLinInterpol(x)
  }
  if (anyNA(y)) {
    y <- imputeLinInterpol(y)
  }
  # Force consistent scalar outputs
  if (length(unique(y)) > 1 && !any(is.na(border))) {
    
    auc <- pracma::trapz(x, y)
    pop <- length(x)
    qscore <- qscoreCalculator(x, y)
    
    MaxInt_val <- ifelse(length(int) >= border[3L], int[border[3L]], NA)
    foundRT_val <- ifelse(length(rt) >= border[3L], rt[border[3L]], NA)
    
    out <- data.frame(
      Component = as.character(d$ID),
      Sample = as.character(sample_name),
      AUC = as.numeric(auc),
      MaxInt = as.numeric(MaxInt_val),
      SNR = as.numeric(qscore[1]),
      peak_cor = as.numeric(qscore[2]),
      foundRT = as.numeric(foundRT_val),
      pop = as.numeric(pop),
      ID = as.character(d$ID),
      NAME = as.character(d$NAME),
      mz = as.numeric(d$`m/z`),
      tr = as.numeric(d$tr),
      stringsAsFactors = FALSE
    )
    
  } else {
    
    out <- data.frame(
      Component = as.character(d$ID),
      Sample = as.character(sample_name),
      AUC = NA_real_,
      MaxInt = NA_real_,
      SNR = NA_real_,
      peak_cor = NA_real_,
      foundRT = NA_real_,
      pop = NA_real_,
      ID = as.character(d$ID),
      NAME = as.character(d$NAME),
      mz = as.numeric(d$`m/z`),
      tr = as.numeric(d$tr),
      stringsAsFactors = FALSE
    )
  }
  
  colnames(out) <- c("Component", "Sample", "AUC", "MaxInt",
                     "SNR", "peak_cor", "foundRT", "pop", "ID", "NAME", "mz", "tr") # Force the names here
  
  return(out)
}


# function for data handling
dataHandling <- function(files, string, QC_pattern, polarity){
  data <- MsExperiment()
  experimentFiles(data) <-
    MsExperimentFiles(mzML = setNames(files, basename(
      tools::file_path_sans_ext(files)
    )))
  
  sampleData(data) <- DataFrame(sample_index = 1:length(files),
                                spectraOrigin = files)
  if (string == "QC_files"){
    sampleData(data)$type <- QC_pattern
  }
  if (string == "files_batch"){
    # Define study and QC samples --> all not QC files are deemed study files
    sampleData(data)$type <- "study"
    sampleData(data)$type[grep(pattern = QC_pattern, files)] <- QC_pattern
  }
  
  # --- CHANGED: build a temp SQLite DB, populate it once, use OfflineSql ---
  db_path <- tempfile(fileext = ".sqlite")
  sp_init <- Spectra(
    experimentFiles(data)[["mzML"]],
    backend  = MsBackendMzR(),
    BPPARAM  = SnowParam(workers = 1L)
  )
  sp <- setBackend(sp_init,
                   MsBackendOfflineSql(),
                   drv    = SQLite(),          # ← driver passed to setBackend
                   dbname = db_path,           # ← dbname passed to setBackend
                   BPPARAM = SnowParam(workers = 1L))
  
  # sp <- Spectra(
  #   experimentFiles(data)[["mzML"]],
  #   backend  = MsBackendMzR()
  #   #BPPARAM  = SnowParam(workers = 1L)
  # )
  if (polarity == "positive") {
    spectra(data) <- filterPolarity(sp, 1)
  } else if (polarity == "negative") {
    spectra(data) <- filterPolarity(sp, 0)
  }
  
  sampleData(data)$raw_file <- normalizePath(files)
  data <- linkSampleData(data, with = "sampleData.raw_file = spectra.dataOrigin")
  
  return (data)
}

# edited null handling (leading to null plots)
standardize_results <- function(df) {
  required_cols <- c(
    "Component", "Sample", "AUC", "MaxInt",
    "SNR", "peak_cor", "foundRT", "pop",
    "ID", "NAME", "mz", "tr"
  )
  
  # 1. Ensure it's a data frame
  df <- as.data.frame(df)
  
  # 2. Add missing columns with correct types to avoid logical NAs in numeric cols
  missing <- setdiff(required_cols, colnames(df))
  for (m in missing) df[[m]] <- NA
  
  # 3. Handle list-columns safely without changing row count
  df[] <- lapply(df, function(x) {
    if (is.list(x)) {
      # Replace empty list elements with NA, take only 1st element of others
      sapply(x, function(el) if (length(el) == 0) NA else el[1])
    } else {
      x
    }
  })
  
  # 4. Strict selection and ordering (Avoids manual renaming later)
  df <- df[, required_cols, drop = FALSE]
  
  colnames(df) <- c("Component", "Sample", "AUC", "MaxInt",
                    "SNR", "peak_cor", "foundRT", "pop", "ID", "NAME", "mz", "tr") # Force the names here
  
  return(df)
}

# edited null handling (leading to null plots)
safe_bind <- function(x) {
  # 1. Filter out non-dataframes or empty ones
  x <- Filter(function(df) is.data.frame(df) && nrow(df) > 0, x)
  if (length(x) == 0) return(NULL)
  
  # 2. Standardize each DF before binding
  x <- lapply(x, standardize_results)
  
  # 3. Combine - since they are now identical in structure, this is safe
  out <- dplyr::bind_rows(x)
  
  # 4. Final Type Conversion (Ensure numeric columns aren't characters)
  numeric_cols <- c("AUC", "MaxInt", "SNR", "peak_cor", "foundRT", "mz", "tr")
  out[numeric_cols] <- lapply(out[numeric_cols], as.numeric)
  
  rownames(out) <- NULL
  return(out)
}


## jo: wouldn't it be better to call the function on a data object instead
## of a file path? The (advanced) user could eventually do some more quality
## checks on the data before?
## Pablo: Definitely something I should do, but I might do it a separate
## function, since for use with the GUI the file input is nice.

tardisPeaks <-
    function(file_path = NULL,
    lcmsData = NULL,
    dbData,
    ppm = 5,
    rtdev = 18,
    mass_range = NULL,
    polarity = "positive",
    output_directory,
    plots_samples = FALSE,
    plots_QC = FALSE,
    diagnostic_plots = TRUE,
    batch_positions,
    QC_pattern = "QC",
    sample_pattern = "",
    rt_alignment = TRUE,
    int_std_id,
    screening_mode = FALSE,
    smoothing = TRUE,
    max_int_filter = NULL) {
        results_samples <-
            data.frame(
                Component = character(0),
                Sample = character(0),
                AUC = numeric(0),
                SNR = numeric(0),
                peak_cor = numeric(0),
                foundRT = numeric(0),
                pop = numeric(0)
            )
        results_QCs <-
            data.frame(
                Component = character(0),
                Sample = character(0),
                AUC = numeric(0),
                SNR = numeric(0),
                peak_cor = numeric(0),
                foundRT = numeric(0),
                pop = numeric(0)
            )
        if (is.null(file_path) == FALSE) {
            files <-
                list.files(file_path, full.names = T, pattern = "mzML|mzXML")
        }
        if (is.null(lcmsData) == FALSE) {
            if (polarity == "positive") {
                spectra(lcmsData) <- filterPolarity(spectra(lcmsData), 1)
            } else if (polarity == "negative") {
                spectra(lcmsData) <- filterPolarity(spectra(lcmsData), 0)
            }
            suppressWarnings(
                lcmsData <- linkSampleData(
                    lcmsData,
                    with = "sampleData.raw_file = spectra.dataOrigin"
                )
            )
        }
        if (is.null(mass_range) == FALSE) {
            dbData <-
                dbData[which(dbData$`m/z` < mass_range[2] &
                    dbData$`m/z` > mass_range[1]), ]
        }
        info_compounds <- dbData
        if (screening_mode == TRUE) {
            results_screening <-
                data.frame(
                    Component = character(0),
                    Sample = character(0),
                    AUC = numeric(0),
                    SNR = numeric(0),
                    peak_cor = numeric(0),
                    foundRT = numeric(0),
                    pop = numeric(0)
                )
            if (is.null(file_path) == FALSE) {
                QC_files <-
                    files[grep(pattern = QC_pattern, files)]
                data_QC <- dataHandling(QC_files, "QC_files", QC_pattern, polarity)
            } else {
                data_QC <- lcmsData[which(sampleData(lcmsData)$type == QC_pattern)]
            }
            if (is.null(mass_range) == FALSE) {
                data_QC <- filterSpectra(data_QC, filterMzRange, mz = mass_range) |>
                    filterSpectra(filterEmptySpectra)
            } else {
                data_QC <- data_QC
            }
            spectra_QC <- data_QC@spectra
            checkScans(spectra_QC)
            data_QC@spectra <- spectra_QC
            ## Create ranges for all compounds
            ranges <- createRanges(data_QC, dbData, ppm, rtdev)
            ## Get mz & rt ranges
            mzRanges <- ranges[[1L]]
            rtRanges <- ranges[[2L]]

            if (rt_alignment == TRUE) {
                ## Get the ranges for the internal standard compounds
                internal_standards_rt <-
                    rtRanges[which(dbData$ID %in% int_std_id), ]
                internal_standards_mz <-
                    mzRanges[which(dbData$ID %in% int_std_id), ]
                dbData_std <- dbData[which(dbData$ID %in% int_std_id), ]
                ## Get QC sample names
                sample_names <-
                    lapply(data_QC@sampleData$spectraOrigin, basename)
                ## Initiate empty vectors
                int_std_foundrt <- c()
                int_std <- c()
                ## Retrieve foundRT of internal standards in QC's,
                ## loop over all samples and all internal standards
                # for (j in 1:dim(internal_standards_rt)[1]) {
                #     rt_list <- list()
                #     int_list <- list()
                #     x_list <- list()
                #     y_list <- list()
                #     for (i in 1:length(sample_names)) {
                #         sample_name <- unlist(sample_names[i])
                #         
                #         res <- smoothingSG(
                #           dbData_std$tr[j],
                #           unique(dataOrigin(spectra_QC))[i],
                #           spectra_QC,
                #           internal_standards_rt[j, ],
                #           internal_standards_mz[j, ],
                #           smoothing
                #         )
                #         rt <- res$rt
                #         int <- res$int
                #         border <- res$border
                # 
                #         # Save found RT for internal standard target
                #         int_std_foundrt <-
                #             cbind(int_std_foundrt, rt[border[3L]]) # this will finally contain
                #         # all found rts for the different internal standards in this sample
                #     }
                # 
                #     ## this will contain all the found rts in all samples
                #     int_std <-
                #         rbind(int_std, int_std_foundrt)
                #     int_std_foundrt <- c()
                # }
                
                # lapply code block
                results_list1 <- lapply(1:nrow(internal_standards_rt), function(j) {
                  local_found_rt <- numeric(length(sample_names))
                  for (i in 1:length(sample_names)) {
                    res <- smoothingSG(
                                dbData_std$tr[j],
                                unique(dataOrigin(spectra_QC))[i],
                                spectra_QC,
                                internal_standards_rt[j, ],
                                internal_standards_mz[j, ],
                                smoothing
                              )
                    local_found_rt[i] <- res$rt[res$border[3L]]  # rt of peak
                  }
                  return(local_found_rt)
                }, cl = cl)   # results_list1: list of lists
                
                int_std <- do.call(rbind, results_list1)
                
                ## Define parameters for retention time adjustment, based on the QC's
                ## and the internal standards
                param <-
                    PeakGroupsParam(
                        minFraction = 0.9,
                        span = 0.5,
                        peakGroupsMatrix = int_std
                    )
                data_QC <- adjustRtime(data_QC, param = param)
                data_QC <- applyAdjustedRtime(data_QC)
            }
            # Find all targets in x QC's

            spectra_QC <- data_QC@spectra
            sample_names <-
                lapply(data_QC@sampleData$spectraOrigin, basename)

            # for (j in 1:dim(rtRanges)[1]) {
            #     rt_list <- list()
            #     int_list <- list()
            #     x_list <- list()
            #     y_list <- list()
            #     for (i in 1:length(sample_names)) {
            #         sample_name <- unlist(sample_names[i])
            #         
            #         res <- smoothingSG(
            #           dbData$tr[j],
            #           unique(dataOrigin(spectra_QC))[i],
            #           spectra_QC,
            #           rtRanges[j, ],
            #           mzRanges[j, ],
            #           smoothing
            #         )
            #         rt <- res$rt
            #         int <- res$int
            #         border <- res$border
            #         
            #         idx <- border[1L]:border[2L]
            #         x <- rt[idx]
            #         y <- int[idx]
            #         rt_list <- c(rt_list, list(rt))
            #         int_list <- c(int_list, list(int))
            #         x_list <- c(x_list, list(x))
            #         y_list <- c(y_list, list(y))
            #         
            #         results_screening <- checkValidPeak(x, y, rt, int, border, sample_name, dbData[j, ], results_screening)
            #     }
            #     # Create and save the plot for the current component
            #     batchnr <- 1
            #     if (diagnostic_plots == TRUE) {
            #         plotDiagnostic(
            #             compound_info,
            #             output_directory,
            #             rt_list,
            #             int_list,
            #             x_list,
            #             y_list,
            #             batchnr,
            #             sample_names
            #         )
            #     }
            # }
            
            # lapply code block (j = internal standards)
            results_list2 <- lapply(1:dim(rtRanges)[1], function(j) {  # for all target compounds
              compound_info <- dbData[j, ]  # id, name, mz, rt
              rt_list <- vector("list", length(sample_names))
              int_list <- vector("list", length(sample_names))
              x_list <- vector("list", length(sample_names))
              y_list <- vector("list", length(sample_names))
              results_screening_row <- vector("list", length(sample_names))
              #output_cache <- vector("list", length(dim(rtRanges)[1]))
              
              compound_results <- list() # store each sample's results
              for (i in 1:length(sample_names)) {
                res <- smoothingSG(
                            dbData$tr[j],
                            unique(dataOrigin(spectra_QC))[i],
                            spectra_QC,
                            rtRanges[j, ],
                            mzRanges[j, ],
                            smoothing
                          )
                rt <- if(is.null(res$rt)) NULL else res$rt
                int <- if(is.null(res$int)) NULL else res$int
                border <- if(is.null(res$border)) NULL else res$border
                
                idx <- border[1L]:border[2L]
                x <- rt[idx]
                y <- int[idx]
                
                rt_list[[i]] <- rt
                int_list[[i]] <- int
                x_list[[i]] <- x
                y_list[[i]] <- y
                
                results_screening_row[[i]] <- checkValidPeak(x,
                                                             y,
                                                             compound_info,
                                                             sample_names[i],
                                                             int,
                                                             rt,
                                                             border)
              }
              #Create and save the plot for the current component
              batchnr <- 1
              if (diagnostic_plots == TRUE) {
                plotDiagnostic(
                  compound_info,
                  output_directory,
                  rt_list,
                  int_list,
                  x_list,
                  y_list,
                  batchnr,
                  sample_names
                  #QC_pattern
                )
              }
              
              return(safe_bind(results_screening_row))
            }, cl=cl)
            
            # Combine and standardize result
            results_screening <- safe_bind(results_list2)
            results_screening <- standardize_results(results_screening)
            
            avg_metrics_table <- results_screening %>%
                group_by(Component) %>%
                summarise_at(vars(-Sample), list(~ if (is.numeric(.)) {
                    mean(., na.rm = TRUE)
                } else {
                    first(.)
                }))
            write.csv(avg_metrics_table,
                file = paste0(output_directory, "qc_screening.csv")
            )
        } else {
            ## Loop over the batches
            for (batchnr in 1:length(batch_positions)) {
                dbData <- info_compounds # need to reset? better to keep updated from
                # last batch?
                results_QCs_batch <-
                    data.frame(
                        Component = character(0),
                        Sample = character(0),
                        AUC = numeric(0),
                        SNR = numeric(0),
                        peak_cor = numeric(0),
                        foundRT = numeric(0),
                        pop = numeric(0)
                    )
                if (is.null(file_path) == FALSE) {
                    files_batch <-
                        files[batch_positions[[batchnr]][1]:batch_positions[[batchnr]][2]]
                    data_batch <- dataHandling(files_batch, "files_batch", QC_pattern, polarity)
                } else {
                    data_batch <- lcmsData[batch_positions[[batchnr]][1]:batch_positions[[batchnr]][2]]
                }
                if (is.null(mass_range) == FALSE) {
                    data_batch <- filterSpectra(data_batch, filterMzRange, mz = mass_range) |>
                        filterSpectra(filterEmptySpectra)
                } else {
                    data_batch <- data_batch
                }
                spectra_batch <- data_batch@spectra
                checkScans(spectra_batch)
                data_batch@spectra <- spectra_batch

                data_QC <-
                    data_batch[which(sampleData(data_batch)$type == "QC")]
                if (is.null(mass_range == FALSE)) {
                    spectra_QC <- data_QC@spectra |>
                        filterMzRange(mass_range) |>
                        filterEmptySpectra()
                } else {
                    spectra_QC <- data_QC@spectra
                }
                ranges <- createRanges(data_QC, dbData, ppm, rtdev)
                mzRanges <- ranges[[1]]
                rtRanges <- ranges[[2]]
                if (rt_alignment == TRUE) {
                    ## Get the ranges for the internal standard compounds
                    internal_standards_rt <-
                        rtRanges[which(dbData$ID %in% int_std_id), ]
                    internal_standards_mz <-
                        mzRanges[which(dbData$ID %in% int_std_id), ]
                    dbData_std <- dbData[which(dbData$ID %in% int_std_id), ]
                    sample_names <-
                        lapply(data_QC@sampleData$spectraOrigin, basename)
                    int_std_foundrt <- c()
                    int_std <- c()
                    # for (j in 1:dim(internal_standards_rt)[1]) {
                    #     rt_list <- list()
                    #     int_list <- list()
                    #     x_list <- list()
                    #     y_list <- list()
                    #     for (i in 1:length(sample_names)) {
                    #         sample_name <- unlist(sample_names[i])
                    #         
                    #         res <- smoothingSG(
                    #           dbData_std$tr[j],
                    #           unique(dataOrigin(spectra_QC))[i],
                    #           spectra_QC,
                    #           internal_standards_rt[j, ],
                    #           internal_standards_mz[j, ],
                    #           smoothing
                    #         )
                    #         rt <- res$rt
                    #         int <- res$int
                    #         border <- res$border
                    #         
                    #         int_std_foundrt <-
                    #             cbind(int_std_foundrt, rt[border[3L]])
                    #     }
                    #     int_std <-
                    #         rbind(int_std, int_std_foundrt)
                    #     int_std_foundrt <- c()
                    # }
                    # lapply code block
                    results_list3 <- lapply(1:dim(internal_standards_rt)[1], function(j) {
                      local_found_rt <- numeric(length(sample_names))
                      
                      for (i in 1:length(sample_names)) {
                        res <- smoothingSG(
                                    dbData_std$tr[j],
                                    unique(dataOrigin(spectra_QC))[i],
                                    spectra_QC,
                                    internal_standards_rt[j, ],
                                    internal_standards_mz[j, ],
                                    smoothing
                                  )
                        local_found_rt[i] <- if(is.null(res$border) || is.null(res$rt)) NULL else res$rt[res$border[3L]]  # rt of peak
                      }
                      
                      return(local_found_rt)
                    }, cl=cl) # results_list3: list of lists
                    int_std <- do.call(rbind, results_list3)
                    
                    param <-
                        PeakGroupsParam(
                            minFraction = 0.9,
                            span = 0.5,
                            peakGroupsMatrix = int_std,
                            subset = which(sampleData(data_batch)$type == "QC")
                        )
                    data_batch <- adjustRtime(data_batch, param = param)
                    data_batch <- applyAdjustedRtime(data_batch)
                }
                ## Now, we try and find ALL compounds in the QC samples and save their
                ## foundRT to search the compounds at that RT in the sample files
                ## Skip this step if there aren't any QC's available.
                data_QC <-
                    data_batch[which(sampleData(data_batch)$type == "QC")]
                if (length(data_QC) != 0) {
                    sample_names <-
                        lapply(data_QC@sampleData$spectraOrigin, basename)
                    if (is.null(mass_range) == FALSE) {
                        spectra_QC <- data_QC@spectra |>
                            filterMzRange(mass_range) |>
                            filterEmptySpectra()
                    } else {
                        spectra_QC <- data_QC@spectra
                    }
                    # for (j in 1:dim(rtRanges)[1]) {
                    #     rt_list <- list()
                    #     int_list <- list()
                    #     x_list <- list()
                    #     y_list <- list()
                    #     for (i in 1:length(sample_names)) {
                    #         sample_name <- unlist(sample_names[i])
                    #         
                    #         res <- smoothingSG(
                    #           dbData$tr[j],
                    #           unique(dataOrigin(spectra_QC))[i],
                    #           spectra_QC,
                    #           rtRanges[j, ],
                    #           mzRanges[j, ],
                    #           smoothing
                    #         )
                    #         rt <- res$rt
                    #         int <- res$int
                    #         border <- res$border
                    #         
                    #         idx <- border[1L]:border[2L]
                    #         x <- rt[idx]
                    #         y <- int[idx]
                    #         rt_list <- c(rt_list, list(rt))
                    #         int_list <- c(int_list, list(int))
                    #         x_list <- c(x_list, list(x))
                    #         y_list <- c(y_list, list(y))
                    #         
                    #         results_QCs_batch <- checkValidPeak(x, y, rt, int, border, sample_name, dbData[j, ], results_QCs_batch)
                    #     }
                    #     if (plots_QC == TRUE) {
                    #         plotQCs(
                    #             compound_info,
                    #             output_directory,
                    #             rt_list,
                    #             int_list,
                    #             x_list,
                    #             y_list,
                    #             batchnr,
                    #             sample_names
                    #         )
                    #     }
                    # }
                    # lapply code block
                    results_list4 <- lapply(1:dim(rtRanges)[1], function(j) {
                      compound_info <- dbData[j, ]
                      rt_list <- vector("list", length(sample_names))
                      int_list <- vector("list", length(sample_names))
                      x_list <- vector("list", length(sample_names))
                      y_list <- vector("list", length(sample_names))
                      results_QCs_batch_row <- vector("list", length(sample_names))
                      
                      for (i in 1:length(sample_names)) {
                        res <- smoothingSG(
                                    dbData$tr[j],
                                    unique(dataOrigin(spectra_QC))[i],
                                    spectra_QC,
                                    rtRanges[j, ],
                                    mzRanges[j, ],
                                    smoothing
                                  )
                        rt <- if(is.null(res$rt)) NULL else res$rt
                        int <- if(is.null(res$int)) NULL else res$int
                        border <- if(is.null(res$border)) NULL else res$border
                        
                        idx <- border[1L]:border[2L]
                        x <- rt[idx]
                        y <- int[idx]
                        rt_list[[i]] <- rt
                        int_list[[i]] <- int
                        x_list[[i]] <- x
                        y_list[[i]] <- y
                        
                        results_QCs_batch_row[[i]] <- checkValidPeak(x,
                                                                     y,
                                                                     dbData[j, ],
                                                                     sample_names[i],
                                                                     int,
                                                                     rt,
                                                                     border)  # is a 1-row dataframe
                      }
                      # Create and save the plot for the current component
                      batchnr <- 1
                      if (plots_QC == TRUE) {
                        plotQCs(
                          compound_info,
                          output_directory,
                          rt_list,
                          int_list,
                          x_list,
                          y_list,
                          batchnr,
                          sample_names
                        )
                      }
                      return(safe_bind(results_QCs_batch_row))
                    }, cl=cl)
                    
                    # Combine and standardize result
                    results_QCs_batch <- safe_bind(results_list4)
                    results_QCs_batch <- standardize_results(results_QCs_batch)
                    
                    #results_QCs <- rbind(results_QCs, results_QCs_batch)
                    ## Replace rtmed with average foundRT from previous results
                    new_rt_avg <- results_QCs_batch %>%
                        group_by(ID) %>%
                        summarise(mean = mean(foundRT), na.rm = TRUE)
                    dbData <- merge(dbData, new_rt_avg, by = "ID")
                    dbData$trold <- dbData$tr
                    dbData$tr <- new_rt_avg$mean
                    ## If no RT is found, restore old RT
                    for (k in 1:dim(dbData)[1]) {
                        if (is.na(dbData$tr[k]) == TRUE) {
                            dbData$tr[k] <- dbData$trold[k]
                        }
                    }
                }
                ## Next do the whole analysis for the samples in the same batch of the
                ## QC's to find ALL the compounds at the corrected RT. (SAMPLES + QC)
                ## Get sample data
                sample_names <-
                    lapply(data_batch@sampleData$spectraOrigin, basename)
                # Create ranges around new RT
                ranges <- createRanges(data_batch, dbData, ppm, rtdev)
                mzRanges <- ranges[[1]]
                rtRanges <- ranges[[2]]
                if (is.null(mass_range) == FALSE) {
                    spectra <- data_batch@spectra |>
                        filterMzRange(mz = mass_range) |>
                        filterEmptySpectra()
                } else {
                    spectra <- data_batch@spectra
                }

                # for (j in 1:dim(rtRanges)[1]) {
                #     rt_list <- list()
                #     int_list <- list()
                #     x_list <- list()
                #     y_list <- list()
                #     for (i in 1:length(sample_names)) {
                #         sample_name <- unlist(sample_names[i])
                #         
                #         res <- smoothingSG(
                #           dbData$tr[j],
                #           unique(dataOrigin(spectra))[i],
                #           spectra,
                #           rtRanges[j, ],
                #           mzRanges[j, ],
                #           smoothing
                #         )
                #         rt <- res$rt
                #         int <- res$int
                #         border <- res$border
                #         
                #         idx <- border[1L]:border[2L]
                #         x <- rt[idx]
                #         y <- int[idx]
                #         rt_list <- c(rt_list, list(rt))
                #         int_list <- c(int_list, list(int))
                #         x_list <- c(x_list, list(x))
                #         y_list <- c(y_list, list(y))
                # 
                #         results_samples <- checkValidPeak(x, y, rt, int, border, sample_name, dbData[j, ], results_samples)
                #     }
                #     if (plots_samples == TRUE) {
                #         plotSamples(
                #             compound_info,
                #             output_directory,
                #             rt_list,
                #             int_list,
                #             x_list,
                #             y_list,
                #             batchnr,
                #             sample_names
                #         )
                #     }
                #     if (diagnostic_plots == TRUE) {
                #         plotDiagnostic(
                #             compound_info,
                #             output_directory,
                #             rt_list,
                #             int_list,
                #             x_list,
                #             y_list,
                #             batchnr,
                #             sample_names
                #         )
                #     }
                # }
                # lapply code block
                results_list5 <- lapply(1:dim(rtRanges)[1], function(j){
                  compound_info <- dbData[j, ]
                  rt_list <- vector("list", length(sample_names))
                  int_list <- vector("list", length(sample_names))
                  x_list <- vector("list", length(sample_names))
                  y_list <- vector("list", length(sample_names))
                  results_samples_row <- vector("list", length(sample_names))
                  
                  for (i in 1:length(sample_names)) {
                    res <- smoothingSG(
                                dbData$tr[j],
                                unique(dataOrigin(spectra))[i],
                                spectra,
                                rtRanges[j, ],
                                mzRanges[j, ],
                                smoothing
                              )
                    rt <- if(is.null(res$rt)) NULL else res$rt
                    int <- if(is.null(res$int)) NULL else res$int
                    border <- if(is.null(res$border)) NULL else res$border
                    
                    idx <- border[1L]:border[2L]
                    x <- rt[idx]
                    y <- int[idx]
                    rt_list[[i]] <- rt
                    int_list[[i]] <- int
                    x_list[[i]] <- x
                    y_list[[i]] <- y
                    
                    results_samples_row[[i]] <- checkValidPeak(x,
                                                               y,
                                                               dbData[j, ],
                                                               sample_names[i],
                                                               int,
                                                               rt,
                                                               border)
                    
                  }
                  if (plots_samples == TRUE) {
                    plotSamples(
                      compound_info,
                      output_directory,
                      rt_list,
                      int_list,
                      x_list,
                      y_list,
                      batchnr,
                      sample_names
                      #STD_pattern
                    )
                  }
                  if (diagnostic_plots == TRUE) {
                    plotDiagnostic(
                      compound_info,
                      output_directory,
                      rt_list,
                      int_list,
                      x_list,
                      y_list,
                      batchnr,
                      sample_names
                      #QC_pattern
                    )
                  }
                  return (safe_bind(results_samples_row))
                }, cl=cl)
                
                # Combine dataframes
                results_samples <- safe_bind(results_list5)
                results_samples <- standardize_results(results_samples)
            }
            results <- results_samples  # 1 batch, same copy
            results <- standardize_results(results)
            # }
            # results <- results_samples

            if (is.null(max_int_filter) == FALSE && max_int_filter != 0) {
                results <- results[which(results$MaxInt >= max_int_filter), ]
            }

            # AUC, int, SNR & peakcor tables for each component peak in every sample
            auc_table <- results %>%
                select(Component, Sample, AUC) %>%
                spread(Sample, AUC, fill = NA, drop = FALSE)

            write.csv(auc_table, file = paste0(output_directory, "auc_table.csv"))
            pop_table <- results %>%
                select(Component, Sample, pop) %>%
                spread(Sample, pop, fill = NA, drop = FALSE)

            write.csv(pop_table, file = paste0(output_directory, "pop_table.csv"))
            SNR_table <- results %>%
                select(Component, Sample, SNR) %>%
                spread(Sample, SNR, fill = NA, drop = FALSE)

            write.csv(SNR_table, file = paste0(output_directory, "snr_table.csv"))
            int_table <- results %>%
                select(Component, Sample, MaxInt) %>%
                spread(Sample, MaxInt, fill = NA, drop = FALSE)

            write.csv(int_table, file = paste0(output_directory, "int_table.csv"))
            peakcor_table <- results %>%
                select(Component, Sample, peak_cor) %>%
                spread(Sample, peak_cor, fill = NA, drop = FALSE)

            write.csv(peakcor_table, file = paste0(output_directory, "peakcor_table.csv"))


            # summarize feature table based on QC's
            avg_metrics_table <- NULL
            if (length(data_QC) != 0) {
                QC_results <- results[grep("QC", results$Sample), ]
                avg_metrics_table <- QC_results %>%
                    group_by(Component) %>%
                    summarise_at(vars(-Sample), list(~ if (is.numeric(.)) {
                        mean(., na.rm = TRUE)
                    } else {
                        first(.)
                    }))
                write_xlsx(avg_metrics_table, paste0(
                    output_directory,
                    "feat_table.xlsx"
                ))
            }

            # save input parameters to .csv

            input_params <- data.frame(
                "ppm" = .collapse_safe(ppm),
                "rtdev" = .collapse_safe(rtdev),
                "mass_range_low" = .collapse_safe(mass_range[1]),
                "mass_range_high" = .collapse_safe(mass_range[2]),
                "polarity" = .collapse_safe(polarity),
                "batch_positions" = .collapse_safe(batch_positions),
                "QC_pattern" = .collapse_safe(QC_pattern),
                "sample_pattern" = .collapse_safe(sample_pattern),
                "int_std_id" = .collapse_safe(int_std_id),
                "screening_mode" = .collapse_safe(screening_mode),
                "rt_alignment" = .collapse_safe(rt_alignment),
                "plots_samples" = .collapse_safe(plots_samples),
                "plots_QC" = .collapse_safe(plots_QC),
                "diagnostic_plots" = .collapse_safe(diagnostic_plots),
                "max_int_filter" = .collapse_safe(max_int_filter),
                "smoothing" = .collapse_safe(smoothing),
                stringsAsFactors = FALSE
            )

            write.csv(t(input_params),
                file = paste0(output_directory, "input_params.csv"),
                row.names = TRUE
            )

            return(list(auc_table, avg_metrics_table))
        }
    }
