# TARDIS Peak Detection

Main function of the TARDIS package that is called in the Shiny app.
Given data files and a list of targeted compounds it returns the area of
those peaks, optional diagnostic plots and several other parameters. See
vignette for a detailed tutorial.

## Usage

``` r
tardisPeaks(
  file_path = NULL,
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
  max_int_filter = NULL
)
```

## Arguments

- file_path:

  `character(1)` Path to the .mzML or .mzXML files containing LC-MS
  data.

- lcmsData:

  `MsExperiment` MsExperiment containing the data to be preprocessed.
  Sampledata should at least include run type and should match the later
  provided "QC" or "sample" pattern.

- dbData:

  Output of
  [`createTargetList()`](https://pablovgd.github.io/TARDIS/reference/createTargetList.md)

- ppm:

  `numeric(1)` Allowed deviance from given m/z of targets in ppm.

- rtdev:

  `numeric(1)` Allowed deviance from given retention time of compound,
  defines search window for the peak picking algorithm.

- mass_range:

  `numeric(2)` If the user uses data with overlapping mass windows, only
  one mass window at the time can be analyzed, specify this here.

- polarity:

  `character(1)` Ionisation mode to be considered, can be either
  "positive" or "negative"

- output_directory:

  `character(1)` Provide directory to store output

- plots_samples:

  `logical(1)` Create plots for all samples

- plots_QC:

  `logical(1)` Create plots for all QCs

- diagnostic_plots:

  `logical(1)` Create diagnostic plots of 5 QCs spread across the runs

- batch_positions:

  `list` Indicate start and end file of each batch, e.g.
  `list(c(1,20),c(21,40))`

- QC_pattern:

  `character(1)` Pattern of QC files

- sample_pattern:

  `character(1)` Pattern of sample files

- rt_alignment:

  `logical(1)` Align retention time based on internal standard compounds
  in the QC samples.

- int_std_id:

  `character` Provide ID's of internal standard compounds for retention
  time alignment

- screening_mode:

  `logical(1)` Run the algorithm over 5 QCs to quickly check retention
  time shifts

- smoothing:

  `logical(1)` Smooth the peaks with `sgolayfilt()`

- max_int_filter:

  `numeric(1)` Disregard peaks with a max. int. lower than this value

## Value

returns `list` with auc table and feature table with summarized stats
per compound. Outputs plots and other tables to output folder.
