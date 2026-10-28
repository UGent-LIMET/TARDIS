# Rename Thermo .raw files Script to read XCalibur sequence list and rename files according to run type (sample, QC, blank, etc...)

Rename Thermo .raw files Script to read XCalibur sequence list and
rename files according to run type (sample, QC, blank, etc...)

## Usage

``` r
renameRawFiles(data_path_raw_files, data_path_list)
```

## Arguments

- data_path_raw_files:

  path to raw files. Important: all runs from the sequence list need to
  have a corresponding .raw file in your input folder.

- data_path_list:

  Path to exported sequence list in .csv format, all columns may be
  exported

## Author

Pablo Vangeenderhuysen
