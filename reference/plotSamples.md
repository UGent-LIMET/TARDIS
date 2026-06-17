# Function to make plots of samples

Function to make plots of samples

## Usage

``` r
plotSamples(
  compound_info,
  output_directory,
  rt_list,
  int_list,
  x_list,
  y_list,
  batchnr,
  sample_names
)
```

## Arguments

- compound_info:

  info for target compounds like name, id, ...

- output_directory:

  `character(1)` directory to save plots.

- rt_list:

  `list` list with retention times of complete window.

- int_list:

  `list` list with intensities in complete window.

- x_list:

  `list` list with retention times of integrated peak.

- y_list:

  `list` list with intensities of integrated peak.

- batchnr:

  `numeric(1)` batch number.

- sample_names:

  names of samples.

## Author

Pablo Vangeenderhuysen
