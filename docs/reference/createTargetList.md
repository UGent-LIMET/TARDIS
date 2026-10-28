# Create target list

Loads a .xlsx or .csv file that contains all the information needed to
detect targeted compounds.

## Usage

``` r
createTargetList(
  input_directory_targets,
  pos_pattern,
  neg_pattern,
  polarity,
  ion_column,
  columns_of_interest
)
```

## Arguments

- input_directory_targets:

  Full path to the file containing the database with targeted compounds.

- pos_pattern:

  Pattern of positive ions to be recognized in ion_column

- neg_pattern:

  Pattern of negative ions to be recognized in ion_column

- polarity:

  Indicates the polarity, can be either "negative" or "positive"

- ion_column:

  Name of the column in which the type of ionisation is found

- columns_of_interest:

  Column names of the columns with information that need to be included,
  can be as many as needed, but needs to include "ID",
  "Name","m/z-value" and "RT (min)"

## Value

A dataframe containing the targeted compounds

## Author

Pablo Vangeenderhuysen
