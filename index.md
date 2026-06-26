# TARDIS ![](reference/figures/tardis.png)

R package for *TArgeted Raw Data Integration In Spectrometry*

## Installation

Make sure `R` (**version \>= 4.4.0**) is installed on your computer:

<https://cloud.r-project.org/index.html>

For **Windows** users, `Rtools` should be installed as well:

<https://cran.r-project.org/bin/windows/Rtools/rtools44/rtools.html>

For **Mac** users, please install Xcode and a GNU Fortran compiler, see:
<https://mac.r-project.org/tools/index.html> For Xcode, you can run this
line in your Mac OS terminal:

    sudo xcode-select --install

For the latest version of `TARDIS`, `BiocManager` **version 3.20** is
required: <https://www.bioconductor.org/install/>

In R, run:

    if (!require("devtools", quietly = TRUE))
        install.packages("devtools")

    if (!require("BiocManager", quietly = TRUE))
        install.packages("BiocManager")

    BiocManager::install("UGent-LIMET/TARDIS",build_vignettes = FALSE,
      dependencies=TRUE)

Since the package contains some example data, connection timeout is
often an issue. You can increase your timeout setting in R using:

    options(timeout = "300")

The number indicates the timeout setting in seconds.

## Usage

Load the package:

    library(TARDIS)

**Read the vignettes** for a tutorial on how to use `TARDIS`

Available online: <https://pablovgd.github.io/TARDIS/>

Or, if you installed with `build_vignettes = TRUE`

    RShowDoc("quick_start",package = "TARDIS")

    RShowDoc("gui_tutorial",package = "TARDIS")

    RShowDoc("case_study",package = "TARDIS")

To launch the GUI in R:

    runTardis()

## Copyright

`TARDIS` is licensed under the
[GPLv3](http://choosealicense.com/licenses/gpl-3.0/)

As a summary, the GPLv3 license requires attribution, inclusion of
copyright and license information, disclosure of source code and
changes. Derivative work must be available under the same terms.

2.  Pablo Vangeenderhuysen (2024)
