# TARDIS <img src="man/figures/tardis.png" width="150" height="150" align = right />        

R package for *TArgeted Raw Data Integration In Spectrometry*

## Installation
Make sure `R` (**version >= 4.4.0**) is installed on your computer:

https://cloud.r-project.org/index.html

For **Windows** users, `Rtools` should be installed as well:

https://cran.r-project.org/bin/windows/Rtools/rtools44/rtools.html

For **Mac** users, please install Xcode and a GNU Fortran compiler, see:
https://mac.r-project.org/tools/index.html
For Xcode, you can run this line in your Mac OS terminal:

```
sudo xcode-select --install
```

In R, run:

```
pak::pak("UGent-LIMET/TARDIS")
```

If asked do you want to continue, type y:

```
? Do you want to continue (Y/n) y
```

Since the package contains some example data, connection timeout is often an 
issue. You can increase your timeout setting in R using:

```
options(timeout = "300")

```

The number indicates the timeout setting in seconds. 

## Usage

Load the package:

```
library(TARDIS)
```


**Read the vignettes** for a tutorial on how to use `TARDIS`

Available online: https://ugent-limet.github.io/TARDIS/

Or, if you installed with `build_vignettes = TRUE`

```
RShowDoc("quick_start",package = "TARDIS")
```

```
RShowDoc("gui_tutorial",package = "TARDIS")
```

```
RShowDoc("case_study",package = "TARDIS")
```

To launch the GUI in R:

```
runTardis()
```

## Copyright

`TARDIS` is licensed under the [GPLv3](http://choosealicense.com/licenses/gpl-3.0/)

As a summary, the GPLv3 license requires attribution, inclusion of copyright and license information, disclosure of source code and changes. Derivative work must be available under the same terms.