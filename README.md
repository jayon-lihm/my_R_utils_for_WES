# R utils for WES

A set of scripts/functions that I use frequently for tumor/normal paired WES analysis

## Installing required packages

Required CRAN and Bioconductor packages are listed in `R_packages_needed.txt`. From the repository root, run:

```r
source("package_install.R")
```

This installs any packages not already present (via CRAN, falling back to Bioconductor for packages not on CRAN).
