## Reads R_packages_needed.txt and installs any packages not already installed.
## CRAN packages are installed via install.packages(); if that fails
## (e.g. the package is only on Bioconductor), it falls back to BiocManager::install().

if (!file.exists("R_packages_needed.txt")) {
  stop("R_packages_needed.txt not found in the current working directory")
}

pkgs <- readLines("R_packages_needed.txt")
pkgs <- sub("#.*", "", pkgs)          # strip comments
pkgs <- trimws(pkgs)
pkgs <- pkgs[pkgs != ""]              # drop blank lines

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

for (pkg in pkgs) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    cat(pkg, "already installed, skipping.\n")
    next
  }
  cat("Installing", pkg, "...\n")
  result <- tryCatch({
    install.packages(pkg)
    requireNamespace(pkg, quietly = TRUE)
  }, error = function(e) FALSE)

  if (!isTRUE(result)) {
    cat(pkg, "not found on CRAN, trying Bioconductor...\n")
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
  }
}
