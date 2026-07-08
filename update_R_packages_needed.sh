#!/usr/bin/env bash
# Scans the repo's R scripts for packages used via library()/require()/
# requireNamespace() calls and pkg:: / pkg::: namespace references, then
# compares them against R_packages_needed.txt. Any package that is used
# but not yet listed gets appended to R_packages_needed.txt under a
# "newly detected" section for the user to file under CRAN/Bioconductor.
set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PKG_FILE="$REPO_DIR/R_packages_needed.txt"

if [ ! -f "$PKG_FILE" ]; then
  echo "Error: $PKG_FILE not found." >&2
  exit 1
fi

# Packages that ship with base R (or are "recommended" packages bundled
# with R) - these never need to be installed, so don't flag them.
BASE_PKGS="base compiler datasets grDevices graphics grid methods parallel splines stats stats4 tools utils tcltk translations MASS boot class cluster codetools foreign KernSmooth lattice mgcv nlme nnet rpart spatial survival Matrix"

R_FILES=()
while IFS= read -r -d '' f; do
  R_FILES+=("$f")
done < <(find "$REPO_DIR" -maxdepth 3 -name "*.R" ! -name "package_install.R" -print0)

if [ "${#R_FILES[@]}" -eq 0 ]; then
  echo "No R scripts found under $REPO_DIR."
  exit 0
fi

USED_PKGS=$(perl -ne '
  while (/\b(?:library|require)\s*\(\s*["\x27]?([A-Za-z][A-Za-z0-9._]*)["\x27]?/g) { print "$1\n"; }
  while (/\brequireNamespace\s*\(\s*["\x27]([A-Za-z][A-Za-z0-9._]*)["\x27]/g) { print "$1\n"; }
  while (/\b([A-Za-z][A-Za-z0-9._]*):::?/g) { print "$1\n"; }
' "${R_FILES[@]}" | sort -u)

FILTERED=$(comm -23 <(echo "$USED_PKGS") <(echo "$BASE_PKGS" | tr ' ' '\n' | sort -u))

KNOWN=$(grep -v '^[[:space:]]*#' "$PKG_FILE" | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//' -e '/^$/d' | sort -u)

NEW_PKGS=$(comm -23 <(echo "$FILTERED") <(echo "$KNOWN"))

if [ -z "$NEW_PKGS" ]; then
  echo "R_packages_needed.txt is up to date - no new packages found in the R scripts."
  exit 0
fi

echo "Found package(s) used in the R scripts but missing from R_packages_needed.txt:"
echo "$NEW_PKGS" | sed 's/^/  - /'

{
  echo ""
  echo "# Newly detected packages ($(date +%F)) - please move each into the CRAN or Bioconductor section above"
  echo "$NEW_PKGS"
} >> "$PKG_FILE"

echo "Appended the package(s) above to $PKG_FILE."
