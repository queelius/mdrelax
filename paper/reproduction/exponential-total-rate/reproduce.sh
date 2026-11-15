#!/bin/sh
# Run from the repository root. Installs this checkout into a fresh R library.
set -eu
test -f DESCRIPTION
test -f paper/reproduction/exponential-total-rate/verify.R
task_library=$(mktemp -d)
trap 'rm -rf "$task_library"' EXIT HUP INT TERM
R CMD INSTALL --library="$task_library" .
R_LIBS_USER="$task_library" Rscript --vanilla paper/reproduction/exponential-total-rate/verify.R "$@"
R_LIBS_USER="$task_library" Rscript --vanilla paper/reproduction/exponential-total-rate/negative-control.R
