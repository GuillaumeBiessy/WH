#!/usr/bin/env bash
# calibrer_nosuggests.sh -- proves that _R_CHECK_DEPENDS_ONLY_=true takes effect under a given R copy.
#   bash calibrer_nosuggests.sh <R home> <scratch dir>
# Builds a witness package whose example calls a package of its Suggests (spelling), then checks it twice.
# Expected (measured 2026-10-25, R 4.5.3 Windows): without the variable "Status: OK"; with it "1 ERROR",
# "there is no package called 'spelling'". By design R keeps testthat / RUnit / tinytest for the tests and the
# VignetteBuilder (tools:::setRlibs), so a package's own tests passing under noSuggests proves nothing on its own.
set -u
RHOME="$1"; S="$2"; R="$RHOME/bin/R.exe"
rm -rf "$S/temoinsug" "$S/sans" "$S/avec"; mkdir -p "$S/temoinsug/R" "$S/temoinsug/man" "$S/sans" "$S/avec"
cd "$S"
cat > temoinsug/DESCRIPTION <<'EOF'
Package: temoinsug
Version: 0.1
Title: Witness
Description: Witness package.
Author: G
Maintainer: G <g@example.org>
License: GPL-3
Suggests: spelling
EOF
echo 'export(f)' > temoinsug/NAMESPACE
echo 'f <- function() 1' > temoinsug/R/f.R
cat > temoinsug/man/f.Rd <<'EOF'
\name{f}
\alias{f}
\title{F}
\description{F.}
\usage{f()}
\examples{
f()
spelling::spell_check_text("hello")
}
EOF
"$R" CMD build temoinsug > build.log 2>&1
"$R" CMD check --no-manual --output=sans temoinsug_0.1.tar.gz > sans.log 2>&1
_R_CHECK_DEPENDS_ONLY_=true "$R" CMD check --no-manual --output=avec temoinsug_0.1.tar.gz > avec.log 2>&1
SANS=$(grep '^Status' sans/temoinsug.Rcheck/00check.log); AVEC=$(grep '^Status' avec/temoinsug.Rcheck/00check.log)
echo "sans: $SANS | avec: $AVEC"
[ "$SANS" = "Status: OK" ] && grep -q "no package called 'spelling'" avec/temoinsug.Rcheck/00check.log \
	&& echo "CALIBRATION_NOSUGGESTS=OK" || { echo "CALIBRATION_NOSUGGESTS=FAIL"; exit 1; }
