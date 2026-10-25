#!/usr/bin/env bash
# run_windows.sh -- the Windows bench: builds the tarballs, calibrates the instrument, checks every configuration.
#   bash run_windows.sh <repo> <outdir> [<ref of the version to check>] [<calibration ref>]
# Defaults: WH as on CRAN (origin/master, package files identical to the 2.0.0 resubmission) and WH 1.0.6 (c413689),
# the version that failed CRAN's MKL check in 2023 -- pure R, so it isolates the BLAS from any compiled code.
# Calibration: 1.0.6 must FAIL tests under mkl2023-avx and PASS under mkl (MKL 2026.1, AVX-512), as measured in
# ebmcd (dev/mkl/note-mkl.md). If it does not, the bench does not see what CRAN saw: read nothing else.
# Heavy enough to go under the machine lock (gate-lock-r): ~3 min per check, 9 checks.
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
. "$HERE/configs.sh"
REPO="$1"; OUT="$2"; REF="${3:-origin/master}"; CALIB="${4:-c413689}"
RBUILD="${RBUILD:-C:/Program Files/R/R-4.5.3/bin/R.exe}"
mkdir -p "$OUT"; rm -f "$OUT/summary.tsv"
build() {   # build <ref> <dir> [R CMD build options] -> prints the tarball path
	rm -rf "$OUT/$2"; mkdir -p "$OUT/$2"
	git -C "$REPO" archive "$1" | tar -x -C "$OUT/$2"
	(cd "$OUT" && "$RBUILD" CMD build "${@:3}" "$2" > "$2.build.log" 2>&1)
	ls -t "$OUT"/*.tar.gz | head -1
}
if [ "$CALIB" != none ]; then
	T106=$(build "$CALIB" src-calib --no-build-vignettes --no-manual)
	for c in mkl2023-avx mkl; do bash "$HERE/check_config.sh" "$c" "$T106" "$OUT/calib"; done
	cp "$OUT/calib/summary.tsv" "$OUT/calib-summary.tsv"
fi
T=$(build "$REF" src-check)
for c in $CONFIGS; do bash "$HERE/check_config.sh" "$c" "$T" "$OUT"; done
echo "RUN_WINDOWS_DONE"
