#!/usr/bin/env bash
# check_config.sh -- R CMD check --as-cran of one tarball under the local substitute of one CRAN configuration.
#   bash check_config.sh <config> <tarball> <outdir>
# <config> is one of $CONFIGS (configs.sh). Writes <outdir>/<config>/ (check directory, check.log, probe.log) and
# appends one line to <outdir>/summary.tsv:
#   config  tarball  check-status  tests-line  blas-banner
# The environment of the configuration is exported in a subshell: it reaches every subprocess of the check
# (installation, examples, tests, vignettes), which all run the private copy's Rterm and load its DLLs.
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
. "$HERE/configs.sh"
CONF="$1"; TARBALL="$2"; OUT="$3"
D="$OUT/$CONF"
rm -rf "$D"; mkdir -p "$D"
# Link the package against the NAMES Rblas.dll / Rlapack.dll of the reference copy. Linked against the copy itself, the
# Rtools linker takes LAPACK symbols from the BLAS DLL when it also exports LAPACK (OpenBLAS), and records the DLL's
# internal name (libopenblas.dll, mkl_rt.3.dll): Windows then resolves that name on the PATH -- measured 2026-10-25,
# WH loaded Inkscape's libopenblas.dll, 40 times slower and not the OpenBLAS under test. Linked by name, the package
# binds at run time to the Rblas.dll / Rlapack.dll the private copy has already loaded: CRAN's -llapack -lRblas.
MV="$D/Makevars.link"
printf 'BLAS_LIBS = -L"%s/bin/x64" -lRblas\nLAPACK_LIBS = -L"%s/bin/x64" -lRlapack\n' \
	"$PRIVE/R-4.5.3-ref" "$PRIVE/R-4.5.3-ref" > "$MV"
(
	"config_$CONF" || exit 9
	export _R_CHECK_CRAN_INCOMING_REMOTE_=false R_MAKEVARS_USER="$MV"
	R="$RHOME/bin/R.exe"; RS="$RHOME/bin/Rscript.exe"
	MKL_VERBOSE=1 OPENBLAS_VERBOSE=2 "$RS" "$HERE/blas_probe.R" > "$D/probe.log" 2>&1
	"$R" CMD check --as-cran --output="$D" "$TARBALL" > "$D/check.log" 2>&1
	echo "rc=$?" >> "$D/check.log"
)
LOG=$(ls "$D"/*.Rcheck/00check.log 2>/dev/null | head -1)
STATUS=$(grep -h '^Status:' "$LOG" 2>/dev/null | tail -1 | sed 's/^Status: //')
TESTS=$(grep -h '\[ FAIL' "$D"/*.Rcheck/tests/testthat.Rout* 2>/dev/null | tail -1 | tr -s ' ')
MKLLINE=$(grep -m1 'MKL_VERBOSE oneMKL' "$D/probe.log" | tr -d '\r')
if [ -n "$MKLLINE" ]; then   # e.g. "oneMKL 2026 Update 1 AVX2 sequential"
	BANNER="$(echo "$MKLLINE" | grep -oE 'oneMKL [0-9.]+ Update [0-9]+') $(echo "$MKLLINE" \
		| grep -oE '\(Intel\(R\) [A-Z0-9-]+\)' | head -1 | sed 's/(Intel(R) //; s/)//') $(echo "$MKLLINE" | awk '{print $NF}')"
else
	BANNER=$(grep -m1 -oE 'Core: [A-Za-z0-9]+' "$D/probe.log" || echo "reference")
fi
# Which BLAS / LAPACK DLL names the compiled package imports: only Rblas.dll / Rlapack.dll are acceptable.
OD="${OBJDUMP:-C:/rtools45/x86_64-w64-mingw32.static.posix/bin/objdump.exe}"
IMPORTS=$(for f in "$D"/*.Rcheck/*/libs/x64/*.dll; do [ -f "$f" ] && "$OD" -p "$f" | grep 'DLL Name' \
	| grep -viE 'api-ms|kernel32|msvcrt|R\.dll|Rcpp|libgcc|libstdc'; done | awk '{print $3}' | sort -u | tr '\n' ' ')
printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$CONF" "$(basename "$TARBALL")" "${STATUS:-NO-STATUS}" "${TESTS:-no-testthat-line}" \
	"$BANNER" "imports: ${IMPORTS:-none}" >> "$OUT/summary.tsv"
tail -1 "$OUT/summary.tsv"
