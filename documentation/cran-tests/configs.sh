# configs.sh -- the local substitutes of CRAN's check configurations, sourced by check_config.sh.
# One function per configuration: it sets RHOME (the private R copy) and exports the environment of the check.
# Sources of each CRAN configuration (read on 2026-10-25):
#   check flavors   https://cran.r-project.org/web/checks/check_flavors.html
#   BLAS checks     https://www.stats.ox.ac.uk/pub/bdr/Rblas/README.txt
#   noSuggests      https://www.stats.ox.ac.uk/pub/bdr/noSuggests/README.txt
# Private R copies (never under C:/Program Files nor %LOCALAPPDATA%), built by the ebmcd repo:
#   dev/mkl/mkl_prive.sh (MKL copies), dev/blas-lapack/install_prive.sh (OpenBLAS); R-4.5.3-openblas-cran is
#   R-4.5.3-ref with Rblas.dll ONLY replaced by OpenBLAS 0.3.29 (MD5 a080acb0...), LAPACK kept as R's reference,
#   as CRAN does (libRblas.so symlinked to OpenBLAS, system LAPACK kept).
PRIVE="${PRIVE:-C:/Users/Guillaume Biessy/R-prive}"

CONFIGS="ref nosuggests openblas mkl mkl-avx2 mkl2023-avx2 mkl2023-avx"

config_ref() {            # r-release-windows-x86_64: R's reference BLAS and LAPACK
	RHOME="$PRIVE/R-4.5.3-ref"
}
config_nosuggests() {     # noSuggests: packages in Suggests unavailable
	RHOME="$PRIVE/R-4.5.3-ref"
	export _R_CHECK_DEPENDS_ONLY_=true
}
config_openblas() {       # OpenBLAS: serial OpenBLAS 0.3.29 as BLAS, reference LAPACK; CRAN CPU = Broadwell
	RHOME="$PRIVE/R-4.5.3-openblas-cran"
	export OPENBLAS_NUM_THREADS=1 OPENBLAS_CORETYPE=Haswell
}
config_mkl() {            # MKL 2026.1 serial, the kernels of this machine (AVX-512)
	RHOME="$PRIVE/R-4.5.3-mkl"
	export MKL_THREADING_LAYER=SEQUENTIAL
}
config_mkl-avx2() {       # MKL 2026.1 serial on CRAN's processor today (2 x Xeon E5-2690 v4, AVX2)
	RHOME="$PRIVE/R-4.5.3-mkl"
	export MKL_THREADING_LAYER=SEQUENTIAL MKL_ENABLE_INSTRUCTIONS=AVX2
}
config_mkl2023-avx2() {   # MKL 2023.2 serial (CRAN until 2026-07) on an AVX2 processor
	RHOME="$PRIVE/R-4.5.3-mkl2023"
	export MKL_THREADING_LAYER=SEQUENTIAL MKL_ENABLE_INSTRUCTIONS=AVX2
}
config_mkl2023-avx() {    # MKL 2023.2 serial on CRAN's processor of 2023-2024 (Xeon E5-2440, AVX): WH 1.x failed here
	RHOME="$PRIVE/R-4.5.3-mkl2023"
	export MKL_THREADING_LAYER=SEQUENTIAL MKL_ENABLE_INSTRUCTIONS=AVX
}
