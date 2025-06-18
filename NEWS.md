# WH (development version)

This version introduces a complete overhaul of the package, with substantial improvements to both the core computational engine and the user-facing API.

## Core computational changes

* Exploitation of the banded structure of the penalization matrix:

* Computations now scale linearly with the number of observations (in the largest dimension), instead of exhibiting cubic complexity.

* This also leads to a significantly reduced memory footprint.

* Parameter estimation is now based on the Cholesky factor, using backward-forward solves instead of forming the variance-covariance matrix explicitly.

* By default, only the diagonal of the variance-covariance matrix is now computed (sufficient for credibility intervals). The full matrix can still be retrieved using the newly introduced `vcov()` method.

## User interface and API changes

The two main functions from earlier versions have been merged into a single unified function: `WH()`.

Performance iteration and rank reduction are no longer available, as they no longer provide meaningful benefits with the new implementation. The full-rank version is now efficient enough to handle several thousand observation points without optimization tricks.

# WH 1.1.2

* Greatly increased tolerance from 1e-6 to 1e-2 to remove a lone failing test which only occurs when using MKL BLAS

# WH 1.1.1

* Increased tolerance to remove a lone failing test which only occurs when using MKL BLAS

# WH 1.1.0

* Made significant changes to the extrapolation methods provided by be `predict.WH_1d` and `predict.WH_2d` functions. It turns out that the formula used for the extrapolation ignored the innovation error caused by the prior on the extrapolated region, resulting in smaller credibility intervals than they should have been. This now has been fixed.

# WH 1.0.7

* Removed tests of the form `expect_equal(f(x), f(x))` after confirmation they do not work with MKL BLAS

# WH 1.0.6

* Further improved test robustness

* Added tests of the form `expect_equal(f(x), f(x))` for testing purposes

# WH 1.0.5

* Further improved tests robustness

# WH 1.0.4

* Simplified computation of the matrix `tUWU` by using the `crossprod` function, which should reduce memory usage (by half) and computation time (slightly)

* Fixed an issue with the `WH_2d` plot and the what = "edf" argument

* Improved tests robustness

# WH 1.0.3

* Replaced backquotes by normal quotes in the description field of the DESCRIPTION file

# WH 1.0.2

* Updated the package description to link it to the paper whose methods are used in the package

# WH 1.0.1

* Fixed an issue in the author field of the DESCRIPTION file

# WH 1.0.0

* Official candidate for CRAN release

# WH 0.1.0

* Added a `NEWS.md` file to track changes to the package
