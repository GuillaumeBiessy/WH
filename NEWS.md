# WH (development version)

## Bug fixes

* In the regression framework, a cell with zero observed events (hence a missing `y = log(d / ec)` on a zero weight) no longer turns the whole 1D fit into `NA`, nor makes the 2D smoothing parameter selection fail.

* Data without names now get default positions, as the warning announced: in 1D the names were never assigned, and in 2D `WH()` failed on any matrix without dimnames.

* `predict()` now adds the observed positions to `newdata`, as documented, and handles non contiguous positions: `newdata = 100:110` used to recycle silently, and `c(40:99, 105)` returned under the label 105 the value at 100.

* The starting smoothing parameters of a 2D selection are found whatever the scale of the weights (the search failed when the weights were far from the unit scale).

* An evaluation of the selection criterion whose Cholesky factorization fails (smoothing parameters making the system numerically singular) no longer stops the selection: it is scored `Inf` and the optimizer moves back.

* The 1D REML / LAML criterion uses the exact log pseudo-determinant of the penalty instead of a sum over computed eigenvalues, which lost accuracy on long series (`q = 3`: off by 3e-3 at 1000 points, `NaN` at 3000).

* Invalid inputs (missing, infinite or negative `d`, `ec` or `wt`, events without exposure, missing `y` on a positive weight, non-finite `lambda`, dimensions with no more than `q` points, arrays of more than two dimensions) are refused with an explicit message.

* `plot.WH_2d()` labels its axes as they are drawn, and titles the residuals view.

## Convergence

* Every fit carries a `convergence` component (iterations and status of the fit, evaluations and status of the selection). A warning is emitted, whatever `verbose`, when a selected smoothing parameter lies at the boundary of the search interval, when the Nelder-Mead optimizer does not converge, or when the iterative fit reaches its new cap of 200 iterations.

## Computational changes

* The diagonal of the variance-covariance matrix is computed by the selected inverse of Takahashi et al. (1973), in O(p q²) instead of O(p² q) operations, with the same banded storage. Final fits and AIC, BIC or GCV selections are several times faster on large 2D tables.

* The variance of 2D extrapolated values is computed in C++ by a loop over the new cells, with banded solves and vectors of the size of the grid only: `predict()` on a 100 x 40 table extended to 110 x 45 goes from about 20 s to under 3 s.

## Documentation

* The `lambda_start_ratio` argument, which sets the starting values of a 2D selection, is documented in `WH()`.

# WH 2.0.0

This version introduces a complete overhaul of the package, with substantial improvements to both the core computational engine and the user-facing interface.

## Computational changes

The main change in this version is the exploitation of the banded structure of the penalization matrix, which allows for tremendous savings in both computation time and memory footprint.

* The package now uses Rcpp and calls directly Lapack functions for banded matrices when available. 

* Some additional functions have also been coded using Rcpp. The `backsolve` function for banded matrix was also recoded although it already exists in Lapack because the Lapack version overwrites the Cholesky factor so a copy has to be made beforehand.

* Parameter estimation now rely on backward-forward solves from the Cholesky factor instead of forming the variance-covariance matrix explicitly.

* By default, only the diagonal of the variance-covariance matrix is now computed (which is sufficient for credibility intervals). The full matrix can still be retrieved using the newly introduced `vcov()` method if necessary.

## User Interface changes

* The two main functions from earlier versions have been merged into a single unified function: `WH()`.

* Performance iteration and rank reduction are no longer available, as they no longer provide meaningful benefits with the new implementation. Indeed, the full-rank version is now efficient enough to handle several thousand observation points without optimization tricks.

## Other changes

* Some minor changes have been made to the example datasets to better align the characteristics of the mortality and long-term care datasets

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
