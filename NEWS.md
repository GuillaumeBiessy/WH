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
