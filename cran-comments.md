## R CMD check results

0 errors | 0 warnings | 0 notes

## Resubmission

This is a resubmission following recent package acceptance but failure to the MKL additional test.

In this version I have:

* Removed tests of the form `expect_equal(f(x), f(x))` added for testing purposes in the 1.0.6 version of the package
