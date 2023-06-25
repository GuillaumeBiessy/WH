## R CMD check results

0 errors | 0 warnings | 0 notes

## Resubmission

This is a resubmission following recent package acceptance but failure to the MKL additional test.

In this version I have:

* Increased test tolerance to avoid numerical test failure.

* Added sanity check test of the form `expect_equal(f(x), f(x))`. Failure to those tests would almost certainly indicate randomness in the output of MKL.
