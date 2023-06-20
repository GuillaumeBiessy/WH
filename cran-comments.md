## R CMD check results

0 errors | 0 warnings | 0 notes

## Resubmission

This is a resubmission following recent package acceptance but failure to additional tests using the alternative BLAS / LAPACK implementation from Intel MKL '2022.3.0'. I suspect this is linked to some variability in the results of the eigendecomposition when performed by MKL.

In this version I have:

* Improved the robustness of the model tests by only considering the model fit, esimated standard deviation and diagnoses while comparing two fits rather than comparing the whole object. This should hopefully prevent test failure when using the MKL implementation of BLAS / LAPACK.

* Improved the presentation of unit tests by embedding them in calls to `test_that`
