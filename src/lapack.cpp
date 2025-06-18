#include <Rcpp.h>
#include <R_ext/Lapack.h>

using namespace Rcpp;

// Those functions call Lapack routines for banded matrices in compact storage

// [[Rcpp::export]]
NumericMatrix cholesky_compact_lapack(NumericMatrix& C, bool inplace = false) {

  NumericMatrix R = inplace ? C : clone(C);

  int n = R.ncol();
  int q = R.nrow() - 1;

  char uplo = 'U';
  int ldab = q + 1;
  int info;

  F77_CALL(dpbtrf)(&uplo, &n, &q, R.begin(), &ldab, &info FCONE);

  if (info != 0) stop("LAPACK dpbtrf failed with info = %d", info);

  return R;
}

// [[Rcpp::export]]
NumericMatrix invert_cholesky_compact_lapack(const NumericMatrix& C, bool transpose = false) {

  int n = C.ncol();
  int q = C.nrow() - 1;

  char uplo  = 'U';
  char trans = transpose ? 'T' : 'N';
  char diag  = 'N';
  int nrhs = n;
  int ldab = q + 1;
  int ldb  = n;
  int info;

  NumericMatrix R = clone(C);
  NumericMatrix X(n, n);
  for (int i = 0; i < n; ++i) X(i, i) = 1.0;

  F77_CALL(dtbtrs)(&uplo, &trans, &diag, &n, &q, &nrhs, R.begin(), &ldab, X.begin(), &ldb, &info FCONE FCONE FCONE);

  if (info != 0) stop("LAPACK dtbtrs failed with info = %d", info);

  return X;
}

// [[Rcpp::export]]
NumericVector backsolve_mat_compact_lapack(const NumericMatrix& C, const NumericMatrix& B, bool transpose = false) {

  int n = C.ncol();
  int q = C.nrow() - 1;

  char uplo  = 'U';
  char trans = transpose ? 'T' : 'N';
  char diag  = 'N';
  int nrhs = B.ncol();
  int ldab = q + 1;
  int ldb  = n;
  int info;

  NumericMatrix R = clone(C);
  NumericMatrix X = clone(B);

  F77_CALL(dtbtrs)(&uplo, &trans, &diag, &n, &q, &nrhs, R.begin(), &ldab, X.begin(), &ldb, &info FCONE FCONE FCONE);

  if (info != 0) stop("LAPACK dtbtrs failed with info = %d", info);

  return X;
}

// [[Rcpp::export]]
NumericVector eigenvalues_compact_lapack(const NumericMatrix& C) {

  int n   = C.ncol();          // Taille de la matrice
  int kd  = C.nrow() - 1;      // Nombre de bandes
  int ldab = kd + 1;
  int ldz  = 1;                // Pas besoin de vecteurs
  int info;

  NumericMatrix ab = clone(C);      // Compact bande, LAPACK-style
  NumericVector w(n);               // Valeurs propres
  NumericVector work(3 * n);        // Espace de travail

  char jobz = 'N';  // 'N' = only eigenvalues
  char uplo = 'U';  // bande supérieure

  F77_CALL(dsbev)(&jobz, &uplo, &n, &kd, ab.begin(), &ldab, w.begin(), nullptr, &ldz, work.begin(), &info FCONE FCONE);

  if (info != 0) stop("dsbev failed with info = %d", info);

  return w;
}

