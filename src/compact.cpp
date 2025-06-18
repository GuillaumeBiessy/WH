#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix create_P_compact_cpp(int p, int q) {

  std::vector<double> b(q + 1);
  b[0] = 1.0;
  for (int k = 1; k <= q; ++k)
    b[k] = - b[k - 1] * (q + 1 - k) / k;

  NumericMatrix C(q + 1, p);
  for (int i = 0; i < p - q; ++i) {
    int k_max = std::min(q, p - 1 - i);
    for (int k = 0; k <= k_max; ++k) {
      double bk = b[k];
      for (int l = k; l <= k_max; ++l) {
        C(q + k - l, i + l) += bk * b[l];
      }
    }
  }
  return C;
}

// [[Rcpp::export]]
NumericVector backsolve_compact_cpp(const NumericMatrix& C, const NumericVector& b, bool transpose = false) {

  int p = C.ncol();
  int q = C.nrow() - 1;
  NumericVector x(p);

  if (transpose) {

    for (int i = 0; i < p; ++i) {
      double sum = 0.0;
      int k = i - q;
      int j_min = std::max(k, 0);
      for (int j = j_min; j < i; ++j) {
        sum += C(j - k, i) * x[j];
      }
      x[i] = (b[i] - sum) / C(q, i);
    }

  } else {

    for (int i = p - 1; i >= 0; --i) {
      double sum = 0.0;
      int k = i + q;
      int j_max = std::min(k, p - 1);
      for (int j = i + 1; j <= j_max; ++j) {
        sum += C(k - j, j) * x[j];
      }
      x[i] = (b[i] - sum) / C(q, i);
    }
  }
  return x;
}

// [[Rcpp::export]]
NumericVector diag_V_compact_cpp(const NumericMatrix& C) {

  int p = C.ncol();
  int q = C.nrow() - 1;

  NumericVector result(p);
  std::vector<double> K_col(p);
  std::vector<double> K_diag(p);

  for (int i = 0; i < p; ++i) {
    double Kii = 1.0 / C(q, i);
    K_diag[i] = Kii;
    result[i] = Kii * Kii;
  }

  for (int j = p - 1; j >= 0; --j) {

    K_col[j] = K_diag[j];

    for (int i = j - 1; i >= 0; --i) {
      int delta = i + q;
      double sum = 0.0;

      int k_max = std::min(j, delta);
      for (int k = i + 1; k <= k_max; ++k) {
        sum += C(delta - k, k) * K_col[k];
      }
      double Kij = - K_diag[i] * sum;
      K_col[i] = Kij;
      result[i] += Kij * Kij;
    }
  }
  return result;
}

// [[Rcpp::export]]
NumericVector get_prod_P_y_compact_cpp(const NumericMatrix& y, const List& C_list) {

  const NumericMatrix& C1 = C_list[0];
  const NumericMatrix& C2 = C_list[1];

  const int p = y.nrow();
  const int q1 = C1.nrow() - 1;
  const int q2 = C2.nrow() - 1;
  const int n_blocs = y.ncol();

  NumericVector z(p * n_blocs, 0.0);

  // Bloc 1 : ligne i fixe
  for (int i = 0; i < p; ++i) {
    const double Pii = C1(q1, i);
    for (int j = 0; j < n_blocs; ++j)
      z[j * p + i] += Pii * y(i, j);

    const int d_left  = std::min(q1, p - 1 - i);
    const int d_right = std::min(q1, i);
    const int d_sym   = std::min(d_left, d_right);

    for (int d = 1; d <= d_sym; ++d) {
      const int k_inf = i - d;
      const int k_sup = i + d;
      const double Pki = C1(q1 - d, i);
      const double Pik = C1(q1 - d, k_sup);
      for (int j = 0; j < n_blocs; ++j)
        z[j * p + i] += Pki * y(k_inf, j) + Pik * y(k_sup, j);
    }

    for (int d = d_sym + 1; d <= d_left; ++d) {
      const int k_sup = i + d;
      const double Pik = C1(q1 - d, k_sup);
      for (int j = 0; j < n_blocs; ++j)
        z[j * p + i] += Pik * y(k_sup, j);
    }

    for (int d = d_sym + 1; d <= d_right; ++d) {
      const int k_inf = i - d;
      const double Pki = C1(q1 - d, i);
      for (int j = 0; j < n_blocs; ++j)
        z[j * p + i] += Pki * y(k_inf, j);
    }
  }

  // Bloc 2 : colonne j fixe
  for (int j = 0; j < n_blocs; ++j) {
    const int offset = j * p;
    const double Pjj = C2(q2, j);

    for (int i = 0; i < p; ++i)
      z[offset + i] += Pjj * y(i, j);

    const int d_left  = std::min(q2, n_blocs - 1 - j);
    const int d_right = std::min(q2, j);
    const int d_sym   = std::min(d_left, d_right);

    for (int d = 1; d <= d_sym; ++d) {
      const int k_inf = j - d;
      const int k_sup = j + d;
      const double Pkj = C2(q2 - d, j);
      const double Pjk = C2(q2 - d, k_sup);
      for (int i = 0; i < p; ++i)
        z[offset + i] += Pkj * y(i, k_inf) + Pjk * y(i, k_sup);
    }

    for (int d = d_sym + 1; d <= d_left; ++d) {
      const int k_sup = j + d;
      const double Pjk = C2(q2 - d, k_sup);
      for (int i = 0; i < p; ++i)
        z[offset + i] += Pjk * y(i, k_sup);
    }

    for (int d = d_sym + 1; d <= d_right; ++d) {
      const int k_inf = j - d;
      const double Pkj = C2(q2 - d, j);
      for (int i = 0; i < p; ++i)
        z[offset + i] += Pkj * y(i, k_inf);
    }
  }
  return z;
}

// [[Rcpp::export]]
NumericMatrix get_prod_P_K_compact_cpp(const NumericMatrix& K, const List& C_list) {

  const NumericMatrix& C1 = C_list[0];  // bloc 1
  const NumericMatrix& C2 = C_list[1];  // bloc 2

  const int n = K.nrow();

  // ----- Bloc 1 -----
  const int q1 = C1.nrow() - 1;
  const int p1 = C1.ncol();
  const int n_blocs1 = n / p1;

  // ----- Bloc 2 -----
  const int q2 = C2.nrow() - 1;
  const int n_blocs2 = C2.ncol();
  const int p2 = n / n_blocs2;

  NumericMatrix R(n, n);

  // ---------- Bloc 1 ----------
  for (int i_bloc = 0; i_bloc < n_blocs1; ++i_bloc) {
    int offset_i = i_bloc * p1;

    // Diagonal block
    {
      int offset_j = offset_i;

      for (int i = 0; i < p1; ++i) {
        int global_i = offset_i + i;
        double Pii = C1(q1, i);

        for (int j = 0; j < p1; ++j) {
          int global_j = offset_j + j;
          if (global_i <= global_j)
            R(global_i, global_j) += Pii * K(global_i, global_j);
        }

        int d_left  = std::min(q1, i);
        int d_right = std::min(q1, p1 - 1 - i);
        int d_sym   = std::min(d_left, d_right);

        for (int d = 1; d <= d_sym; ++d) {
          int gi_km = offset_i + (i - d);
          int gi_kp = offset_i + (i + d);
          double Pki = C1(q1 - d, i);
          double Pik = C1(q1 - d, i + d);

          for (int j = 0; j < p1; ++j) {
            int global_j = offset_j + j;
            if (gi_km <= global_j)
              R(global_i, global_j) += Pki * K(gi_km, global_j);
            if (gi_kp <= global_j)
              R(global_i, global_j) += Pik * K(gi_kp, global_j);
          }
        }

        for (int d = d_sym + 1; d <= d_left; ++d) {
          int gi_km = offset_i + (i - d);
          double Pki = C1(q1 - d, i);

          for (int j = 0; j < p1; ++j) {
            int global_j = offset_j + j;
            if (gi_km <= global_j)
              R(global_i, global_j) += Pki * K(gi_km, global_j);
          }
        }

        for (int d = d_sym + 1; d <= d_right; ++d) {
          int gi_kp = offset_i + (i + d);
          double Pik = C1(q1 - d, i + d);

          for (int j = 0; j < p1; ++j) {
            int global_j = offset_j + j;
            if (gi_kp <= global_j)
              R(global_i, global_j) += Pik * K(gi_kp, global_j);
          }
        }
      }
    }

    // Off-diagonal blocks
    for (int j_bloc = i_bloc + 1; j_bloc < n_blocs1; ++j_bloc) {
      int offset_j = j_bloc * p1;

      for (int i = 0; i < p1; ++i) {
        int global_i = offset_i + i;
        double Pii = C1(q1, i);

        for (int j = 0; j < p1; ++j) {
          int global_j = offset_j + j;
          R(global_i, global_j) += Pii * K(global_i, global_j);
        }

        int d_left  = std::min(q1, i);
        int d_right = std::min(q1, p1 - 1 - i);
        int d_sym   = std::min(d_left, d_right);

        for (int d = 1; d <= d_sym; ++d) {
          int gi_km = offset_i + (i - d);
          int gi_kp = offset_i + (i + d);
          double Pki = C1(q1 - d, i);
          double Pik = C1(q1 - d, i + d);

          for (int j = 0; j < p1; ++j) {
            int global_j = offset_j + j;
            R(global_i, global_j) += Pki * K(gi_km, global_j) + Pik * K(gi_kp, global_j);
          }
        }

        for (int d = d_sym + 1; d <= d_left; ++d) {
          int gi_km = offset_i + (i - d);
          double Pki = C1(q1 - d, i);

          for (int j = 0; j < p1; ++j) {
            int global_j = offset_j + j;
            R(global_i, global_j) += Pki * K(gi_km, global_j);
          }
        }

        for (int d = d_sym + 1; d <= d_right; ++d) {
          int gi_kp = offset_i + (i + d);
          double Pik = C1(q1 - d, i + d);

          for (int j = 0; j < p1; ++j) {
            int global_j = offset_j + j;
            R(global_i, global_j) += Pik * K(gi_kp, global_j);
          }
        }
      }
    }
  }

  // ---------- Bloc 2 ----------
  for (int i = 0; i < n_blocs2; ++i) {
    int offset_i = i * p2;

    int d_left  = std::min(q2, i);
    int d_right = std::min(q2, n_blocs2 - 1 - i);
    int d_sym   = std::min(d_left, d_right);

    double coeff_diag = C2(q2, i);
    for (int row = 0; row < p2; ++row) {
      int row_i = offset_i + row;
      for (int col = row_i; col < n; ++col) {
        R(row_i, col) += coeff_diag * K(row_i, col);
      }
    }

    for (int d = 1; d <= d_sym; ++d) {
      int k1 = i - d;
      int k2 = i + d;
      int offset_k1 = k1 * p2;
      int offset_k2 = k2 * p2;
      double coeff_k1 = C2(q2 - d, i);
      double coeff_k2 = C2(q2 - d, k2);

      for (int row = 0; row < p2; ++row) {
        int row_i = offset_i + row;
        int row_k1 = offset_k1 + row;
        int row_k2 = offset_k2 + row;

        for (int col = row_k1; col < n; ++col)
          R(row_i, col) += coeff_k1 * K(row_k1, col);

        for (int col = row_k2; col < n; ++col)
          R(row_i, col) += coeff_k2 * K(row_k2, col);
      }
    }

    for (int d = d_sym + 1; d <= d_left; ++d) {
      int k1 = i - d;
      int offset_k1 = k1 * p2;
      double coeff = C2(q2 - d, i);

      for (int row = 0; row < p2; ++row) {
        int row_i = offset_i + row;
        int row_k1 = offset_k1 + row;

        for (int col = row_k1; col < n; ++col)
          R(row_i, col) += coeff * K(row_k1, col);
      }
    }

    for (int d = d_sym + 1; d <= d_right; ++d) {
      int k2 = i + d;
      int offset_k2 = k2 * p2;
      double coeff = C2(q2 - d, k2);

      for (int row = 0; row < p2; ++row) {
        int row_i = offset_i + row;
        int row_k2 = offset_k2 + row;

        for (int col = row_k2; col < n; ++col)
          R(row_i, col) += coeff * K(row_k2, col);
      }
    }
  }
  return R;
}

// [[Rcpp::export]]
NumericMatrix submatrix_compact_cpp(const NumericMatrix& C, const IntegerVector& drop_indices) {

  int q = C.nrow() - 1;
  int n = C.ncol();

  LogicalVector drop_flags(n, false);
  for (int idx : drop_indices) {
    drop_flags[idx - 1] = true;  // 1-based → 0-based
  }

  std::vector<int> keep_cols;
  for (int j = 0; j < n; ++j) {
    if (!drop_flags[j]) keep_cols.push_back(j);
  }

  int n_new = keep_cols.size();
  NumericMatrix C_new(q + 1, n_new);

  // Pré-calcul de i_shift
  std::vector<int> i_shift(n + 1, 0);
  int pos = 0;
  for (int i = 0; i <= n; ++i) {
    while (pos < drop_indices.size() && drop_indices[pos] <= i)
      ++pos;
    i_shift[i] = pos;
  }

  for (int jj = 0; jj < n_new; ++jj) {
    int j = keep_cols[jj];
    for (int r = 0; r <= q; ++r) {

      int i = j - (q - r);
      if (i < 0 || drop_flags[i]) continue;

      int i_new = i + 1 - i_shift[i + 1];
      int r_new = q - (jj + 1 - i_new);

      if (r_new >= 0 && r_new <= q) {
        C_new(r_new, jj) = C(r, j);
      }
    }
  }
  return C_new;
}


