
<!-- README.md is generated from README.Rmd. Please edit that file -->

# A Modern Take on Whittaker-Henderson Smoothing

<!-- badges: start -->

[![R-CMD-check](https://github.com/GuillaumeBiessy/WH/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/GuillaumeBiessy/WH/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

## What is Whittaker-Henderson smoothing ?

### Origin

The Whittaker-Henderson (WH) smoothing is a graduation method which
attenuates the impact of sample fluctuations. Initially introduced by
Whittaker (1922) for the construction of mortality tables and improved
by the work of Henderson (1924), it remains to date one of the most
popular smoothing method among actuaries working on the modelling of
person insurance risks such as death, disability, long-term care and
unemployement. Whittaker-Henderson smoothing generalizes to
two-dimension smoothing but requires in all cases that the observation
be equally spaced on each dimension.

### The one-dimension case

Let
![y](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;y "y")
be an observation vector and
![w \ge 0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;w%20%5Cge%200 "w \ge 0")
a vector of positive weights, both of size
![n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n "n").
The estimator associated with WH smoothing is:

![\hat{y} = \underset{\theta}{\text{argmin}}\\{F(y,w,\theta) + R\_{\lambda,q}(\theta)\\}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Chat%7By%7D%20%3D%20%5Cunderset%7B%5Ctheta%7D%7B%5Ctext%7Bargmin%7D%7D%5C%7BF%28y%2Cw%2C%5Ctheta%29%20%2B%20R_%7B%5Clambda%2Cq%7D%28%5Ctheta%29%5C%7D "\hat{y} = \underset{\theta}{\text{argmin}}\{F(y,w,\theta) + R_{\lambda,q}(\theta)\}")

where:

- ![F(y,w,\theta) = \underset{i = 1}{\overset{n}{\sum}} w_i(y_i - \theta_i)^2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;F%28y%2Cw%2C%5Ctheta%29%20%3D%20%5Cunderset%7Bi%20%3D%201%7D%7B%5Coverset%7Bn%7D%7B%5Csum%7D%7D%20w_i%28y_i%20-%20%5Ctheta_i%29%5E2 "F(y,w,\theta) = \underset{i = 1}{\overset{n}{\sum}} w_i(y_i - \theta_i)^2")
  is a fidelity criterion and

- ![R(\theta,\lambda,q) = \lambda \underset{i = 1}{\overset{n - q}{\sum}} (\Delta^q\theta)\_i^2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;R%28%5Ctheta%2C%5Clambda%2Cq%29%20%3D%20%5Clambda%20%5Cunderset%7Bi%20%3D%201%7D%7B%5Coverset%7Bn%20-%20q%7D%7B%5Csum%7D%7D%20%28%5CDelta%5Eq%5Ctheta%29_i%5E2 "R(\theta,\lambda,q) = \lambda \underset{i = 1}{\overset{n - q}{\sum}} (\Delta^q\theta)_i^2")
  a regularity criterion

![\Delta^q](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5CDelta%5Eq "\Delta^q")
represents in the above expression the difference operator of order
![q](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;q "q")
such that for all
![i\in\[1,n - q\]](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;i%5Cin%5B1%2Cn%20-%20q%5D "i\in[1,n - q]"):

![(\Delta^q\theta)\_i = \underset{k = 0}{\overset{q}{\sum}} \begin{pmatrix}q \\\\ k\end{pmatrix} (- 1)^{q - k} \theta\_{i + k}.](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%28%5CDelta%5Eq%5Ctheta%29_i%20%3D%20%5Cunderset%7Bk%20%3D%200%7D%7B%5Coverset%7Bq%7D%7B%5Csum%7D%7D%20%5Cbegin%7Bpmatrix%7Dq%20%5C%5C%20k%5Cend%7Bpmatrix%7D%20%28-%201%29%5E%7Bq%20-%20k%7D%20%5Ctheta_%7Bi%20%2B%20k%7D. "(\Delta^q\theta)_i = \underset{k = 0}{\overset{q}{\sum}} \begin{pmatrix}q \\ k\end{pmatrix} (- 1)^{q - k} \theta_{i + k}.")

Let us define
![W = \text{Diag}(w)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;W%20%3D%20%5Ctext%7BDiag%7D%28w%29 "W = \text{Diag}(w)")
the diagonal matrix of weights and
![D\_{n,q}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;D_%7Bn%2Cq%7D "D_{n,q}")
the matrix of differences of order
![q](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;q "q"),
of dimensions
![(n - q,n)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%28n%20-%20q%2Cn%29 "(n - q,n)"),
such that
![(D\_{n,q}\theta)\_i = (\Delta^q\theta)\_i](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%28D_%7Bn%2Cq%7D%5Ctheta%29_i%20%3D%20%28%5CDelta%5Eq%5Ctheta%29_i "(D_{n,q}\theta)_i = (\Delta^q\theta)_i").
Actually, only differences of order
![q\in\\{1,2\\}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;q%5Cin%5C%7B1%2C2%5C%7D "q\in\{1,2\}")
ar of practical interest. The associated difference matrices are
represented below:

![D_1 = \begin{pmatrix}
1 & - 1 &  & 0 \\\\
& \ddots & \ddots & \\\\
0 & & 1 & - 1
\end{pmatrix}
\quad\quad
D_2 = \begin{pmatrix}
1 & - 2 & 1 & & 0 \\\\
& \ddots & \ddots & \ddots & \\\\
0 & & 1 & - 2 & 1
\end{pmatrix}.](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;D_1%20%3D%20%5Cbegin%7Bpmatrix%7D%0A1%20%26%20-%201%20%26%20%20%26%200%20%5C%5C%0A%26%20%5Cddots%20%26%20%5Cddots%20%26%20%5C%5C%0A0%20%26%20%26%201%20%26%20-%201%0A%5Cend%7Bpmatrix%7D%0A%5Cquad%5Cquad%0AD_2%20%3D%20%5Cbegin%7Bpmatrix%7D%0A1%20%26%20-%202%20%26%201%20%26%20%26%200%20%5C%5C%0A%26%20%5Cddots%20%26%20%5Cddots%20%26%20%5Cddots%20%26%20%5C%5C%0A0%20%26%20%26%201%20%26%20-%202%20%26%201%0A%5Cend%7Bpmatrix%7D. "D_1 = \begin{pmatrix}
1 & - 1 &  & 0 \\
& \ddots & \ddots & \\
0 & & 1 & - 1
\end{pmatrix}
\quad\quad
D_2 = \begin{pmatrix}
1 & - 2 & 1 & & 0 \\
& \ddots & \ddots & \ddots & \\
0 & & 1 & - 2 & 1
\end{pmatrix}.")

The fidelity and regularity criteria may be rewritten using matrix
notations:

![\begin{aligned}
F(y,w,\theta) &= \underset{i = 1}{\overset{n}{\sum}} w_i(y_i - \theta_i)^2 = \Vert\sqrt{W}(y - \theta)\Vert^2 = (y - \theta)^TW(y - \theta) \\\\
R(\theta,\lambda,q) &= \lambda \underset{i = 1}{\overset{n - q}{\sum}} (\Delta^q\theta)\_i^2 = \lambda\Vert D\_{n,q}\theta\Vert^2 = \lambda\theta^TD\_{n,q}^TD\_{n,q}\theta
\end{aligned}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cbegin%7Baligned%7D%0AF%28y%2Cw%2C%5Ctheta%29%20%26%3D%20%5Cunderset%7Bi%20%3D%201%7D%7B%5Coverset%7Bn%7D%7B%5Csum%7D%7D%20w_i%28y_i%20-%20%5Ctheta_i%29%5E2%20%3D%20%5CVert%5Csqrt%7BW%7D%28y%20-%20%5Ctheta%29%5CVert%5E2%20%3D%20%28y%20-%20%5Ctheta%29%5ETW%28y%20-%20%5Ctheta%29%20%5C%5C%0AR%28%5Ctheta%2C%5Clambda%2Cq%29%20%26%3D%20%5Clambda%20%5Cunderset%7Bi%20%3D%201%7D%7B%5Coverset%7Bn%20-%20q%7D%7B%5Csum%7D%7D%20%28%5CDelta%5Eq%5Ctheta%29_i%5E2%20%3D%20%5Clambda%5CVert%20D_%7Bn%2Cq%7D%5Ctheta%5CVert%5E2%20%3D%20%5Clambda%5Ctheta%5ETD_%7Bn%2Cq%7D%5ETD_%7Bn%2Cq%7D%5Ctheta%0A%5Cend%7Baligned%7D "\begin{aligned}
F(y,w,\theta) &= \underset{i = 1}{\overset{n}{\sum}} w_i(y_i - \theta_i)^2 = \Vert\sqrt{W}(y - \theta)\Vert^2 = (y - \theta)^TW(y - \theta) \\
R(\theta,\lambda,q) &= \lambda \underset{i = 1}{\overset{n - q}{\sum}} (\Delta^q\theta)_i^2 = \lambda\Vert D_{n,q}\theta\Vert^2 = \lambda\theta^TD_{n,q}^TD_{n,q}\theta
\end{aligned}")

and the associated estimator becomes:

![\hat{y} = \underset{\theta}{\text{argmin}} \left\lbrace(y - \theta)^TW(y - \theta) + \theta^TP\_\lambda\theta\right\rbrace](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Chat%7By%7D%20%3D%20%5Cunderset%7B%5Ctheta%7D%7B%5Ctext%7Bargmin%7D%7D%20%5Cleft%5Clbrace%28y%20-%20%5Ctheta%29%5ETW%28y%20-%20%5Ctheta%29%20%2B%20%5Ctheta%5ETP_%5Clambda%5Ctheta%5Cright%5Crbrace "\hat{y} = \underset{\theta}{\text{argmin}} \left\lbrace(y - \theta)^TW(y - \theta) + \theta^TP_\lambda\theta\right\rbrace")

by noting in the one-dimension case
![P\_\lambda = \lambda D\_{n,q}^TD\_{n,q}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;P_%5Clambda%20%3D%20%5Clambda%20D_%7Bn%2Cq%7D%5ETD_%7Bn%2Cq%7D "P_\lambda = \lambda D_{n,q}^TD_{n,q}").

### The two-dimension case

In the two-dimension case, we start from an observation matrix
![Y](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;Y "Y")
and an associated matrix
![\Omega](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5COmega "\Omega")
of positive weights, both of dimensions
![n_x \times n_z](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n_x%20%5Ctimes%20n_z "n_x \times n_z").

The estimator associated with Whittaker-Henderson smoothing in this case
reads:

![\widehat{Y} = \underset{\Theta}{\text{argmin}}\\{F(Y,\Omega, \Theta) + R\_{\lambda,q}(\Theta)\\}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cwidehat%7BY%7D%20%3D%20%5Cunderset%7B%5CTheta%7D%7B%5Ctext%7Bargmin%7D%7D%5C%7BF%28Y%2C%5COmega%2C%20%5CTheta%29%20%2B%20R_%7B%5Clambda%2Cq%7D%28%5CTheta%29%5C%7D "\widehat{Y} = \underset{\Theta}{\text{argmin}}\{F(Y,\Omega, \Theta) + R_{\lambda,q}(\Theta)\}")

where:

![\begin{aligned}
F(Y,\Omega, \Theta) &= \underset{i = 1}{\overset{n_x}{\sum}}\underset{j = 1}{\overset{n_z}{\sum}} \Omega\_{i,j}(Y\_{i,j} - \Theta\_{i,j})^2\quad \text{is a fidelity criterion} \\\\
R(\Theta,\lambda,q) &= \lambda_x \underset{j = 1}{\overset{n_z}{\sum}}\underset{i = 1}{\overset{n_x - q_x}{\sum}} (\Delta^{q_x}\Theta\_{\bullet,j})\_i^2 + \lambda_z \underset{i = 1}{\overset{n_x}{\sum}}\underset{j = 1}{\overset{n_z - q_z}{\sum}} (\Delta^{q_z}\Theta\_{i,\bullet})\_j^2 \quad \text{is a regularity criterion.}
\end{aligned}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cbegin%7Baligned%7D%0AF%28Y%2C%5COmega%2C%20%5CTheta%29%20%26%3D%20%5Cunderset%7Bi%20%3D%201%7D%7B%5Coverset%7Bn_x%7D%7B%5Csum%7D%7D%5Cunderset%7Bj%20%3D%201%7D%7B%5Coverset%7Bn_z%7D%7B%5Csum%7D%7D%20%5COmega_%7Bi%2Cj%7D%28Y_%7Bi%2Cj%7D%20-%20%5CTheta_%7Bi%2Cj%7D%29%5E2%5Cquad%20%5Ctext%7Bis%20a%20fidelity%20criterion%7D%20%5C%5C%0AR%28%5CTheta%2C%5Clambda%2Cq%29%20%26%3D%20%5Clambda_x%20%5Cunderset%7Bj%20%3D%201%7D%7B%5Coverset%7Bn_z%7D%7B%5Csum%7D%7D%5Cunderset%7Bi%20%3D%201%7D%7B%5Coverset%7Bn_x%20-%20q_x%7D%7B%5Csum%7D%7D%20%28%5CDelta%5E%7Bq_x%7D%5CTheta_%7B%5Cbullet%2Cj%7D%29_i%5E2%20%2B%20%5Clambda_z%20%5Cunderset%7Bi%20%3D%201%7D%7B%5Coverset%7Bn_x%7D%7B%5Csum%7D%7D%5Cunderset%7Bj%20%3D%201%7D%7B%5Coverset%7Bn_z%20-%20q_z%7D%7B%5Csum%7D%7D%20%28%5CDelta%5E%7Bq_z%7D%5CTheta_%7Bi%2C%5Cbullet%7D%29_j%5E2%20%5Cquad%20%5Ctext%7Bis%20a%20regularity%20criterion.%7D%0A%5Cend%7Baligned%7D "\begin{aligned}
F(Y,\Omega, \Theta) &= \underset{i = 1}{\overset{n_x}{\sum}}\underset{j = 1}{\overset{n_z}{\sum}} \Omega_{i,j}(Y_{i,j} - \Theta_{i,j})^2\quad \text{is a fidelity criterion} \\
R(\Theta,\lambda,q) &= \lambda_x \underset{j = 1}{\overset{n_z}{\sum}}\underset{i = 1}{\overset{n_x - q_x}{\sum}} (\Delta^{q_x}\Theta_{\bullet,j})_i^2 + \lambda_z \underset{i = 1}{\overset{n_x}{\sum}}\underset{j = 1}{\overset{n_z - q_z}{\sum}} (\Delta^{q_z}\Theta_{i,\bullet})_j^2 \quad \text{is a regularity criterion.}
\end{aligned}")

The latter criterion may be decomposed as the sum of:

- a one-dimension regularity criterion applied to all rows of
  ![\Theta](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5CTheta "\Theta")
  and

- a one-dimension regularity criterion applied to all columns of
  ![\Theta](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5CTheta "\Theta").

Once again it is more convenient to adopt matrix notations by defining
![y = \text{vec}(Y)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;y%20%3D%20%5Ctext%7Bvec%7D%28Y%29 "y = \text{vec}(Y)"),
![w = \text{vec}(\Omega)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;w%20%3D%20%5Ctext%7Bvec%7D%28%5COmega%29 "w = \text{vec}(\Omega)"),
![\theta = \text{vec}(\Theta)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctheta%20%3D%20%5Ctext%7Bvec%7D%28%5CTheta%29 "\theta = \text{vec}(\Theta)")
the vectors obtained by concatenating the columns of
![Y](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;Y "Y"),
![\Omega](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5COmega "\Omega")
and
![\Theta](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5CTheta "\Theta")
respectively and by noting
![W = \text{Diag}(w)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;W%20%3D%20%5Ctext%7BDiag%7D%28w%29 "W = \text{Diag}(w)")
and
![n = n_x \times n_z](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n%20%3D%20n_x%20%5Ctimes%20n_z "n = n_x \times n_z").

The fidelity and regularity criteria may this be expressed as funuctions
of
![y](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;y "y"),
![w](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;w "w")
and
![\theta](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctheta "\theta"),
using matrix notations:

![\begin{aligned}
F(y,w, \theta) &= \underset{i = 1}{\overset{n}{\sum}}w_i(y_i - \theta_i)^2 = (y - \theta)^TW(y - \theta) \\\\
R(\theta,\lambda,q) &= \theta^{T}(\lambda_x I\_{n_z} \otimes D\_{n_x,q_x}^{T}D\_{n_x,q_x} + \lambda_z D\_{n_z,q_z}^{T}D\_{n_z,q_z} \otimes I\_{n_x}) \theta
\end{aligned}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cbegin%7Baligned%7D%0AF%28y%2Cw%2C%20%5Ctheta%29%20%26%3D%20%5Cunderset%7Bi%20%3D%201%7D%7B%5Coverset%7Bn%7D%7B%5Csum%7D%7Dw_i%28y_i%20-%20%5Ctheta_i%29%5E2%20%3D%20%28y%20-%20%5Ctheta%29%5ETW%28y%20-%20%5Ctheta%29%20%5C%5C%0AR%28%5Ctheta%2C%5Clambda%2Cq%29%20%26%3D%20%5Ctheta%5E%7BT%7D%28%5Clambda_x%20I_%7Bn_z%7D%20%5Cotimes%20D_%7Bn_x%2Cq_x%7D%5E%7BT%7DD_%7Bn_x%2Cq_x%7D%20%2B%20%5Clambda_z%20D_%7Bn_z%2Cq_z%7D%5E%7BT%7DD_%7Bn_z%2Cq_z%7D%20%5Cotimes%20I_%7Bn_x%7D%29%20%5Ctheta%0A%5Cend%7Baligned%7D "\begin{aligned}
F(y,w, \theta) &= \underset{i = 1}{\overset{n}{\sum}}w_i(y_i - \theta_i)^2 = (y - \theta)^TW(y - \theta) \\
R(\theta,\lambda,q) &= \theta^{T}(\lambda_x I_{n_z} \otimes D_{n_x,q_x}^{T}D_{n_x,q_x} + \lambda_z D_{n_z,q_z}^{T}D_{n_z,q_z} \otimes I_{n_x}) \theta
\end{aligned}")

which also leads to the equation:

![\hat{y} = \underset{\theta}{\text{argmin}} \left\lbrace(y - \theta)^TW(y - \theta) + \theta^TP\_\lambda\theta\right\rbrace](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Chat%7By%7D%20%3D%20%5Cunderset%7B%5Ctheta%7D%7B%5Ctext%7Bargmin%7D%7D%20%5Cleft%5Clbrace%28y%20-%20%5Ctheta%29%5ETW%28y%20-%20%5Ctheta%29%20%2B%20%5Ctheta%5ETP_%5Clambda%5Ctheta%5Cright%5Crbrace "\hat{y} = \underset{\theta}{\text{argmin}} \left\lbrace(y - \theta)^TW(y - \theta) + \theta^TP_\lambda\theta\right\rbrace")

except in this case:

![P\_\lambda = \lambda_x I\_{n_z} \otimes D\_{n_x,q_x}^{T}D\_{n_x,q_x} + \lambda_z D\_{n_z,q_z}^{T}D\_{n_z,q_z} \otimes I\_{n_x}.](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;P_%5Clambda%20%3D%20%5Clambda_x%20I_%7Bn_z%7D%20%5Cotimes%20D_%7Bn_x%2Cq_x%7D%5E%7BT%7DD_%7Bn_x%2Cq_x%7D%20%2B%20%5Clambda_z%20D_%7Bn_z%2Cq_z%7D%5E%7BT%7DD_%7Bn_z%2Cq_z%7D%20%5Cotimes%20I_%7Bn_x%7D. "P_\lambda = \lambda_x I_{n_z} \otimes D_{n_x,q_x}^{T}D_{n_x,q_x} + \lambda_z D_{n_z,q_z}^{T}D_{n_z,q_z} \otimes I_{n_x}.")

## Installation

You can install the development version of WH from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("GuillaumeBiessy/WH", build_vignettes = TRUE)
```

## How to use the package ?

The `WH` package features two main functions `WH_1d` and `WH_2d`
corresponding to the one-dimension and two-dimension cases respectively.
Two arguments are mandatory for those functions:

- The vector (or matrix in the two-dimension case) `d` corresponding to
  the number of observed events of interest by age (or by age and
  duration in the two-dimension case). `d` should have named elements
  (or rows and columns) for the model results to be extrapolated.

- The vector (or matrix in the two-dimension case) `ec` corresponding to
  the portfolio central exposure by age (or by age and duration in the
  two-dimension case) whose dimensions should match those of `d`. The
  contribution of each individual to the portfolio central exposure
  corresponds to the time the individual was actually observed with
  corresponding age (and duration in the two-dimension cas). It always
  ranges from 0 to 1 and is affected by individuals leaving the
  portfolio, no matter the cause, as well as censoring and truncating
  phenomena.

Additional arguments are described in the documentation of those
functions.

The package also embed two fictive agregated datasets to illustrate how
to use it:

- `portfolio_mortality` contains the agregated number of deaths and
  associated central exposure by age for an annuity portfolio.

- `portfolio_LTC` contains the agregated number of deaths and associated
  central exposure by age and duration (in years) since the onset of LTC
  for the annuitant database of a long-term care portfolio.

``` r
# One-dimension case
keep <- which(portfolio_mort$ec > 0) # observations with no data
d <- portfolio_mort$d[keep]
ec <- portfolio_mort$ec[keep]

WH_1d_fit <- WH_1d(d, ec)
Using outer iteration / Brent method
```

``` r
# Two_dimension case
keep_age <- which(rowSums(portfolio_LTC$ec) > 1e2)
keep_duration <- which(colSums(portfolio_LTC$ec) > 1e2)

d  <- portfolio_LTC$d[keep_age, keep_duration]
ec <- portfolio_LTC$ec[keep_age, keep_duration]

WH_2d_fit <- WH_2d(d, ec)
Using performance iteration / Nelder-Mead method
```

Functions `WH_1d` and `WH_2d` output objects of class `"WH_1d"` and
`"WH_2d"` to which additional functions (including generic S3 methods)
may be applied:

- The `print` function provides a glimpse of the fitted results

``` r
WH_1d_fit
An object fitted using the WH_1D function
Initial data contains 74 data points:
  Observation positions:  19  to  92 
Optimal smoothing parameter selected: 23368 
Associated degrees of freedom: 4.6 
WH_2d_fit
An object fitted using the WH_2D function
Initial data contains 176 data points:
  First  dimension:  74  to  89 
  Second dimension:  0  to  10 
Optimal smoothing parameters selected: 126.2   9.4 
Associated degrees of freedom: 9.2 
```

- The `plot` function generates rough plots of the model fit, the
  associated standard deviation, the model residuals or the associated
  degrees of freedom. See the `plot.WH_1d` and `plot.WH_2d` functions
  help for more details.

``` r
plot(WH_1d_fit)
```

<img src="man/figures/README-plot-1.png" width="100%" style="display: block; margin: auto;" />

``` r
plot(WH_1d_fit, "res")
```

<img src="man/figures/README-plot-2.png" width="100%" style="display: block; margin: auto;" />

``` r
plot(WH_1d_fit, "edf")
```

<img src="man/figures/README-plot-3.png" width="100%" style="display: block; margin: auto;" />

``` r

plot(WH_2d_fit)
```

<img src="man/figures/README-plot-4.png" width="100%" style="display: block; margin: auto;" />

``` r
plot(WH_2d_fit, "std_y_hat")
```

<img src="man/figures/README-plot-5.png" width="100%" style="display: block; margin: auto;" />

- The `predict` function generates an extrapolation of the model. It
  requires a `newdata` argument, a named list with one or two elements
  corresponding to the positions of the new observations. In the
  two-dimension case constraints are used so that the predicted values
  matches the fitted values for the initial observations (see Carballo,
  Durbán, and Lee 2021 to understand why this is required).

``` r
WH_1d_fit |> predict(newdata = 18:99) |> plot()
```

<img src="man/figures/README-predict-1.png" width="100%" style="display: block; margin: auto;" />

``` r
WH_2d_fit |> predict(newdata = list(age = 50:99,
                                    duration = 0:19)) |> plot()
```

<img src="man/figures/README-predict-2.png" width="100%" style="display: block; margin: auto;" />

- Finally the `output_to_df` converts an `"WH_1d"` or `"WH_2d"` object
  into a `data.frame`. Information about the fit is lost in the process.
  This may be useful to produce better visualizations, for example using
  the ggplot2 package.

``` r
WH_1d_df <- WH_1d_fit |> output_to_df()
WH_2d_df <- WH_2d_fit |> output_to_df()
```

## Further WH smoothing theory

### Explicit solution

An explicit solution to the smoothing equation is obtained by computing
the derivate according to
![\theta](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctheta "\theta")
of the minimized quantity in the equation. In the cas
![w](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;w "w")
contains enough non-0 weights:

![\hat{y} = (W + P\_\lambda)^{- 1}Wy.](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Chat%7By%7D%20%3D%20%28W%20%2B%20P_%5Clambda%29%5E%7B-%201%7DWy. "\hat{y} = (W + P_\lambda)^{- 1}Wy.")

### Role of the penalization

In the smoothing equation,
![(y - \theta)^{T}W(y - \theta)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%28y%20-%20%5Ctheta%29%5E%7BT%7DW%28y%20-%20%5Ctheta%29 "(y - \theta)^{T}W(y - \theta)")
is a fidelity criterion
![\theta^{T}P\_\lambda\theta](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctheta%5E%7BT%7DP_%5Clambda%5Ctheta "\theta^{T}P_\lambda\theta")
a regularity criterion. The relative importance of those criterions is
controlled by the smoothing parameter (or parameter vectors in the
two-dimension case)
![\lambda](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Clambda "\lambda").

In the one-dimension case, the penalization matrix may be rewritten:

![\theta^{T}P\_\lambda\theta = \begin{cases}\lambda\underset{i = 1}{\overset{n - 1}{\sum}}(\theta\_{i + 1} - \theta_i)^2 & \text{si }q = 1\\\\ \lambda\underset{i = 1}{\overset{n - 2}{\sum}}(\[\theta\_{i + 2} - \theta\_{i + 1}\] - \[\theta\_{i + 1} - \theta_i\])^2 & \text{si }q = 2
\end{cases}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctheta%5E%7BT%7DP_%5Clambda%5Ctheta%20%3D%20%5Cbegin%7Bcases%7D%5Clambda%5Cunderset%7Bi%20%3D%201%7D%7B%5Coverset%7Bn%20-%201%7D%7B%5Csum%7D%7D%28%5Ctheta_%7Bi%20%2B%201%7D%20-%20%5Ctheta_i%29%5E2%20%26%20%5Ctext%7Bsi%20%7Dq%20%3D%201%5C%5C%20%5Clambda%5Cunderset%7Bi%20%3D%201%7D%7B%5Coverset%7Bn%20-%202%7D%7B%5Csum%7D%7D%28%5B%5Ctheta_%7Bi%20%2B%202%7D%20-%20%5Ctheta_%7Bi%20%2B%201%7D%5D%20-%20%5B%5Ctheta_%7Bi%20%2B%201%7D%20-%20%5Ctheta_i%5D%29%5E2%20%26%20%5Ctext%7Bsi%20%7Dq%20%3D%202%0A%5Cend%7Bcases%7D "\theta^{T}P_\lambda\theta = \begin{cases}\lambda\underset{i = 1}{\overset{n - 1}{\sum}}(\theta_{i + 1} - \theta_i)^2 & \text{si }q = 1\\ \lambda\underset{i = 1}{\overset{n - 2}{\sum}}([\theta_{i + 2} - \theta_{i + 1}] - [\theta_{i + 1} - \theta_i])^2 & \text{si }q = 2
\end{cases}")

It is easily seen that:

- In the case where
  ![q_x = 1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;q_x%20%3D%201 "q_x = 1")
  then this is equal to 0 only in the case where, for all
  ![i\in\[1,n_x - 1\]](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;i%5Cin%5B1%2Cn_x%20-%201%5D "i\in[1,n_x - 1]"),
  ![\theta\_{i + 1} - \theta_i = 0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctheta_%7Bi%20%2B%201%7D%20-%20%5Ctheta_i%20%3D%200 "\theta_{i + 1} - \theta_i = 0")
  or equivalently for all
  ![i\in\[1,n_x\]](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;i%5Cin%5B1%2Cn_x%5D "i\in[1,n_x]"),
  ![\theta_i = \theta_1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctheta_i%20%3D%20%5Ctheta_1 "\theta_i = \theta_1").

- In the case where
  ![q_x = 2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;q_x%20%3D%202 "q_x = 2"),
  this is equal to 0 only in the case where, for all
  ![i\in\[1,n_x - 2\]](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;i%5Cin%5B1%2Cn_x%20-%202%5D "i\in[1,n_x - 2]"),
  ![\theta\_{i + 2} - \theta\_{i + 1} = \theta\_{i + 1} - \theta_i](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctheta_%7Bi%20%2B%202%7D%20-%20%5Ctheta_%7Bi%20%2B%201%7D%20%3D%20%5Ctheta_%7Bi%20%2B%201%7D%20-%20%5Ctheta_i "\theta_{i + 2} - \theta_{i + 1} = \theta_{i + 1} - \theta_i"),
  or equivalently for all
  ![i\in\[1,n_x\]](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;i%5Cin%5B1%2Cn_x%5D "i\in[1,n_x]"),
  ![\theta_i = \theta_1 + (i - 1)(\theta_2 - \theta_1)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctheta_i%20%3D%20%5Ctheta_1%20%2B%20%28i%20-%201%29%28%5Ctheta_2%20-%20%5Ctheta_1%29 "\theta_i = \theta_1 + (i - 1)(\theta_2 - \theta_1)").

In the sense of the penalization of order
![q_x](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;q_x "q_x"),
the space of observation vectors that are totally smooth is therefore
the space of polynomial functions of degrees at most
![q - 1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;q%20-%201 "q - 1").
Those results carry on to two-dimension smoothing where a totally smooth
observation matrix
![Y](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;Y "Y")
corresponds to constant / aligned observations on the rows / columns of
the observation matrix, depending on the values of
![q_x](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;q_x "q_x")
and
![q_z](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;q_z "q_z").

### Confidence / credibility intervals

Let
![B](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;B "B")
a matrix such that
![P\_\lambda = B^TB](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;P_%5Clambda%20%3D%20B%5ETB "P_\lambda = B^TB").
The smoothing equation may be rearranged:

![(y - \theta)^TW(y - \theta) + \theta^TP\_\lambda\theta = \left\Vert\begin{pmatrix}\sqrt{W}y \\\\ 0\end{pmatrix} - \begin{pmatrix}\sqrt{W} \\\\ B\end{pmatrix}\theta\right\Vert^2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%28y%20-%20%5Ctheta%29%5ETW%28y%20-%20%5Ctheta%29%20%2B%20%5Ctheta%5ETP_%5Clambda%5Ctheta%20%3D%20%5Cleft%5CVert%5Cbegin%7Bpmatrix%7D%5Csqrt%7BW%7Dy%20%5C%5C%200%5Cend%7Bpmatrix%7D%20-%20%5Cbegin%7Bpmatrix%7D%5Csqrt%7BW%7D%20%5C%5C%20B%5Cend%7Bpmatrix%7D%5Ctheta%5Cright%5CVert%5E2 "(y - \theta)^TW(y - \theta) + \theta^TP_\lambda\theta = \left\Vert\begin{pmatrix}\sqrt{W}y \\ 0\end{pmatrix} - \begin{pmatrix}\sqrt{W} \\ B\end{pmatrix}\theta\right\Vert^2")

which corresponds to an weighted mean square problem. In particular, if
![y \sim \mathcal{N}(\theta,\sigma^2W^{- 1})](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;y%20%5Csim%20%5Cmathcal%7BN%7D%28%5Ctheta%2C%5Csigma%5E2W%5E%7B-%201%7D%29 "y \sim \mathcal{N}(\theta,\sigma^2W^{- 1})")
then
![\sqrt{W}y \sim \mathcal{N}(\sqrt{W}\theta,\sigma^2I_n)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Csqrt%7BW%7Dy%20%5Csim%20%5Cmathcal%7BN%7D%28%5Csqrt%7BW%7D%5Ctheta%2C%5Csigma%5E2I_n%29 "\sqrt{W}y \sim \mathcal{N}(\sqrt{W}\theta,\sigma^2I_n)")
and the framework of regression may be used. Note that
![\theta = \mathbb{E}(y)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctheta%20%3D%20%5Cmathbb%7BE%7D%28y%29 "\theta = \mathbb{E}(y)")
is in this case both the parameter vector and the underlying law WH
smoothing is trying to estimate.

In the cas of 0 weights, the matrix
![W^{- 1}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;W%5E%7B-%201%7D "W^{- 1}")
is not properly defined. By defining
![w\_\ast](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;w_%5Cast "w_\ast")
the subvector of stricly positive weights and
![y\_\ast](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;y_%5Cast "y_\ast")
the associated observation vector and noting
![W\_\ast = \text{Diag}(w\_\ast)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;W_%5Cast%20%3D%20%5Ctext%7BDiag%7D%28w_%5Cast%29 "W_\ast = \text{Diag}(w_\ast)"),
the previous quantity may be rewritten:

![\left\Vert\begin{pmatrix}\sqrt{W}y \\\\ 0\end{pmatrix} - \begin{pmatrix}\sqrt{W} \\\\ B\end{pmatrix}\theta\right\Vert^2 = \left\Vert\begin{pmatrix}\sqrt{W\_\ast}y\_\ast \\\\ 0\end{pmatrix} - \begin{pmatrix}\sqrt{W\_\ast} \\\\ B\end{pmatrix}\theta\right\Vert^2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cleft%5CVert%5Cbegin%7Bpmatrix%7D%5Csqrt%7BW%7Dy%20%5C%5C%200%5Cend%7Bpmatrix%7D%20-%20%5Cbegin%7Bpmatrix%7D%5Csqrt%7BW%7D%20%5C%5C%20B%5Cend%7Bpmatrix%7D%5Ctheta%5Cright%5CVert%5E2%20%3D%20%5Cleft%5CVert%5Cbegin%7Bpmatrix%7D%5Csqrt%7BW_%5Cast%7Dy_%5Cast%20%5C%5C%200%5Cend%7Bpmatrix%7D%20-%20%5Cbegin%7Bpmatrix%7D%5Csqrt%7BW_%5Cast%7D%20%5C%5C%20B%5Cend%7Bpmatrix%7D%5Ctheta%5Cright%5CVert%5E2 "\left\Vert\begin{pmatrix}\sqrt{W}y \\ 0\end{pmatrix} - \begin{pmatrix}\sqrt{W} \\ B\end{pmatrix}\theta\right\Vert^2 = \left\Vert\begin{pmatrix}\sqrt{W_\ast}y_\ast \\ 0\end{pmatrix} - \begin{pmatrix}\sqrt{W_\ast} \\ B\end{pmatrix}\theta\right\Vert^2")

and the problem may now be expressed properly in terms of
![y\_\ast](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;y_%5Cast "y_\ast")
and
![w\_\ast](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;w_%5Cast "w_\ast").
Nevertheless to keep things simple we stick with
![y](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;y "y")
and
![w](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;w "w")
and consider that
![W^{- 1}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;W%5E%7B-%201%7D "W^{- 1}")
contains infinite values associated with initial 0 weights in
![w](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;w "w").

The previous explicit solution to the smoothing equation garanties the
normality of
![\hat{y}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Chat%7By%7D "\hat{y}").
In the framework of regression, the law of error propagation yields:

![\text{Var}(\hat{y}) = \text{Var}\[(W + P\_\lambda)^{- 1}Wy\] = (W + P\_\lambda)^{- 1}W\text{Var}(y) W(W + P\_\lambda)^{- 1} = \sigma^2 (W + P\_\lambda)^{- 1}W(W + P\_\lambda)^{- 1}.](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctext%7BVar%7D%28%5Chat%7By%7D%29%20%3D%20%5Ctext%7BVar%7D%5B%28W%20%2B%20P_%5Clambda%29%5E%7B-%201%7DWy%5D%20%3D%20%28W%20%2B%20P_%5Clambda%29%5E%7B-%201%7DW%5Ctext%7BVar%7D%28y%29%20W%28W%20%2B%20P_%5Clambda%29%5E%7B-%201%7D%20%3D%20%5Csigma%5E2%20%28W%20%2B%20P_%5Clambda%29%5E%7B-%201%7DW%28W%20%2B%20P_%5Clambda%29%5E%7B-%201%7D. "\text{Var}(\hat{y}) = \text{Var}[(W + P_\lambda)^{- 1}Wy] = (W + P_\lambda)^{- 1}W\text{Var}(y) W(W + P_\lambda)^{- 1} = \sigma^2 (W + P_\lambda)^{- 1}W(W + P_\lambda)^{- 1}.")

Unfortunately
![\mathbb{E}(\hat{y}) = (W + P\_\lambda)^{- 1}W\mathbb{E}(y) \ne \mathbb{E}(y)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbb%7BE%7D%28%5Chat%7By%7D%29%20%3D%20%28W%20%2B%20P_%5Clambda%29%5E%7B-%201%7DW%5Cmathbb%7BE%7D%28y%29%20%5Cne%20%5Cmathbb%7BE%7D%28y%29 "\mathbb{E}(\hat{y}) = (W + P_\lambda)^{- 1}W\mathbb{E}(y) \ne \mathbb{E}(y)")
if
![\lambda \ne 0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Clambda%20%5Cne%200 "\lambda \ne 0").
The penalization introduces a *smoothing bias* and those results may not
be used to compute a confidence interval for
![\theta](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctheta "\theta").

An alternative is to adopt a bayesian interpretation of the smoothing
equation and to interpretate the penalization as an *a priori* of the
form
![\exp(- \theta^{T}P\_\lambda\theta / 2\sigma^2)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cexp%28-%20%5Ctheta%5E%7BT%7DP_%5Clambda%5Ctheta%20%2F%202%5Csigma%5E2%29 "\exp(- \theta^{T}P_\lambda\theta / 2\sigma^2)")
on
![\theta](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctheta "\theta"),
which boils down to assume that
![\theta \sim \mathcal{N}(0, \sigma^2P\_\lambda^{-})](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctheta%20%5Csim%20%5Cmathcal%7BN%7D%280%2C%20%5Csigma%5E2P_%5Clambda%5E%7B-%7D%29 "\theta \sim \mathcal{N}(0, \sigma^2P_\lambda^{-})")
(where
![A^-](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;A%5E- "A^-")
corresponds to the pseudo-inverse of a matrix
![A](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;A "A")).
Bayes formula then yields:

![\begin{aligned}
f(\theta \| y) &\propto f_y(y \| \theta) f\_\theta(\theta) \\\\
&\propto \exp\left(- \frac{1}{2\sigma^2}(y - \theta)^{T}W(y - \theta)\right)\exp\left(-\frac{1}{2\sigma^2}\theta^{T}P\_\lambda\theta\right) \\\\
&\propto \exp\left(- \frac{1}{2\sigma^2}\left\[(y - \theta)^{T}W(y - \theta) +\theta^{T}P\_\lambda\theta\right\]\right)
\end{aligned}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cbegin%7Baligned%7D%0Af%28%5Ctheta%20%7C%20y%29%20%26%5Cpropto%20f_y%28y%20%7C%20%5Ctheta%29%20f_%5Ctheta%28%5Ctheta%29%20%5C%5C%0A%26%5Cpropto%20%5Cexp%5Cleft%28-%20%5Cfrac%7B1%7D%7B2%5Csigma%5E2%7D%28y%20-%20%5Ctheta%29%5E%7BT%7DW%28y%20-%20%5Ctheta%29%5Cright%29%5Cexp%5Cleft%28-%5Cfrac%7B1%7D%7B2%5Csigma%5E2%7D%5Ctheta%5E%7BT%7DP_%5Clambda%5Ctheta%5Cright%29%20%5C%5C%0A%26%5Cpropto%20%5Cexp%5Cleft%28-%20%5Cfrac%7B1%7D%7B2%5Csigma%5E2%7D%5Cleft%5B%28y%20-%20%5Ctheta%29%5E%7BT%7DW%28y%20-%20%5Ctheta%29%20%2B%5Ctheta%5E%7BT%7DP_%5Clambda%5Ctheta%5Cright%5D%5Cright%29%0A%5Cend%7Baligned%7D "\begin{aligned}
f(\theta | y) &\propto f_y(y | \theta) f_\theta(\theta) \\
&\propto \exp\left(- \frac{1}{2\sigma^2}(y - \theta)^{T}W(y - \theta)\right)\exp\left(-\frac{1}{2\sigma^2}\theta^{T}P_\lambda\theta\right) \\
&\propto \exp\left(- \frac{1}{2\sigma^2}\left[(y - \theta)^{T}W(y - \theta) +\theta^{T}P_\lambda\theta\right]\right)
\end{aligned}")

The mode of the posterior distribution
![\hat{\theta} = \text{argmax} \[f(\theta \| y)\]](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Chat%7B%5Ctheta%7D%20%3D%20%5Ctext%7Bargmax%7D%20%5Bf%28%5Ctheta%20%7C%20y%29%5D "\hat{\theta} = \text{argmax} [f(\theta | y)]"),
also known as *maximum a posteriori* (MAP) therefore matches the
solution
![\hat{y}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Chat%7By%7D "\hat{y}")
of the smoothing equation. Furthermore, a Taylor extension of
![\ln f(\theta \| y)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cln%20f%28%5Ctheta%20%7C%20y%29 "\ln f(\theta | y)")
in
![\theta = \hat{y}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctheta%20%3D%20%5Chat%7By%7D "\theta = \hat{y}")
allows the posterior distribution to be recognized as
![\mathcal{N}(\hat{y}, \sigma^2(W + P\_\lambda)^{- 1})](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathcal%7BN%7D%28%5Chat%7By%7D%2C%20%5Csigma%5E2%28W%20%2B%20P_%5Clambda%29%5E%7B-%201%7D%29 "\mathcal{N}(\hat{y}, \sigma^2(W + P_\lambda)^{- 1})").

An unbiased estimator of
![\sigma^2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Csigma%5E2 "\sigma^2")
is then given by:

![\hat{\sigma}^2 = \frac{\Vert\sqrt{W}(y - \hat{y})\Vert^2}{n\_\* - \text{edf}}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Chat%7B%5Csigma%7D%5E2%20%3D%20%5Cfrac%7B%5CVert%5Csqrt%7BW%7D%28y%20-%20%5Chat%7By%7D%29%5CVert%5E2%7D%7Bn_%2A%20-%20%5Ctext%7Bedf%7D%7D "\hat{\sigma}^2 = \frac{\Vert\sqrt{W}(y - \hat{y})\Vert^2}{n_* - \text{edf}}")

where
![n\_\*](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n_%2A "n_*")
corresponds to the number of non 0 weights and
![\text{edf} = \text{tr}(H) = \text{tr}\[(W + P\_\lambda)^{- 1}W\]](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctext%7Bedf%7D%20%3D%20%5Ctext%7Btr%7D%28H%29%20%3D%20%5Ctext%7Btr%7D%5B%28W%20%2B%20P_%5Clambda%29%5E%7B-%201%7DW%5D "\text{edf} = \text{tr}(H) = \text{tr}[(W + P_\lambda)^{- 1}W]").

This result allows
![100(1 -\alpha)\\%](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;100%281%20-%5Calpha%29%5C%25 "100(1 -\alpha)\%")
credibility intervals of the form
![\[\hat{y} \pm t\_{n - \text{edf}}(\frac{1 - \alpha}{2})\sqrt{\text{diag}\left\lbrace(W + P\_\lambda)^{- 1}\right\rbrace}\]](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5B%5Chat%7By%7D%20%5Cpm%20t_%7Bn%20-%20%5Ctext%7Bedf%7D%7D%28%5Cfrac%7B1%20-%20%5Calpha%7D%7B2%7D%29%5Csqrt%7B%5Ctext%7Bdiag%7D%5Cleft%5Clbrace%28W%20%2B%20P_%5Clambda%29%5E%7B-%201%7D%5Cright%5Crbrace%7D%5D "[\hat{y} \pm t_{n - \text{edf}}(\frac{1 - \alpha}{2})\sqrt{\text{diag}\left\lbrace(W + P_\lambda)^{- 1}\right\rbrace}]")
to be computed where
![t\_{n\_\* - \text{edf}}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;t_%7Bn_%2A%20-%20%5Ctext%7Bedf%7D%7D "t_{n_* - \text{edf}}")
is the distribution function for the student law of parameter
![n\_\* - \text{edf}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n_%2A%20-%20%5Ctext%7Bedf%7D "n_* - \text{edf}").
Let us note that the use of student distribution instead of the normal
distribution is linked with the presence of an unknown
![\sigma^2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Csigma%5E2 "\sigma^2")
parameter.

## What should I use as observations and weights ?

The previous section show that applying WH smoothing to a couple
![(y,w)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%28y%2Cw%29 "(y,w)")
such that
![y \sim \mathcal{N}(\theta, \sigma^2W^{- 1})](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;y%20%5Csim%20%5Cmathcal%7BN%7D%28%5Ctheta%2C%20%5Csigma%5E2W%5E%7B-%201%7D%29 "y \sim \mathcal{N}(\theta, \sigma^2W^{- 1})")
where
![\theta](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctheta "\theta")
is the underlying law we want to estimate,
![\sigma^2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Csigma%5E2 "\sigma^2")
is an overdispersion parameter to be estimated, and
![W = \text{Diag}(w)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;W%20%3D%20%5Ctext%7BDiag%7D%28w%29 "W = \text{Diag}(w)"),
then, using a bayesian interpration, credibility intervals were
available for the posterior distribution of
![\theta \| y](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctheta%20%7C%20y "\theta | y").
In this section, in the framework of survival analysis models, we
exhibit a candidate for
![(y,w)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%28y%2Cw%29 "(y,w)").

### Survival analysis framework - one-dimension case

Let us consider the observation of
![m](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;m "m")
individuals in the cas of a longitudinal study under left-truncating and
right-censoring phenomena. Let us assume a single law must be estimated
(for example a mortality law) and that it depends on a single covariate
named
![x](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;x "x")
(for example the age of the insured). This law may be entirely
characterized by providing one of the 3 following quantities:

- The cumulative distribution function
  ![F(x)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;F%28x%29 "F(x)")
  or the survival function
  ![S(x) = 1 - F(x)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;S%28x%29%20%3D%201%20-%20F%28x%29 "S(x) = 1 - F(x)"),

- The associated density function
  ![f(x) = - \frac{\text{d}}{\text{d}x}S(x)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;f%28x%29%20%3D%20-%20%5Cfrac%7B%5Ctext%7Bd%7D%7D%7B%5Ctext%7Bd%7Dx%7DS%28x%29 "f(x) = - \frac{\text{d}}{\text{d}x}S(x)"),

- The hazard rate
  ![\mu(x) = - \frac{\text{d}}{\text{d}x}\text{ln} S(x)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmu%28x%29%20%3D%20-%20%5Cfrac%7B%5Ctext%7Bd%7D%7D%7B%5Ctext%7Bd%7Dx%7D%5Ctext%7Bln%7D%20S%28x%29 "\mu(x) = - \frac{\text{d}}{\text{d}x}\text{ln} S(x)")

Let us assume that the underlying law depends on a parameter vector
![\beta](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cbeta "\beta")
that is to be estimated using maximum likelihood. The likelihood
associated with the observation of the individuals is:

![\mathcal{L}(\beta) = \underset{i = 1}{\overset{m}{\prod}} \left\[\frac{f(x_i + t_i,\beta)}{S(x_i,\beta)}\right\]^{\delta_i}\left\[\frac{S(x_i + t_i,\beta)}{S(x_i,\beta)}\right\]^{1 - \delta_i}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathcal%7BL%7D%28%5Cbeta%29%20%3D%20%5Cunderset%7Bi%20%3D%201%7D%7B%5Coverset%7Bm%7D%7B%5Cprod%7D%7D%20%5Cleft%5B%5Cfrac%7Bf%28x_i%20%2B%20t_i%2C%5Cbeta%29%7D%7BS%28x_i%2C%5Cbeta%29%7D%5Cright%5D%5E%7B%5Cdelta_i%7D%5Cleft%5B%5Cfrac%7BS%28x_i%20%2B%20t_i%2C%5Cbeta%29%7D%7BS%28x_i%2C%5Cbeta%29%7D%5Cright%5D%5E%7B1%20-%20%5Cdelta_i%7D "\mathcal{L}(\beta) = \underset{i = 1}{\overset{m}{\prod}} \left[\frac{f(x_i + t_i,\beta)}{S(x_i,\beta)}\right]^{\delta_i}\left[\frac{S(x_i + t_i,\beta)}{S(x_i,\beta)}\right]^{1 - \delta_i}")

where for each
![x_i](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;x_i "x_i")
represents the age at the start of the observation,
![t_i](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;t_i "t_i")
is the observation duration and
![\delta_i](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cdelta_i "\delta_i")
is 1 if the event of interest (for example death) has been observed and
0 otherwise. Those 3 elements should be computed by taking into account
the date of subscribing and a possible termination date (for example
beauce of policy lapse) for each individual, possible presence of a
waiting period and exclusion of specific periods of time because of
incomplete data or medical underwriting effects. The observation period
may therefore be shorter than the actual period of presence of the
individuals in the portfolio.

Maximization of the previous likelihood may be rewritten by taking the
logarithm and using the relations:

![\begin{aligned}
S(x) & = \exp\left(\underset{u = 0}{\overset{x}{\int}}\mu(u)\text{d}u\right) & f(x) & = \mu(x)S(x)
\end{aligned}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cbegin%7Baligned%7D%0AS%28x%29%20%26%20%3D%20%5Cexp%5Cleft%28%5Cunderset%7Bu%20%3D%200%7D%7B%5Coverset%7Bx%7D%7B%5Cint%7D%7D%5Cmu%28u%29%5Ctext%7Bd%7Du%5Cright%29%20%26%20f%28x%29%20%26%20%3D%20%5Cmu%28x%29S%28x%29%0A%5Cend%7Baligned%7D "\begin{aligned}
S(x) & = \exp\left(\underset{u = 0}{\overset{x}{\int}}\mu(u)\text{d}u\right) & f(x) & = \mu(x)S(x)
\end{aligned}")

which yields, after a few simplifications:

![\ell(\beta) = \underset{i = 1}{\overset{m}{\sum}} \left\[\delta_i \ln\mu(x_i + t_i,\beta) - \underset{u = 0}{\overset{t_i}{\int}}\mu(x_i + u,\beta)\text{d}u\right\]](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cell%28%5Cbeta%29%20%3D%20%5Cunderset%7Bi%20%3D%201%7D%7B%5Coverset%7Bm%7D%7B%5Csum%7D%7D%20%5Cleft%5B%5Cdelta_i%20%5Cln%5Cmu%28x_i%20%2B%20t_i%2C%5Cbeta%29%20-%20%5Cunderset%7Bu%20%3D%200%7D%7B%5Coverset%7Bt_i%7D%7B%5Cint%7D%7D%5Cmu%28x_i%20%2B%20u%2C%5Cbeta%29%5Ctext%7Bd%7Du%5Cright%5D "\ell(\beta) = \underset{i = 1}{\overset{m}{\sum}} \left[\delta_i \ln\mu(x_i + t_i,\beta) - \underset{u = 0}{\overset{t_i}{\int}}\mu(x_i + u,\beta)\text{d}u\right]")

Let us assume that the hazard rate is a piecewise constant function on
one-year interval between integer ages or more formally:
![\mu(x + \epsilon) = \mu(x)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmu%28x%20%2B%20%5Cepsilon%29%20%3D%20%5Cmu%28x%29 "\mu(x + \epsilon) = \mu(x)")
for all
![x \in \mathbb{N}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;x%20%5Cin%20%5Cmathbb%7BN%7D "x \in \mathbb{N}")
and
![\epsilon \in \[0,1\[](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cepsilon%20%5Cin%20%5B0%2C1%5B "\epsilon \in [0,1[").

Let us further note that if
![\mathbf{1}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7B1%7D "\mathbf{1}")
is the index function, then for all
![0 \le a \< x\_{\max}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;0%20%5Cle%20a%20%3C%20x_%7B%5Cmax%7D "0 \le a < x_{\max}"),
![\underset{x = x\_{\min}}{\overset{x\_{\max}}{\sum}} \mathbf{1}(x \le a \< x + 1) = 1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cunderset%7Bx%20%3D%20x_%7B%5Cmin%7D%7D%7B%5Coverset%7Bx_%7B%5Cmax%7D%7D%7B%5Csum%7D%7D%20%5Cmathbf%7B1%7D%28x%20%5Cle%20a%20%3C%20x%20%2B%201%29%20%3D%201 "\underset{x = x_{\min}}{\overset{x_{\max}}{\sum}} \mathbf{1}(x \le a < x + 1) = 1")
by noting
![x\_{\min} = \min(x)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;x_%7B%5Cmin%7D%20%3D%20%5Cmin%28x%29 "x_{\min} = \min(x)")
and
![x\_{\max} = \max(x)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;x_%7B%5Cmax%7D%20%3D%20%5Cmax%28x%29 "x_{\max} = \max(x)").
The previous likelihood therefore becomes:

![\ell(\beta) = \underset{i = 1}{\overset{m}{\sum}} \left\[\underset{x = x\_{\min}}{\overset{x\_{\max}}{\sum}} \delta_i\mathbf{1}(x \le x_i + t_i \< x + 1)  \ln\mu(x_i + t_i,\beta) - \underset{u = 0}{\overset{t_i}{\int}}\underset{x = x\_{\min}}{\overset{x\_{\max}}{\sum}} \mathbf{1}(x \le x_i + u \< x + 1)\mu(x_i + u,\beta)\text{d}u\right\]](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cell%28%5Cbeta%29%20%3D%20%5Cunderset%7Bi%20%3D%201%7D%7B%5Coverset%7Bm%7D%7B%5Csum%7D%7D%20%5Cleft%5B%5Cunderset%7Bx%20%3D%20x_%7B%5Cmin%7D%7D%7B%5Coverset%7Bx_%7B%5Cmax%7D%7D%7B%5Csum%7D%7D%20%5Cdelta_i%5Cmathbf%7B1%7D%28x%20%5Cle%20x_i%20%2B%20t_i%20%3C%20x%20%2B%201%29%20%20%5Cln%5Cmu%28x_i%20%2B%20t_i%2C%5Cbeta%29%20-%20%5Cunderset%7Bu%20%3D%200%7D%7B%5Coverset%7Bt_i%7D%7B%5Cint%7D%7D%5Cunderset%7Bx%20%3D%20x_%7B%5Cmin%7D%7D%7B%5Coverset%7Bx_%7B%5Cmax%7D%7D%7B%5Csum%7D%7D%20%5Cmathbf%7B1%7D%28x%20%5Cle%20x_i%20%2B%20u%20%3C%20x%20%2B%201%29%5Cmu%28x_i%20%2B%20u%2C%5Cbeta%29%5Ctext%7Bd%7Du%5Cright%5D "\ell(\beta) = \underset{i = 1}{\overset{m}{\sum}} \left[\underset{x = x_{\min}}{\overset{x_{\max}}{\sum}} \delta_i\mathbf{1}(x \le x_i + t_i < x + 1)  \ln\mu(x_i + t_i,\beta) - \underset{u = 0}{\overset{t_i}{\int}}\underset{x = x_{\min}}{\overset{x_{\max}}{\sum}} \mathbf{1}(x \le x_i + u < x + 1)\mu(x_i + u,\beta)\text{d}u\right]")

The piecewise constant assumption yields
![\mathbf{1}(x \le x_i + t_i \< x + 1) \ln\mu(x_i + t_i,\beta) = \mathbf{1}(x \le x_i + t_i \< x + 1) \ln\mu(x,\beta)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7B1%7D%28x%20%5Cle%20x_i%20%2B%20t_i%20%3C%20x%20%2B%201%29%20%5Cln%5Cmu%28x_i%20%2B%20t_i%2C%5Cbeta%29%20%3D%20%5Cmathbf%7B1%7D%28x%20%5Cle%20x_i%20%2B%20t_i%20%3C%20x%20%2B%201%29%20%5Cln%5Cmu%28x%2C%5Cbeta%29 "\mathbf{1}(x \le x_i + t_i < x + 1) \ln\mu(x_i + t_i,\beta) = \mathbf{1}(x \le x_i + t_i < x + 1) \ln\mu(x,\beta)")
and
![\mathbf{1}(x \le x_i + u \< x + 1)\mu(x_i + u,\beta) = \mathbf{1}(x \le x_i + u \< x + 1) \ln\mu(x,\beta)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7B1%7D%28x%20%5Cle%20x_i%20%2B%20u%20%3C%20x%20%2B%201%29%5Cmu%28x_i%20%2B%20u%2C%5Cbeta%29%20%3D%20%5Cmathbf%7B1%7D%28x%20%5Cle%20x_i%20%2B%20u%20%3C%20x%20%2B%201%29%20%5Cln%5Cmu%28x%2C%5Cbeta%29 "\mathbf{1}(x \le x_i + u < x + 1)\mu(x_i + u,\beta) = \mathbf{1}(x \le x_i + u < x + 1) \ln\mu(x,\beta)").
The two sums may then be interverted, giving:

![\begin{aligned}
\ell(\beta) &= \underset{x = x\_{\min}}{\overset{x\_{\max}}{\sum}} \left\[\ln\mu(x,\beta) d(x) - \mu(x,\beta) e_c(x)\right\] \quad \text{where} \\\\
d(x) & = \underset{i = 1}{\overset{m}{\sum}} \delta_i \mathbf{1}(x \le x_i + t_i \< x + 1)  \quad \text{and} \\\\
e_c(x) & = \underset{i = 1}{\overset{m}{\sum}}\underset{u = 0}{\overset{t_i}{\int}}\mathbf{1}(x \le x_i + u \< x + 1)\text{d}u = \underset{i = 1}{\overset{m}{\sum}} \left\[\min(t_i, x - x_i + 1) - \max(0, x - x_i)\right\]^+
\end{aligned}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cbegin%7Baligned%7D%0A%5Cell%28%5Cbeta%29%20%26%3D%20%5Cunderset%7Bx%20%3D%20x_%7B%5Cmin%7D%7D%7B%5Coverset%7Bx_%7B%5Cmax%7D%7D%7B%5Csum%7D%7D%20%5Cleft%5B%5Cln%5Cmu%28x%2C%5Cbeta%29%20d%28x%29%20-%20%5Cmu%28x%2C%5Cbeta%29%20e_c%28x%29%5Cright%5D%20%5Cquad%20%5Ctext%7Bwhere%7D%20%5C%5C%0Ad%28x%29%20%26%20%3D%20%5Cunderset%7Bi%20%3D%201%7D%7B%5Coverset%7Bm%7D%7B%5Csum%7D%7D%20%5Cdelta_i%20%5Cmathbf%7B1%7D%28x%20%5Cle%20x_i%20%2B%20t_i%20%3C%20x%20%2B%201%29%20%20%5Cquad%20%5Ctext%7Band%7D%20%5C%5C%0Ae_c%28x%29%20%26%20%3D%20%5Cunderset%7Bi%20%3D%201%7D%7B%5Coverset%7Bm%7D%7B%5Csum%7D%7D%5Cunderset%7Bu%20%3D%200%7D%7B%5Coverset%7Bt_i%7D%7B%5Cint%7D%7D%5Cmathbf%7B1%7D%28x%20%5Cle%20x_i%20%2B%20u%20%3C%20x%20%2B%201%29%5Ctext%7Bd%7Du%20%3D%20%5Cunderset%7Bi%20%3D%201%7D%7B%5Coverset%7Bm%7D%7B%5Csum%7D%7D%20%5Cleft%5B%5Cmin%28t_i%2C%20x%20-%20x_i%20%2B%201%29%20-%20%5Cmax%280%2C%20x%20-%20x_i%29%5Cright%5D%5E%2B%0A%5Cend%7Baligned%7D "\begin{aligned}
\ell(\beta) &= \underset{x = x_{\min}}{\overset{x_{\max}}{\sum}} \left[\ln\mu(x,\beta) d(x) - \mu(x,\beta) e_c(x)\right] \quad \text{where} \\
d(x) & = \underset{i = 1}{\overset{m}{\sum}} \delta_i \mathbf{1}(x \le x_i + t_i < x + 1)  \quad \text{and} \\
e_c(x) & = \underset{i = 1}{\overset{m}{\sum}}\underset{u = 0}{\overset{t_i}{\int}}\mathbf{1}(x \le x_i + u < x + 1)\text{d}u = \underset{i = 1}{\overset{m}{\sum}} \left[\min(t_i, x - x_i + 1) - \max(0, x - x_i)\right]^+
\end{aligned}")

where
![d(x)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;d%28x%29 "d(x)")
and
![e_c(x)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;e_c%28x%29 "e_c(x)")
corresponds respectively to the number of observed deaths between age
![x](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;x "x")
and
![x + 1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;x%20%2B%201 "x + 1")
and to the sum of the observation duration between those dates, by
noting
![a^+ = \max(a, 0)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;a%5E%2B%20%3D%20%5Cmax%28a%2C%200%29 "a^+ = \max(a, 0)").

### Extension to the two-dimension case

The extension of the previous approach to the two-dimension case only
requires minor adjustements to the previous proposition. Let us note
![z\_{\min} = \min(z)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;z_%7B%5Cmin%7D%20%3D%20%5Cmin%28z%29 "z_{\min} = \min(z)")
and
![z\_\text{max} = \max(z)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;z_%5Ctext%7Bmax%7D%20%3D%20%5Cmax%28z%29 "z_\text{max} = \max(z)").
The piecewise constant assumption needs to be extended to each of the
two dimensions. Formally we now assume that
![\mu(x + \epsilon, z + \xi) = \mu(x, z)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmu%28x%20%2B%20%5Cepsilon%2C%20z%20%2B%20%5Cxi%29%20%3D%20%5Cmu%28x%2C%20z%29 "\mu(x + \epsilon, z + \xi) = \mu(x, z)")
for all
![x, z \in \mathbb{N}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;x%2C%20z%20%5Cin%20%5Cmathbb%7BN%7D "x, z \in \mathbb{N}")
and
![\epsilon, \xi \in \[0,1\[](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cepsilon%2C%20%5Cxi%20%5Cin%20%5B0%2C1%5B "\epsilon, \xi \in [0,1[").
The sums on
![x](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;x "x")
are replaced by double sums on the values of both
![x](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;x "x")
and
![z](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;z "z")
and the likelihood becomes:

![\begin{aligned}
\ell(\beta) &= \underset{x = x\_{\min}}{\overset{x\_{\max}}{\sum}} \underset{z = z\_{\min}}{\overset{z\_{\max}}{\sum}}\left\[\ln\mu(x,z,\beta) d(x,z) - \mu(x,z,\beta) e_c(x,z)\right\] \quad \text{where} \\\\
d(x,z) & = \underset{i = 1}{\overset{m}{\sum}} \delta_i \mathbf{1}(x \le x_i + t_i \< x + 1) \mathbf{1}(z \le z_i + t_i \< z + 1)  \quad \text{and}\\\\
e_c(x,z) & = \underset{i = 1}{\overset{m}{\sum}}\underset{u = 0}{\overset{t_i}{\int}}\mathbf{1}(x \le x_i + u \< x + 1)\mathbf{1}(z \le z_i + u \< z + 1)\text{d}u \\\\
& = \underset{i = 1}{\overset{m}{\sum}} \left\[\min(t_i, x + 1 - x_i, z + 1 - z_i) - \max(0, x - x_i, z - z_i)\right\]^+.
\end{aligned}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cbegin%7Baligned%7D%0A%5Cell%28%5Cbeta%29%20%26%3D%20%5Cunderset%7Bx%20%3D%20x_%7B%5Cmin%7D%7D%7B%5Coverset%7Bx_%7B%5Cmax%7D%7D%7B%5Csum%7D%7D%20%5Cunderset%7Bz%20%3D%20z_%7B%5Cmin%7D%7D%7B%5Coverset%7Bz_%7B%5Cmax%7D%7D%7B%5Csum%7D%7D%5Cleft%5B%5Cln%5Cmu%28x%2Cz%2C%5Cbeta%29%20d%28x%2Cz%29%20-%20%5Cmu%28x%2Cz%2C%5Cbeta%29%20e_c%28x%2Cz%29%5Cright%5D%20%5Cquad%20%5Ctext%7Bwhere%7D%20%5C%5C%0Ad%28x%2Cz%29%20%26%20%3D%20%5Cunderset%7Bi%20%3D%201%7D%7B%5Coverset%7Bm%7D%7B%5Csum%7D%7D%20%5Cdelta_i%20%5Cmathbf%7B1%7D%28x%20%5Cle%20x_i%20%2B%20t_i%20%3C%20x%20%2B%201%29%20%5Cmathbf%7B1%7D%28z%20%5Cle%20z_i%20%2B%20t_i%20%3C%20z%20%2B%201%29%20%20%5Cquad%20%5Ctext%7Band%7D%5C%5C%0Ae_c%28x%2Cz%29%20%26%20%3D%20%5Cunderset%7Bi%20%3D%201%7D%7B%5Coverset%7Bm%7D%7B%5Csum%7D%7D%5Cunderset%7Bu%20%3D%200%7D%7B%5Coverset%7Bt_i%7D%7B%5Cint%7D%7D%5Cmathbf%7B1%7D%28x%20%5Cle%20x_i%20%2B%20u%20%3C%20x%20%2B%201%29%5Cmathbf%7B1%7D%28z%20%5Cle%20z_i%20%2B%20u%20%3C%20z%20%2B%201%29%5Ctext%7Bd%7Du%20%5C%5C%0A%26%20%3D%20%5Cunderset%7Bi%20%3D%201%7D%7B%5Coverset%7Bm%7D%7B%5Csum%7D%7D%20%5Cleft%5B%5Cmin%28t_i%2C%20x%20%2B%201%20-%20x_i%2C%20z%20%2B%201%20-%20z_i%29%20-%20%5Cmax%280%2C%20x%20-%20x_i%2C%20z%20-%20z_i%29%5Cright%5D%5E%2B.%0A%5Cend%7Baligned%7D "\begin{aligned}
\ell(\beta) &= \underset{x = x_{\min}}{\overset{x_{\max}}{\sum}} \underset{z = z_{\min}}{\overset{z_{\max}}{\sum}}\left[\ln\mu(x,z,\beta) d(x,z) - \mu(x,z,\beta) e_c(x,z)\right] \quad \text{where} \\
d(x,z) & = \underset{i = 1}{\overset{m}{\sum}} \delta_i \mathbf{1}(x \le x_i + t_i < x + 1) \mathbf{1}(z \le z_i + t_i < z + 1)  \quad \text{and}\\
e_c(x,z) & = \underset{i = 1}{\overset{m}{\sum}}\underset{u = 0}{\overset{t_i}{\int}}\mathbf{1}(x \le x_i + u < x + 1)\mathbf{1}(z \le z_i + u < z + 1)\text{d}u \\
& = \underset{i = 1}{\overset{m}{\sum}} \left[\min(t_i, x + 1 - x_i, z + 1 - z_i) - \max(0, x - x_i, z - z_i)\right]^+.
\end{aligned}")

### Likelihood equations

The log-likelihood in the one-dimension or two-dimension cases may be
expressed on the common vectorial form
![\ell(\beta) = \ln\mu(\beta)^{T}d - \mu(\beta)^{T}e_c](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cell%28%5Cbeta%29%20%3D%20%5Cln%5Cmu%28%5Cbeta%29%5E%7BT%7Dd%20-%20%5Cmu%28%5Cbeta%29%5E%7BT%7De_c "\ell(\beta) = \ln\mu(\beta)^{T}d - \mu(\beta)^{T}e_c")
where
![d](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;d "d")
and
![e_c](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;e_c "e_c")
corresponds respectively to the vector of expected events and associated
central exposure.

In the particular case of log-linear models, the hazard rate may be
defined as
![\ln\mu(\beta) = X\beta](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cln%5Cmu%28%5Cbeta%29%20%3D%20X%5Cbeta "\ln\mu(\beta) = X\beta")
with
![X](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;X "X")
a matrix of dimensions
![n \times p](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n%20%5Ctimes%20p "n \times p")
and full-rank
![p \le n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;p%20%5Cle%20n "p \le n").
The use of the logarithmic link ensures that
![\mu = \exp(X\beta)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmu%20%3D%20%5Cexp%28X%5Cbeta%29 "\mu = \exp(X\beta)")
is always positive. Derivatives of the log-likelihood function are for
this model:

![\frac{\partial \ell}{\partial \beta} = X^{T}\left\[d - \exp(X\beta) \odot e_c\right\] \quad \text{and} \quad \frac{\partial^2 \ell}{\partial\beta^2} = - X^{T}W\_{\beta}X \quad \text{where} \quad W\_{\beta} = \text{Diag}(\exp(X\beta) \odot e_c).](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cfrac%7B%5Cpartial%20%5Cell%7D%7B%5Cpartial%20%5Cbeta%7D%20%3D%20X%5E%7BT%7D%5Cleft%5Bd%20-%20%5Cexp%28X%5Cbeta%29%20%5Codot%20e_c%5Cright%5D%20%5Cquad%20%5Ctext%7Band%7D%20%5Cquad%20%5Cfrac%7B%5Cpartial%5E2%20%5Cell%7D%7B%5Cpartial%5Cbeta%5E2%7D%20%3D%20-%20X%5E%7BT%7DW_%7B%5Cbeta%7DX%20%5Cquad%20%5Ctext%7Bwhere%7D%20%5Cquad%20W_%7B%5Cbeta%7D%20%3D%20%5Ctext%7BDiag%7D%28%5Cexp%28X%5Cbeta%29%20%5Codot%20e_c%29. "\frac{\partial \ell}{\partial \beta} = X^{T}\left[d - \exp(X\beta) \odot e_c\right] \quad \text{and} \quad \frac{\partial^2 \ell}{\partial\beta^2} = - X^{T}W_{\beta}X \quad \text{where} \quad W_{\beta} = \text{Diag}(\exp(X\beta) \odot e_c).")

Let us note that those likelihood equations are exactly those that would
have been obtained by treating the central exposure as a deterministic
quantity and by assuming that the number of deaths follows a Poisson GLM
of parameter
![\mu(\beta)\times e_c](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmu%28%5Cbeta%29%5Ctimes%20e_c "\mu(\beta)\times e_c").
The model above therefore behave as a Poisson GLM (John Ashworth Nelder
and Wedderburn 1972).

### Consequences for WH smoothing

Whittaker-Henderson focuses on the particular case where
![X = I_n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;X%20%3D%20I_n "X = I_n")
and
![\beta = \theta](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cbeta%20%3D%20%5Ctheta "\beta = \theta").
In that case the previous equations yield an explicit solution
![\beta = \ln(d) - \ln(e_c)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cbeta%20%3D%20%5Cln%28d%29%20-%20%5Cln%28e_c%29 "\beta = \ln(d) - \ln(e_c)")
and thus
![\mu = d / e_c](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmu%20%3D%20d%20%2F%20e_c "\mu = d / e_c").

Using the asymptotical properties of the maximum likelihood estimator,
we obtain
![\hat{\beta} \sim \mathcal{N}(\beta, W\_{\hat{\beta}}\\,^{- 1})](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Chat%7B%5Cbeta%7D%20%5Csim%20%5Cmathcal%7BN%7D%28%5Cbeta%2C%20W_%7B%5Chat%7B%5Cbeta%7D%7D%5C%2C%5E%7B-%201%7D%29 "\hat{\beta} \sim \mathcal{N}(\beta, W_{\hat{\beta}}\,^{- 1})")
where diagonal elements of
![W\_{\hat{\beta}}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;W_%7B%5Chat%7B%5Cbeta%7D%7D "W_{\hat{\beta}}")
are simply
![e_c \exp(\hat{\beta}) = e_c d / e_c = d](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;e_c%20%5Cexp%28%5Chat%7B%5Cbeta%7D%29%20%3D%20e_c%20d%20%2F%20e_c%20%3D%20d "e_c \exp(\hat{\beta}) = e_c d / e_c = d").

Thus it has been shown that in the survival analysis framework
introduced, asymptotically
![\ln(d / e_c) \sim \mathcal{N}(\ln(\mu), W^{- 1})](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cln%28d%20%2F%20e_c%29%20%5Csim%20%5Cmathcal%7BN%7D%28%5Cln%28%5Cmu%29%2C%20W%5E%7B-%201%7D%29 "\ln(d / e_c) \sim \mathcal{N}(\ln(\mu), W^{- 1})")
where
![W = \text{Diag}(d)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;W%20%3D%20%5Ctext%7BDiag%7D%28d%29 "W = \text{Diag}(d)").
This justify applying WH smoothing to the observation vector
![y = \ln(d / e_c)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;y%20%3D%20%5Cln%28d%20%2F%20e_c%29 "y = \ln(d / e_c)")
with weights
![w = d](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;w%20%3D%20d "w = d").
The overdispersion parameter is in this cas simply
![\sigma^2 = 1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Csigma%5E2%20%3D%201 "\sigma^2 = 1").
We obtain credibility intervals for the posterior distribution
![\theta \| y](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctheta%20%7C%20y "\theta | y")
of the form:
![\[\hat{y} \pm \Phi(\frac{1 - \alpha}{2})\sqrt{\text{diag}\left\lbrace(W + P\_\lambda)^{- 1}\right\rbrace}\]](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5B%5Chat%7By%7D%20%5Cpm%20%5CPhi%28%5Cfrac%7B1%20-%20%5Calpha%7D%7B2%7D%29%5Csqrt%7B%5Ctext%7Bdiag%7D%5Cleft%5Clbrace%28W%20%2B%20P_%5Clambda%29%5E%7B-%201%7D%5Cright%5Crbrace%7D%5D "[\hat{y} \pm \Phi(\frac{1 - \alpha}{2})\sqrt{\text{diag}\left\lbrace(W + P_\lambda)^{- 1}\right\rbrace}]").

### Generalization to penalized maximum likelihood

The previous approach relies on the asymptotic properties of the maximum
likelihood estimator. When few data is available, those properties may
not hold. An alternative approach is to apply the penalization from the
WH smoothing directly to the previous likelihood function and thus
maximize the penalizaed likelihood
![\ell_P(\beta) = \ell(\beta) - \beta^{T}P\_\lambda\beta / 2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cell_P%28%5Cbeta%29%20%3D%20%5Cell%28%5Cbeta%29%20-%20%5Cbeta%5E%7BT%7DP_%5Clambda%5Cbeta%20%2F%202 "\ell_P(\beta) = \ell(\beta) - \beta^{T}P_\lambda\beta / 2").
This can be seen as a generalization of WH smoothing to non-gaussian
likelihood. Still assuming a log-linear model is used, derivatives of
the log-likelihood functions become:

![\frac{\partial \ell_P}{\partial \beta} = X^{T}\left\[d - \exp(X\beta) \odot e_c\right\] - P\_\lambda\beta \quad \text{and} \quad \frac{\partial^2 \ell_P}{\partial\beta^2} = - (X^{T}W\_{\beta}X + P\_\lambda) \quad \text{where} \quad W\_{\beta} = \text{Diag}(\exp(X\beta) \odot e_c).](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cfrac%7B%5Cpartial%20%5Cell_P%7D%7B%5Cpartial%20%5Cbeta%7D%20%3D%20X%5E%7BT%7D%5Cleft%5Bd%20-%20%5Cexp%28X%5Cbeta%29%20%5Codot%20e_c%5Cright%5D%20-%20P_%5Clambda%5Cbeta%20%5Cquad%20%5Ctext%7Band%7D%20%5Cquad%20%5Cfrac%7B%5Cpartial%5E2%20%5Cell_P%7D%7B%5Cpartial%5Cbeta%5E2%7D%20%3D%20-%20%28X%5E%7BT%7DW_%7B%5Cbeta%7DX%20%2B%20P_%5Clambda%29%20%5Cquad%20%5Ctext%7Bwhere%7D%20%5Cquad%20W_%7B%5Cbeta%7D%20%3D%20%5Ctext%7BDiag%7D%28%5Cexp%28X%5Cbeta%29%20%5Codot%20e_c%29. "\frac{\partial \ell_P}{\partial \beta} = X^{T}\left[d - \exp(X\beta) \odot e_c\right] - P_\lambda\beta \quad \text{and} \quad \frac{\partial^2 \ell_P}{\partial\beta^2} = - (X^{T}W_{\beta}X + P_\lambda) \quad \text{where} \quad W_{\beta} = \text{Diag}(\exp(X\beta) \odot e_c).")

Unlike the previously encountered likelihood, those equations does not
have an explicit solution, even when
![X = I_n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;X%20%3D%20I_n "X = I_n"),
because both
![X\beta](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;X%5Cbeta "X\beta")
and
![\exp(X\beta)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cexp%28X%5Cbeta%29 "\exp(X\beta)")
appear in the equations. Using Newton algorithm, a series of estimators
![(\beta_k)\_{k \ge 0}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%28%5Cbeta_k%29_%7Bk%20%5Cge%200%7D "(\beta_k)_{k \ge 0}")
may be built so that it converges to the penalized maximum likelihood
estimator
![\hat{\beta} = \underset{\beta}{\text{argmin}}\\:l_P(\beta)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Chat%7B%5Cbeta%7D%20%3D%20%5Cunderset%7B%5Cbeta%7D%7B%5Ctext%7Bargmin%7D%7D%5C%3Al_P%28%5Cbeta%29 "\hat{\beta} = \underset{\beta}{\text{argmin}}\:l_P(\beta)").

Those estimators are defined by:

![\begin{aligned}
\beta\_{k + 1} &= \beta_k - \left(\left.\frac{\partial^2 \ell_P}{\partial\beta^2}\right\|\_{\beta = \beta_k}\right)^{- 1} \left.\frac{\partial \ell_P}{\partial\beta}\right\|\_{\beta = \beta_k} \\\\ 
&= \beta_k + (X^{T}W_kX + P\_\lambda)^{- 1} \left\[X^{T}\left(d - \exp(X\beta_k) \odot e_c\right) - P\_\lambda \beta_k\right\] \\\\ 
&= \Psi_k X^{T}W_k z_k
\end{aligned}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cbegin%7Baligned%7D%0A%5Cbeta_%7Bk%20%2B%201%7D%20%26%3D%20%5Cbeta_k%20-%20%5Cleft%28%5Cleft.%5Cfrac%7B%5Cpartial%5E2%20%5Cell_P%7D%7B%5Cpartial%5Cbeta%5E2%7D%5Cright%7C_%7B%5Cbeta%20%3D%20%5Cbeta_k%7D%5Cright%29%5E%7B-%201%7D%20%5Cleft.%5Cfrac%7B%5Cpartial%20%5Cell_P%7D%7B%5Cpartial%5Cbeta%7D%5Cright%7C_%7B%5Cbeta%20%3D%20%5Cbeta_k%7D%20%5C%5C%20%0A%26%3D%20%5Cbeta_k%20%2B%20%28X%5E%7BT%7DW_kX%20%2B%20P_%5Clambda%29%5E%7B-%201%7D%20%5Cleft%5BX%5E%7BT%7D%5Cleft%28d%20-%20%5Cexp%28X%5Cbeta_k%29%20%5Codot%20e_c%5Cright%29%20-%20P_%5Clambda%20%5Cbeta_k%5Cright%5D%20%5C%5C%20%0A%26%3D%20%5CPsi_k%20X%5E%7BT%7DW_k%20z_k%0A%5Cend%7Baligned%7D "\begin{aligned}
\beta_{k + 1} &= \beta_k - \left(\left.\frac{\partial^2 \ell_P}{\partial\beta^2}\right|_{\beta = \beta_k}\right)^{- 1} \left.\frac{\partial \ell_P}{\partial\beta}\right|_{\beta = \beta_k} \\ 
&= \beta_k + (X^{T}W_kX + P_\lambda)^{- 1} \left[X^{T}\left(d - \exp(X\beta_k) \odot e_c\right) - P_\lambda \beta_k\right] \\ 
&= \Psi_k X^{T}W_k z_k
\end{aligned}")

by noting
![\eta_k = X\beta_k](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ceta_k%20%3D%20X%5Cbeta_k "\eta_k = X\beta_k"),
![\Psi_k = (X^{T}W_kX + P\_\lambda)^{- 1}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5CPsi_k%20%3D%20%28X%5E%7BT%7DW_kX%20%2B%20P_%5Clambda%29%5E%7B-%201%7D "\Psi_k = (X^{T}W_kX + P_\lambda)^{- 1}"),
![W_k = \text{Diag}(\exp(\eta_k) \odot e_c)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;W_k%20%3D%20%5Ctext%7BDiag%7D%28%5Cexp%28%5Ceta_k%29%20%5Codot%20e_c%29 "W_k = \text{Diag}(\exp(\eta_k) \odot e_c)")
and
![z_k = \eta_k + W_k^{- 1}\[d - \exp(\eta_k) \odot e_c\]](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;z_k%20%3D%20%5Ceta_k%20%2B%20W_k%5E%7B-%201%7D%5Bd%20-%20%5Cexp%28%5Ceta_k%29%20%5Codot%20e_c%5D "z_k = \eta_k + W_k^{- 1}[d - \exp(\eta_k) \odot e_c]").
Setting
![\eta_0 = \ln d - \ln e_c](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ceta_0%20%3D%20%5Cln%20d%20-%20%5Cln%20e_c "\eta_0 = \ln d - \ln e_c")
yields an adequate starting point to the algorithm (the starting
![\beta_0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cbeta_0 "\beta_0")
does not have to be provided). This implies that
![W_0 = \text{Diag}(d)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;W_0%20%3D%20%5Ctext%7BDiag%7D%28d%29 "W_0 = \text{Diag}(d)")
and
![\mathbf{z}\_0 = \ln d - \ln e_c](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7Bz%7D_0%20%3D%20%5Cln%20d%20-%20%5Cln%20e_c "\mathbf{z}_0 = \ln d - \ln e_c")
which corresponds to the observations and weights used in the regression
framework. The update of
![\beta_k](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cbeta_k "\beta_k")
in the Newton optimization step boils down to applying WH smoothing to
the couple
(![z_k](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;z_k "z_k"),
![W_k](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;W_k "W_k")).
The first estimator
![\beta_1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cbeta_1 "\beta_1")
computed is therefore the solution of WH smoothing in the regression
framework. Further iterations relies on the *working vector*
![z_k](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;z_k "z_k")
and associated weights
![W_k](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;W_k "W_k")
that are adjusted at each step based on the previous step results.

The posterior distribution of
![\beta \| y](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cbeta%20%7C%20y "\beta | y")
may be asymptotically approached by
![\mathcal{N}(\hat{\beta}, (X^TW\_{\hat{\beta}}X + P\_\lambda)^{- 1})](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathcal%7BN%7D%28%5Chat%7B%5Cbeta%7D%2C%20%28X%5ETW_%7B%5Chat%7B%5Cbeta%7D%7DX%20%2B%20P_%5Clambda%29%5E%7B-%201%7D%29 "\mathcal{N}(\hat{\beta}, (X^TW_{\hat{\beta}}X + P_\lambda)^{- 1})"),
which allows
![100(1 -\alpha)\\%](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;100%281%20-%5Calpha%29%5C%25 "100(1 -\alpha)\%")
credibility intervales of the form
![\left\[\hat{\beta} \pm \Phi(1 - \alpha / 2) \sqrt{\text{diag}\lbrace(X^TW\_{\hat{\beta}}X + P\_\lambda)^{- 1}\rbrace}\right\]](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cleft%5B%5Chat%7B%5Cbeta%7D%20%5Cpm%20%5CPhi%281%20-%20%5Calpha%20%2F%202%29%20%5Csqrt%7B%5Ctext%7Bdiag%7D%5Clbrace%28X%5ETW_%7B%5Chat%7B%5Cbeta%7D%7DX%20%2B%20P_%5Clambda%29%5E%7B-%201%7D%5Crbrace%7D%5Cright%5D "\left[\hat{\beta} \pm \Phi(1 - \alpha / 2) \sqrt{\text{diag}\lbrace(X^TW_{\hat{\beta}}X + P_\lambda)^{- 1}\rbrace}\right]")
where
![\Phi](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5CPhi "\Phi")
is the cumulative distribution function of the normal distribution.

## How is the optimal smoothing parameter determined ?

### Short answer

The optimal smoothing parameter is determined according to a statistical
criterion. There are two main types of criteria that are adequate in
this case:

- Prediction error criteria aim at minimizing the (asymptotic)
  prediction error. Such criteria include Akkake Information Criterion
  (AIC), Ordinary Cross-Validation (OCV) and Global Cross-Validation
  (GCV). The Bayesian Information Criterion (BIC) is very close to AIC
  and therefore may be added to this category although its
  interpretation is very different.

- Likelihood-based criteria such as profile likelihood or restricted
  maximum likelihood (REML), a variant that accounts for the complexity
  of the model. Such criteria, while being less interesting in the
  asymptotic case, proved to perform better in most real-life situations
  with finite size samples.

At the time of writing the WH package allows REML (the default), AIC,
BIC and GCV criteria to be selected. The `optimize` function is used in
the one-dimension case (see Brent 1973) and the `optim` function with
the Nelder-Mead algorithm (see John A. Nelder and Mead 1965) in the
two-dimension case. Both functions come from the stats package and
perform adequate optimization in the absence of derivatives.

### Long answer

See the upcoming paper !

## References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-brent1973optimize" class="csl-entry">

Brent, Richard P. 1973. “Algorithms for Minimization Without
Derivatives, Chap. 4.” Prentice-Hall, Englewood Cliffs, NJ.

</div>

<div id="ref-carballo2021prediction" class="csl-entry">

Carballo, Alba, Marı́a Durbán, and Dae-Jin Lee. 2021. “Out-of-Sample
Prediction in Multidimensional p-Spline Models.” *Mathematics* 9 (15):
1761.

</div>

<div id="ref-henderson1924new" class="csl-entry">

Henderson, Robert. 1924. “A New Method of Graduation.” *Transactions of
the Actuarial Society of America* 25: 29–40.

</div>

<div id="ref-nelder1965optim" class="csl-entry">

Nelder, John A, and Roger Mead. 1965. “A Simplex Method for Function
Minimization.” *The Computer Journal* 7 (4): 308–13.

</div>

<div id="ref-nelder1972glm" class="csl-entry">

Nelder, John Ashworth, and Robert WM Wedderburn. 1972. “Generalized
Linear Models.” *Journal of the Royal Statistical Society: Series A
(General)* 135 (3): 370–84.

</div>

<div id="ref-whittaker1922new" class="csl-entry">

Whittaker, Edmund T. 1922. “On a New Method of Graduation.” *Proceedings
of the Edinburgh Mathematical Society* 41: 63–75.

</div>

</div>
