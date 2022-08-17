
<!-- README.md is generated from README.Rmd. Please edit that file -->

# WH

<!-- badges: start -->
<!-- badges: end -->

## Notations

Dans cette note, des caractères gras seront utilisés pour désigner les
vecteurs et des lettres majuscules pour les matrices. Si
![\mathbf{y}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7By%7D "\mathbf{y}")
est un vecteur et
![A](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;A "A")
une matrice, on désignera par
![\text{Var}(\mathbf{y})](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctext%7BVar%7D%28%5Cmathbf%7By%7D%29 "\text{Var}(\mathbf{y})")
la matrice de variance-covariance associée à
![\mathbf{y}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7By%7D "\mathbf{y}"),
![\textbf{diag}(A)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7Bdiag%7D%28A%29 "\textbf{diag}(A)")
la diagonale de la matrice
![A](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;A "A")
et
![\text{Diag}(\mathbf{y})](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctext%7BDiag%7D%28%5Cmathbf%7By%7D%29 "\text{Diag}(\mathbf{y})")
la matrice diagonale telle que
![\textbf{diag}(\text{Diag}(\mathbf{y})) = y](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7Bdiag%7D%28%5Ctext%7BDiag%7D%28%5Cmathbf%7By%7D%29%29%20%3D%20y "\textbf{diag}(\text{Diag}(\mathbf{y})) = y"),
par
![\text{tr}(A)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctext%7Btr%7D%28A%29 "\text{tr}(A)")
la somme des valeurs diagonales de
![A](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;A "A"),
par
![A^{T}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;A%5E%7BT%7D "A^{T}")
sa transposée et dans le cas où
![A](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;A "A")
est inversible, par
![A^{- 1}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;A%5E%7B-%201%7D "A^{- 1}")
son inverse, par
![A^{- T}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;A%5E%7B-%20T%7D "A^{- T}")
l’inverse de sa transposée et par
![\|A\|](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%7CA%7C "|A|")
le produit des valeurs propres de
![A](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;A "A").
Le produit de Kronecker de deux matrices
![A](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;A "A")
et
![B](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;B "B")
sera noté
![A \otimes B](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;A%20%5Cotimes%20B "A \otimes B")
et
![A \odot B](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;A%20%5Codot%20B "A \odot B")
désignera le produit de Hadamard *i.e.* le produit élément par élément.
Enfin,
![\text{vec}(A)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctext%7Bvec%7D%28A%29 "\text{vec}(A)")
désignera le vecteur construit en mettant bout à bout les colonnes de
![A](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;A "A").

## Origine du lissage

Le lissage de Whittaker-Henderson est une méthode de régularisation qui
corrige les observations obtenues pour tenir compte des fluctuations
d’échantillonnage. Proposé initialement par Whittaker (1922) pour la
construction des tables de mortalité et enrichi par les travaux de
Henderson (1924), il demeure à ce jour l’une des méthodes les plus
populaires parmi les actuaires pour la construction de tables
d’expérience pour les risques d’assurance de personne : décès, arrêt de
travail, dépendance, perte d’emploi. Le lissage de Whittaker-Henderson
se généralise en effet aux tables en deux dimensions mais nécessite dans
tous les cas de disposer d’observations discrètes régulièrement
espacées.

## Lissage dans le cas unidimensionnel

Soit
![y](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;y "y")
un vecteur d’observations et
![w \ge 0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;w%20%5Cge%200 "w \ge 0")
un vecteur de poids positifs ou nuls, tous deux de taille
![n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n "n").
L’estimateur associé au lissage de Whittaker-Henderson s’écrit :

![\hat{y} = \underset{\theta}{\text{argmin}}\\{F(y,w,\theta) + R\_{\lambda,q}(\theta)\\}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Chat%7By%7D%20%3D%20%5Cunderset%7B%5Ctheta%7D%7B%5Ctext%7Bargmin%7D%7D%5C%7BF%28y%2Cw%2C%5Ctheta%29%20%2B%20R_%7B%5Clambda%2Cq%7D%28%5Ctheta%29%5C%7D "\hat{y} = \underset{\theta}{\text{argmin}}\{F(y,w,\theta) + R_{\lambda,q}(\theta)\}")

où :

-   ![F(y,w,\theta) = \underset{i = 1}{\overset{n}{\sum}} w_i(y_i - \theta_i)^2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;F%28y%2Cw%2C%5Ctheta%29%20%3D%20%5Cunderset%7Bi%20%3D%201%7D%7B%5Coverset%7Bn%7D%7B%5Csum%7D%7D%20w_i%28y_i%20-%20%5Ctheta_i%29%5E2 "F(y,w,\theta) = \underset{i = 1}{\overset{n}{\sum}} w_i(y_i - \theta_i)^2")
    représente un critère de fidélité aux observations

-   ![R(\theta,\lambda,q) = \lambda \underset{i = 1}{\overset{n - q}{\sum}} (\Delta^q\theta)\_i^2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;R%28%5Ctheta%2C%5Clambda%2Cq%29%20%3D%20%5Clambda%20%5Cunderset%7Bi%20%3D%201%7D%7B%5Coverset%7Bn%20-%20q%7D%7B%5Csum%7D%7D%20%28%5CDelta%5Eq%5Ctheta%29_i%5E2 "R(\theta,\lambda,q) = \lambda \underset{i = 1}{\overset{n - q}{\sum}} (\Delta^q\theta)_i^2")
    est un critère de régularité

![\Delta^q](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5CDelta%5Eq "\Delta^q")
représente dans cette dernière expression l’opérateur de différence
avant d’ordre
![q](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;q "q")
tel que pour tout
![i\in\[1,n - q\]](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;i%5Cin%5B1%2Cn%20-%20q%5D "i\in[1,n - q]")
:

![(\Delta^q\theta)\_i = \underset{k = 0}{\overset{q}{\sum}} \begin{pmatrix}q \\\\ k\end{pmatrix}(- 1)^{q - k} \theta\_{i + k}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%28%5CDelta%5Eq%5Ctheta%29_i%20%3D%20%5Cunderset%7Bk%20%3D%200%7D%7B%5Coverset%7Bq%7D%7B%5Csum%7D%7D%20%5Cbegin%7Bpmatrix%7Dq%20%5C%5C%20k%5Cend%7Bpmatrix%7D%28-%201%29%5E%7Bq%20-%20k%7D%20%5Ctheta_%7Bi%20%2B%20k%7D "(\Delta^q\theta)_i = \underset{k = 0}{\overset{q}{\sum}} \begin{pmatrix}q \\ k\end{pmatrix}(- 1)^{q - k} \theta_{i + k}")

Définissons
![W = \text{Diag}(w)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;W%20%3D%20%5Ctext%7BDiag%7D%28w%29 "W = \text{Diag}(w)")
la matrice diagonale des poids et
![D\_{n,q}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;D_%7Bn%2Cq%7D "D_{n,q}")
la matrice des différences d’ordre
![q](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;q "q"),
de dimensions
![(n - q,n)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%28n%20-%20q%2Cn%29 "(n - q,n)"),
telle que
![(D\_{n,q}\theta)\_i = (\Delta^q\theta)\_i](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%28D_%7Bn%2Cq%7D%5Ctheta%29_i%20%3D%20%28%5CDelta%5Eq%5Ctheta%29_i "(D_{n,q}\theta)_i = (\Delta^q\theta)_i").
En pratique, on se limitera aux différences d’ordre
![q\in\\{1,2\\}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;q%5Cin%5C%7B1%2C2%5C%7D "q\in\{1,2\}")
représentées ci-dessous :

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

Les critère de régularité et de fidélité se réécrivent sous forme
matricielle :

![\begin{aligned}
F(y,w,\theta) &= \underset{i = 1}{\overset{n}{\sum}} w_i(y_i - \theta_i)^2 = \\\|\sqrt{W}(y - \theta)\\\|^2 = (y - \theta)^TW(y - \theta) \\\\
R(\theta,\lambda,q) &= \lambda \underset{i = 1}{\overset{n - q}{\sum}} (\Delta^q\theta)\_i^2 = \lambda\\\|D\_{n,q}\theta\\\|^2 = \lambda\theta^TD\_{n,q}^TD\_{n,q}\theta
\end{aligned}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cbegin%7Baligned%7D%0AF%28y%2Cw%2C%5Ctheta%29%20%26%3D%20%5Cunderset%7Bi%20%3D%201%7D%7B%5Coverset%7Bn%7D%7B%5Csum%7D%7D%20w_i%28y_i%20-%20%5Ctheta_i%29%5E2%20%3D%20%5C%7C%5Csqrt%7BW%7D%28y%20-%20%5Ctheta%29%5C%7C%5E2%20%3D%20%28y%20-%20%5Ctheta%29%5ETW%28y%20-%20%5Ctheta%29%20%5C%5C%0AR%28%5Ctheta%2C%5Clambda%2Cq%29%20%26%3D%20%5Clambda%20%5Cunderset%7Bi%20%3D%201%7D%7B%5Coverset%7Bn%20-%20q%7D%7B%5Csum%7D%7D%20%28%5CDelta%5Eq%5Ctheta%29_i%5E2%20%3D%20%5Clambda%5C%7CD_%7Bn%2Cq%7D%5Ctheta%5C%7C%5E2%20%3D%20%5Clambda%5Ctheta%5ETD_%7Bn%2Cq%7D%5ETD_%7Bn%2Cq%7D%5Ctheta%0A%5Cend%7Baligned%7D "\begin{aligned}
F(y,w,\theta) &= \underset{i = 1}{\overset{n}{\sum}} w_i(y_i - \theta_i)^2 = \|\sqrt{W}(y - \theta)\|^2 = (y - \theta)^TW(y - \theta) \\
R(\theta,\lambda,q) &= \lambda \underset{i = 1}{\overset{n - q}{\sum}} (\Delta^q\theta)_i^2 = \lambda\|D_{n,q}\theta\|^2 = \lambda\theta^TD_{n,q}^TD_{n,q}\theta
\end{aligned}")

et l’estimateur associé au lissage devient :

![\hat{y} = \underset{\theta}{\text{argmin}} \left\lbrace(y - \theta)^TW(y - \theta) + \theta^TP\_\lambda\theta\right\rbrace](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Chat%7By%7D%20%3D%20%5Cunderset%7B%5Ctheta%7D%7B%5Ctext%7Bargmin%7D%7D%20%5Cleft%5Clbrace%28y%20-%20%5Ctheta%29%5ETW%28y%20-%20%5Ctheta%29%20%2B%20%5Ctheta%5ETP_%5Clambda%5Ctheta%5Cright%5Crbrace "\hat{y} = \underset{\theta}{\text{argmin}} \left\lbrace(y - \theta)^TW(y - \theta) + \theta^TP_\lambda\theta\right\rbrace")

en notant dans le cas unidimensionnel
![P\_\lambda = \lambda D\_{n,q}^TD\_{n,q}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;P_%5Clambda%20%3D%20%5Clambda%20D_%7Bn%2Cq%7D%5ETD_%7Bn%2Cq%7D "P_\lambda = \lambda D_{n,q}^TD_{n,q}").

## Lissage dans le cas bidimensionnel

Dans le cas bidimensionnel, l’on considère une matrice
![Y](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;Y "Y")
d’observations et une matrice
![\Omega](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5COmega "\Omega")
de poids positifs ou nuls, toutes deux de dimensions
![n_x \times n_z](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n_x%20%5Ctimes%20n_z "n_x \times n_z").

L’estimateur associé au lissage de Whittaker-Henderson s’écrit :

![\widehat{Y} = \underset{\Theta}{\text{argmin}}\\{F(Y,\Omega, \Theta) + R\_{\lambda,q}(\Theta)\\}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cwidehat%7BY%7D%20%3D%20%5Cunderset%7B%5CTheta%7D%7B%5Ctext%7Bargmin%7D%7D%5C%7BF%28Y%2C%5COmega%2C%20%5CTheta%29%20%2B%20R_%7B%5Clambda%2Cq%7D%28%5CTheta%29%5C%7D "\widehat{Y} = \underset{\Theta}{\text{argmin}}\{F(Y,\Omega, \Theta) + R_{\lambda,q}(\Theta)\}")

où

-   ![F(Y,\Omega, \Theta) = \underset{i = 1}{\overset{n_x}{\sum}}\underset{j = 1}{\overset{n_z}{\sum}} \Omega\_{i,j}(Y\_{i,j} - \Theta\_{i,j})^2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;F%28Y%2C%5COmega%2C%20%5CTheta%29%20%3D%20%5Cunderset%7Bi%20%3D%201%7D%7B%5Coverset%7Bn_x%7D%7B%5Csum%7D%7D%5Cunderset%7Bj%20%3D%201%7D%7B%5Coverset%7Bn_z%7D%7B%5Csum%7D%7D%20%5COmega_%7Bi%2Cj%7D%28Y_%7Bi%2Cj%7D%20-%20%5CTheta_%7Bi%2Cj%7D%29%5E2 "F(Y,\Omega, \Theta) = \underset{i = 1}{\overset{n_x}{\sum}}\underset{j = 1}{\overset{n_z}{\sum}} \Omega_{i,j}(Y_{i,j} - \Theta_{i,j})^2")
    représente un critère de fidélité aux observations

-   ![R(\Theta,\lambda,q) = \lambda_x \underset{j = 1}{\overset{n_z}{\sum}}\underset{i = 1}{\overset{n_x - q_x}{\sum}} (\Delta^{q_x}\Theta\_{\bullet,j})\_i^2 + \lambda_z \underset{i = 1}{\overset{n_x}{\sum}}\underset{j = 1}{\overset{n_z - q_z}{\sum}} (\Delta^{q_z}\Theta\_{i,\bullet})\_j^2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;R%28%5CTheta%2C%5Clambda%2Cq%29%20%3D%20%5Clambda_x%20%5Cunderset%7Bj%20%3D%201%7D%7B%5Coverset%7Bn_z%7D%7B%5Csum%7D%7D%5Cunderset%7Bi%20%3D%201%7D%7B%5Coverset%7Bn_x%20-%20q_x%7D%7B%5Csum%7D%7D%20%28%5CDelta%5E%7Bq_x%7D%5CTheta_%7B%5Cbullet%2Cj%7D%29_i%5E2%20%2B%20%5Clambda_z%20%5Cunderset%7Bi%20%3D%201%7D%7B%5Coverset%7Bn_x%7D%7B%5Csum%7D%7D%5Cunderset%7Bj%20%3D%201%7D%7B%5Coverset%7Bn_z%20-%20q_z%7D%7B%5Csum%7D%7D%20%28%5CDelta%5E%7Bq_z%7D%5CTheta_%7Bi%2C%5Cbullet%7D%29_j%5E2 "R(\Theta,\lambda,q) = \lambda_x \underset{j = 1}{\overset{n_z}{\sum}}\underset{i = 1}{\overset{n_x - q_x}{\sum}} (\Delta^{q_x}\Theta_{\bullet,j})_i^2 + \lambda_z \underset{i = 1}{\overset{n_x}{\sum}}\underset{j = 1}{\overset{n_z - q_z}{\sum}} (\Delta^{q_z}\Theta_{i,\bullet})_j^2")
    est un critère de régularité bidimensionnel s’écrivant comme la
    somme :

    -   d’un critère de régularité unidimensionnel appliqué à toutes les
        lignes de
        ![\Theta](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5CTheta "\Theta")
        et

    -   d’un critère de régularité unidimensionnel appliqué à toutes les
        colonnes de
        ![\Theta](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5CTheta "\Theta").

Là encore il est souhaitable d’adopter une forme matricielle en
définissant
![y = \text{vec}(Y)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;y%20%3D%20%5Ctext%7Bvec%7D%28Y%29 "y = \text{vec}(Y)"),
![w = \text{vec}(\Omega)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;w%20%3D%20%5Ctext%7Bvec%7D%28%5COmega%29 "w = \text{vec}(\Omega)"),
![\theta = \text{vec}(\Theta)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctheta%20%3D%20%5Ctext%7Bvec%7D%28%5CTheta%29 "\theta = \text{vec}(\Theta)")
les vecteurs obtenus en mettant bout à bout les colonnes des matrices
![Y](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;Y "Y"),
![\Omega](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5COmega "\Omega")
et
![\Theta](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5CTheta "\Theta")
respectivement et en notant
![W = \text{Diag}(w)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;W%20%3D%20%5Ctext%7BDiag%7D%28w%29 "W = \text{Diag}(w)")
et
![n = n_x \times n_z](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n%20%3D%20n_x%20%5Ctimes%20n_z "n = n_x \times n_z").

Les critère de régularité et de fidélité se réécrivent sous forme
matricielle, en fonction de
![y](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;y "y"),
![w](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;w "w")
et
![\theta](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctheta "\theta")
:

![\begin{aligned}
F(y,w, \theta) &= \underset{i = 1}{\overset{n}{\sum}}w_i(y_i - \theta_i)^2 = (y - \theta)^TW(y - \theta) \\\\
R(\theta,\lambda,q) &= \theta^{T}(\lambda_x I\_{n_z} \otimes D\_{n_x,q_x}^{T}D\_{n_x,q_x} + \lambda_z D\_{n_z,q_z}^{T}D\_{n_z,q_z} \otimes I\_{n_x}) \theta
\end{aligned}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cbegin%7Baligned%7D%0AF%28y%2Cw%2C%20%5Ctheta%29%20%26%3D%20%5Cunderset%7Bi%20%3D%201%7D%7B%5Coverset%7Bn%7D%7B%5Csum%7D%7Dw_i%28y_i%20-%20%5Ctheta_i%29%5E2%20%3D%20%28y%20-%20%5Ctheta%29%5ETW%28y%20-%20%5Ctheta%29%20%5C%5C%0AR%28%5Ctheta%2C%5Clambda%2Cq%29%20%26%3D%20%5Ctheta%5E%7BT%7D%28%5Clambda_x%20I_%7Bn_z%7D%20%5Cotimes%20D_%7Bn_x%2Cq_x%7D%5E%7BT%7DD_%7Bn_x%2Cq_x%7D%20%2B%20%5Clambda_z%20D_%7Bn_z%2Cq_z%7D%5E%7BT%7DD_%7Bn_z%2Cq_z%7D%20%5Cotimes%20I_%7Bn_x%7D%29%20%5Ctheta%0A%5Cend%7Baligned%7D "\begin{aligned}
F(y,w, \theta) &= \underset{i = 1}{\overset{n}{\sum}}w_i(y_i - \theta_i)^2 = (y - \theta)^TW(y - \theta) \\
R(\theta,\lambda,q) &= \theta^{T}(\lambda_x I_{n_z} \otimes D_{n_x,q_x}^{T}D_{n_x,q_x} + \lambda_z D_{n_z,q_z}^{T}D_{n_z,q_z} \otimes I_{n_x}) \theta
\end{aligned}")

et l’on retrouve la même forme que pour le cas unidimensionnel avec dans
le cas bidimensionnel
![P\_\lambda = \lambda_x I\_{n_z} \otimes D\_{n_x,q_x}^{T}D\_{n_x,q_x} + \lambda_z D\_{n_z,q_z}^{T}D\_{n_z,q_z} \otimes I\_{n_x}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;P_%5Clambda%20%3D%20%5Clambda_x%20I_%7Bn_z%7D%20%5Cotimes%20D_%7Bn_x%2Cq_x%7D%5E%7BT%7DD_%7Bn_x%2Cq_x%7D%20%2B%20%5Clambda_z%20D_%7Bn_z%2Cq_z%7D%5E%7BT%7DD_%7Bn_z%2Cq_z%7D%20%5Cotimes%20I_%7Bn_x%7D "P_\lambda = \lambda_x I_{n_z} \otimes D_{n_x,q_x}^{T}D_{n_x,q_x} + \lambda_z D_{n_z,q_z}^{T}D_{n_z,q_z} \otimes I_{n_x}").

## Installation

You can install the development version of WH from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("GuillaumeBiessy/WH")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(WH)
## basic example code
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this. You could also
use GitHub Actions to re-render `README.Rmd` every time you push. An
example workflow can be found here:
<https://github.com/r-lib/actions/tree/v1/examples>.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-henderson1924new" class="csl-entry">

Henderson, Robert. 1924. “A New Method of Graduation.” *Transactions of
the Actuarial Society of America* 25: 29–40.

</div>

<div id="ref-whittaker1922new" class="csl-entry">

Whittaker, Edmund T. 1922. “On a New Method of Graduation.” *Proceedings
of the Edinburgh Mathematical Society* 41: 63–75.

</div>

</div>
