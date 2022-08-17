
<!-- README.md is generated from README.Rmd. Please edit that file -->

# WH

<!-- badges: start -->
<!-- badges: end -->

## Notations

Dans cette note, des caractères gras seront utilisés pour désigner les
vecteurs et des lettres majuscules pour les matrices. Si $\mathbf{y}$
est un vecteur et $A$ une matrice, on désignera par
$\text{Var}(\mathbf{y})$ la matrice de variance-covariance associée à
$\mathbf{y}$, $\textbf{diag}(A)$ la diagonale de la matrice $A$ et
$\text{Diag}(\mathbf{y})$ la matrice diagonale telle que
$\textbf{diag}(\text{Diag}(\mathbf{y})) = y$, par $\text{tr}(A)$ la
somme des valeurs diagonales de $A$, par $A^{T}$ sa transposée et dans
le cas où $A$ est inversible, par $A^{- 1}$ son inverse, par $A^{- T}$
l’inverse de sa transposée et par $|A|$ le produit des valeurs propres
de $A$. Le produit de Kronecker de deux matrices $A$ et $B$ sera noté
$A \otimes B$ et $A \odot B$ désignera le produit de Hadamard *i.e.* le
produit élément par élément. Enfin, $\text{vec}(A)$ désignera le vecteur
construit en mettant bout à bout les colonnes de $A$.

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

Soit $y$ un vecteur d’observations et $w \ge 0$ un vecteur de poids
positifs ou nuls, tous deux de taille $n$. L’estimateur associé au
lissage de Whittaker-Henderson s’écrit :

$$\hat{y} = \underset{\theta}{\text{argmin}}\{F(y,w,\theta) + R_{\lambda,q}(\theta)\}$$

où :

-   $F(y,w,\theta) = \underset{i = 1}{\overset{n}{\sum}} w_i(y_i - \theta_i)^2$
    représente un critère de fidélité aux observations

-   $R(\theta,\lambda,q) = \lambda \underset{i = 1}{\overset{n - q}{\sum}} (\Delta^q\theta)_i^2$
    est un critère de régularité

$\Delta^q$ représente dans cette dernière expression l’opérateur de
différence avant d’ordre $q$ tel que pour tout $i\in[1,n - q]$ :

$$(\Delta^q\theta)_i = \underset{k = 0}{\overset{q}{\sum}} \begin{pmatrix}q \\ k\end{pmatrix}(- 1)^{q - k} \theta_{i + k}.$$

Définissons $W = \text{Diag}(w)$ la matrice diagonale des poids et
$D_{n,q}$ la matrice des différences d’ordre $q$, de dimensions
$(n - q,n)$, telle que $(D_{n,q}\theta)_i = (\Delta^q\theta)_i$. En
pratique, on se limitera aux différences d’ordre $q\in\{1,2\}$
représentées ci-dessous :

$$
D_1 = \begin{pmatrix}
1 & - 1 &  & 0 \\
& \ddots & \ddots & \\
0 & & 1 & - 1
\end{pmatrix}
\quad\quad
D_2 = \begin{pmatrix}
1 & - 2 & 1 & & 0 \\
& \ddots & \ddots & \ddots & \\
0 & & 1 & - 2 & 1
\end{pmatrix}.
$$

Les critère de régularité et de fidélité se réécrivent sous forme
matricielle :

$$
\begin{aligned}
F(y,w,\theta) &= \underset{i = 1}{\overset{n}{\sum}} w_i(y_i - \theta_i)^2 = \Vert \sqrt{W}(y - \theta)\Vert^2 = (y - \theta)^TW(y - \theta) \\
R(\theta,\lambda,q) &= \lambda \underset{i = 1}{\overset{n - q}{\sum}} (\Delta^q\theta)_i^2 = \lambda\Vert D_{n,q}\theta\Vert^2 = \lambda\theta^TD_{n,q}^TD_{n,q}\theta
\end{aligned}
$$

et l’estimateur associé au lissage devient :

$$\hat{y} = \underset{\theta}{\text{argmin}} \left\lbrace(y - \theta)^TW(y - \theta) + \theta^TP_\lambda\theta\right\rbrace$$

en notant dans le cas unidimensionnel
$P_\lambda = \lambda D_{n,q}^TD_{n,q}$.

## Lissage dans le cas bidimensionnel

Dans le cas bidimensionnel, l’on considère une matrice $Y$
d’observations et une matrice $\Omega$ de poids positifs ou nuls, toutes
deux de dimensions $n_x \times n_z$.

L’estimateur associé au lissage de Whittaker-Henderson s’écrit :

$$\widehat{Y} = \underset{\Theta}{\text{argmin}}\{F(Y,\Omega, \Theta) + R_{\lambda,q}(\Theta)\}$$

où

-   $F(Y,\Omega, \Theta) = \underset{i = 1}{\overset{n_x}{\sum}}\underset{j = 1}{\overset{n_z}{\sum}} \Omega_{i,j}(Y_{i,j} - \Theta_{i,j})^2$
    représente un critère de fidélité aux observations

-   $R(\Theta,\lambda,q) = \lambda_x \underset{j = 1}{\overset{n_z}{\sum}}\underset{i = 1}{\overset{n_x - q_x}{\sum}} (\Delta^{q_x}\Theta_{\bullet,j})_i^2 + \lambda_z \underset{i = 1}{\overset{n_x}{\sum}}\underset{j = 1}{\overset{n_z - q_z}{\sum}} (\Delta^{q_z}\Theta_{i,\bullet})_j^2$
    est un critère de régularité bidimensionnel s’écrivant comme la
    somme :

    -   d’un critère de régularité unidimensionnel appliqué à toutes les
        lignes de $\Theta$ et

    -   d’un critère de régularité unidimensionnel appliqué à toutes les
        colonnes de $\Theta$.

Là encore il est souhaitable d’adopter une forme matricielle en
définissant $y = \text{vec}(Y)$, $w = \text{vec}(\Omega)$,
$\theta = \text{vec}(\Theta)$ les vecteurs obtenus en mettant bout à
bout les colonnes des matrices $Y$, $\Omega$ et $\Theta$ respectivement
et en notant $W = \text{Diag}(w)$ et $n = n_x \times n_z$.

Les critère de régularité et de fidélité se réécrivent sous forme
matricielle, en fonction de $y$, $w$ et $\theta$ :

$$
\begin{aligned}
F(y,w, \theta) &= \underset{i = 1}{\overset{n}{\sum}}w_i(y_i - \theta_i)^2 = (y - \theta)^TW(y - \theta) \\
R(\theta,\lambda,q) &= \theta^{T}(\lambda_x I_{n_z} \otimes D_{n_x,q_x}^{T}D_{n_x,q_x} + \lambda_z D_{n_z,q_z}^{T}D_{n_z,q_z} \otimes I_{n_x}) \theta
\end{aligned}
$$

et l’on retrouve la même forme que pour le cas unidimensionnel avec dans
le cas bidimensionnel
$P_\lambda = \lambda_x I_{n_z} \otimes D_{n_x,q_x}^{T}D_{n_x,q_x} + \lambda_z D_{n_z,q_z}^{T}D_{n_z,q_z} \otimes I_{n_x}$.

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

## Bibliographie

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
