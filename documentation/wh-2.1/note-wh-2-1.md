---
title: "WH 2.1 : ce qui change par rapport à la 2.0.0"
subtitle: "Huit bugs corrigés, une algèbre bande plus rapide qui garde son empreinte mémoire, la convergence signalée et une suite de tests qui compare chaque ajustement à la version publiée."
surtitre: "WH — note technique"
date: "26 octobre 2026"
pied: "WH 2.0.0.9000 — branche dev — octobre 2026"
---

## La question

Que faut-il changer dans WH 2.0.0, la version publiée sur le CRAN, pour en faire une 2.1 qui corrige ses défauts,
reprenne les avancées numériques d'ebmcd et de WHfull et porte une suite de tests légère ? Cette note présente
chaque modification de la version de développement `2.0.0.9000`, avec sa mesure contre la 2.0.0.

Le cadre fixé le 2026-10-26 : base d'intégration `dev`, `master` restant la 2.0.0 du CRAN jusqu'à la soumission ;
avertissements de convergence toujours émis ; écarts d'arrondi admis jusqu'à 1e-10 relatif ; `lambda_start_ratio`
documenté. Verbatim de la validation : « Oui comme proposé. Tu me confirmes que ta méthodologie garde le stockage
par bandes et l'empreinte mémoire limitée qui en résulte ? Il faudra que tu me présentes la méthode récurrente de
Takahashi ». Toute l'algèbre reste en stockage bande ; aucune matrice dense de la taille du problème n'est formée,
hors `vcov()`, dense par nature et inchangé.

Le découpage avec WHfull (PR #11) est préservé : ni Newton ni Fellner-Schall, la sélection des paramètres de
lissage reste celle de la 2.0.0 (Brent en 1D, Nelder-Mead en 2D). Les gains C++ d'ebmcd (OpenMP, noyaux denses
X'WX, choix du BLAS) n'ont pas d'objet ici : WH n'a pas de matrice de modèle à assembler.

| Domaine | Ce qui change | Effet sur les résultats d'un utilisateur de la 2.0.0 |
|---|---|---|
| Bugs | 8 corrections, chacune avec un test qui échoue sur la 2.0.0 | seulement là où la 2.0.0 plantait ou rendait faux |
| Variance | inverse sélective de Takahashi | écart ≤ 2e-14 à λ fixé |
| Extrapolation 2D | boucle C++ sur les cellules nouvelles | écart ≤ 3e-14 |
| Critère 1D | log-déterminant exact de la pénalité | écart ≤ 1e-12 sur le critère |
| Convergence | champ `convergence`, avertissements | nouveaux avertissements, aucun sur les jeux du paquet |
| Tests | 76 → 380 attentes, 8 fichiers | — |

## Les bugs corrigés

Chaque correction a son bloc dans `tests/testthat/test-bugs.R` ou `test-inputs.R`, rouge sur la 2.0.0 installée et
vert sur la version de développement (constaté le 2026-10-26).

| Bug de la 2.0.0 | Symptôme | Correction |
|---|---|---|
| Cadre régression : `wt * y` avec `y` manquant sur un poids nul | une seule cellule à zéro décès rend tout l'ajustement 1D `NA`, sans erreur ; en 2D la sélection plante (`portfolio_LTC` compte 42 cellules à zéro) | un poids nul annule son observation, même manquante |
| Noms par défaut jamais affectés (1D) ; `seq_along(dim(d)[[1]])` vaut 1 (2D) | `print` affiche des positions vides, `predict` échoue ; `WH()` plante sur toute matrice sans dimnames | noms affectés ; même défaut corrigé dans WHfull |
| `predict()` suppose `newdata` contigu et contenant les observations | `newdata = 100:110` recycle en silence ; dans `c(40:99, 105)` la valeur rendue en 105 est celle de 100 | calcul sur la grille régulière qui couvre observations et `newdata`, rendu sur leur réunion |
| Départ de λ en 2D par `uniroot` sans extension | plantage dès que les poids s'écartent de l'échelle unité | intervalle étendu seulement si la racine n'y est pas (départs ordinaires inchangés), repris de WHfull |
| Évaluation dont le Cholesky échoue | à poids × 1e-8, Nelder-Mead explore un λ qui rend le système numériquement singulier : toute la sélection s'arrête | l'évaluation vaut `Inf`, l'optimiseur recule ; le nombre de ces évaluations est publié |
| Log-déterminant 1D par valeurs propres | q = 3 : faux de 3e-3 à 1 000 points, `NaN` à 3 000 | forme exacte (section « Le log-déterminant exact en 1D ») |
| Entrées invalides | `NA`, négatifs, exposition nulle avec décès : « valeur manquante là où TRUE / FALSE est requis » ; un tableau 3D pris pour une matrice | refus explicite qui nomme l'argument |
| `plot.WH_2d` | l'axe horizontal porte la durée mais s'appelle « x » ; la vue des résidus n'a pas de titre | étiquettes conformes au dessin, titre ajouté |

## La récurrence de Takahashi

### Le problème

Les écarts-types des valeurs lissées sont les racines de la diagonale de la variance $V = A^{-1}$, où
$A = W + P_\lambda$ est la matrice du système, symétrique, définie positive et en bande de largeur $q$. WH en
dispose sous la forme de son facteur de Cholesky $R$, triangulaire supérieur de même bande : $A = R' R$ ($R'$ désigne la transposée de $R$). La 2.0.0
calculait $R^{-1}$ colonne par colonne pour en sommer les carrés, en $O(p^2 q)$ opérations pour $p$ valeurs lissées :
le coût croît comme le carré de la taille de la table.

::: definition
L'**inverse sélective** d'une matrice bande est la partie de son inverse qui tient dans la bande de son facteur de
Cholesky. L'inverse elle-même est pleine, mais cette partie se calcule sans jamais former le reste.
:::

### L'idée

De $V = R^{-1} (R^{-1})'$ on tire $R\,V = (R^{-1})'$. Le membre de droite est triangulaire inférieur, de diagonale
$1 / r_{ii}$ : la ligne $i$ de $R\,V$ est donc nulle à droite de la diagonale, et vaut $1 / r_{ii}$ sur elle. Comme
la ligne $i$ de $R$ ne porte que les coefficients $r_{ik}$ pour $i \le k \le i + q$, cela donne, pour $i < j \le i + q$,

$$V_{ij} = -\frac{1}{r_{ii}} \sum_{k=i+1}^{i+q} \, r_{ik} \, V_{kj} \tag{1}$$

et, sur la diagonale,

$$V_{ii} = \frac{1}{r_{ii}^2} - \frac{1}{r_{ii}} \sum_{k=i+1}^{i+q} \, r_{ik} \, V_{ki} \tag{2}$$

Dans (1) et (2), les indices $k$ et $j$ sont tous deux entre $i + 1$ et $i + q$, donc à moins de $q$ l'un de l'autre :
chaque $V_{kj}$ appelé est dans la bande, et il a déjà été calculé si l'on remonte de la dernière ligne à la
première. C'est la récurrence de Takahashi, Fagan et Chen (1973), formalisée pour les matrices bande par Hutchinson
et de Hoog (1985).

### Coût et mémoire

Chaque ligne demande $q$ coefficients de $q$ termes : $O(p\,q^2)$ opérations au lieu de $O(p^2 q)$. Le résultat
occupe le même stockage compact que $R$, $(q + 1) \times p$ nombres : l'empreinte mémoire est celle de la
factorisation. Le noyau `band_inverse_compact_cpp` (`src/compact.cpp`) rend la bande de $V$ ; `diag_V_compact_cpp`,
appelé par chaque ajustement final et chaque évaluation des critères AIC, BIC et GCV, en lit la diagonale.

### En 2D

La récurrence s'applique à toute matrice bande, donc au système 2D : pour une table $n_1 \times n_2$, $p = n_1 n_2$ et
la bande vaut $q_2 \, n_1$, celle que la factorisation de Cholesky utilise déjà. Le gain théorique est de l'ordre de
$p / q = n_2 / q_2$, et la variance coûte désormais le même ordre que la factorisation elle-même.

| Noyau seul, matrice bande simulée | 2.0.0 | 2.1 | Gain |
|---|--:|--:|--:|
| p = 1 800, bande 120 | 0,15 s | 0,05 s | 3 |
| p = 4 000, bande 200 | 1,66 s | 0,25 s | 6,6 |
| p = 4 000, bande 80 | 0,64 s | 0,03 s | 21 |

Exactitude : écart relatif à `chol2inv` dense de 1,2e-14 sur `portfolio_LTC` (contre 2e-15 pour le noyau de la
2.0.0), bande entière vérifiée à 1e-10 dans `test-kernels.R`.

## La variance d'extrapolation en 2D

Pour prolonger une table 2D, WH suit Carballo et al. (2021) : les valeurs nouvelles valent $-P_{22}^{-1} P_{21} \hat y$
et leur variance est

$$V_{\text{nouv}} = A_e \, V \, A_e' + P_{22}^{-1}, \qquad A_e = -P_{22}^{-1} P_{21} \tag{3}$$

où $P_{22}$ et $P_{21}$ sont les blocs de la pénalité sur la grille étendue. La 2.0.0 en calculait la diagonale par
une boucle R sur les $p$ cellules observées. La 2.1 boucle en C++ sur les $m$ cellules nouvelles : pour la cellule
$k$, $x = P_{22}^{-1} e_k$ (deux résolutions bande avec le facteur de $P_{22}$), $w = P_{12}\,x$ (produit par la
pénalité), puis

$$(A_e V A_e')_{kk} = ||(R^{-1})' w||^2 \tag{4}$$

par une résolution bande avec le facteur de l'ajustement. Seuls des vecteurs de la taille de la grille sont
alloués : une formulation matricielle aurait formé une matrice dense $p \times m$, 30 Mo sur la grille ci-dessous,
écartée pour tenir l'empreinte mémoire. Vérifié contre la formule dense (3) à 1e-8 dans `test-kernels.R`.

## Le log-déterminant exact en 1D

Le critère REML (ou LAML) contient le logarithme du produit des valeurs propres non nulles de la pénalité
$P = D' D$, $D$ étant la matrice des différences d'ordre $q$ sur $n$ points. La 2.0.0 le calculait en sommant les
logarithmes des valeurs propres rendues par LAPACK, les $q$ plus petites étant forcées à zéro. Or les plus petites
valeurs propres non nulles, de l'ordre de $(\pi / n)^{2q}$, perdent leur précision relative sur une longue série :
à q = 3, la somme est fausse de 3e-3 à 1 000 points et devient `NaN` à 3 000, une valeur propre sortant négative.

Ce produit vaut $\det(D D')$, qui a une forme close :

$$\log \det(D D') = \sum_{k=0}^{q-1} \, \left[\log \binom{n+k}{2k+1} - \log \binom{2k}{k}\right] \tag{5}$$

[vérifiée, non démontrée]{.pastille .hyp} La formule (5) est vérifiée en arithmétique entière exacte (élimination de
Bareiss) pour q = 1 à 6 et n ≤ q + 9, et à 3,5e-10 près contre le déterminant dense jusqu'à n = 60 ; je n'en ai pas
écrit la démonstration. Le 1D ne calcule plus aucune valeur propre ; le 2D garde les siennes, nécessaires à la somme
de Kronecker des deux pénalités. Écart sur le critère des jeux du paquet : ≤ 1e-12.

## La convergence signalée

Chaque ajustement porte un champ `convergence` : `converged`, `inner` (itérations et statut du PIRLS en
vraisemblance) et, quand λ est sélectionné, `outer` (évaluations du critère, évaluations dont le Cholesky a échoué,
λ au bord ou non, code de retour de l'optimiseur). Un avertissement est émis quel que soit `verbose` dans trois cas :

- un λ sélectionné à moins d'une décade d'une borne de l'intervalle de recherche [1e-6, 1e12] ;
- Nelder-Mead rend un code de retour non nul, que la 2.0.0 ignorait ;
- le PIRLS atteint son nouveau plafond de 200 itérations (la 2.0.0 bouclait sans fin à `tol_dev = 0`).

::: ecart
**Bord de l'intervalle : une décade, pas la tolérance de WHfull.** WHfull tient λ pour collé à la borne à moins de
2,4e-4 en log. Sur un cas dont l'optimum est à l'infini (log-taux exactement linéaire, q = 2), le critère de WH
devient plat à la précision machine bien avant la borne, et Brent s'arrête à λ = 5,2e11, à 0,65 en log de 1e12 : la
tolérance de WHfull ne le voit pas. Au-delà de 1e11, l'ajustement n'est de toute façon plus que le polynôme non
pénalisé. Aucun ajustement des jeux du paquet, de la vignette ni du README ne déclenche d'avertissement.
:::

## Mesures avant / après

### Temps

Grilles 2D simulées (décès de Poisson, log-taux additif), `dev/bench/bench.R`, 2.0.0 et 2.1 lancées en alternance
dans le même job sur la machine partagée, chargée par d'autres calculs : fourchettes sur deux passes.

| Grille, étape | 2.0.0 | 2.1 | Gain |
|---|--:|--:|--:|
| 100 × 40 (p = 4 000), `predict()` vers 110 × 45 | 27,7 à 30,2 s | 3,2 à 4,2 s | 7 à 9 |
| 100 × 40, ajustement final (λ fixé) | 5,6 à 7,8 s | 0,9 à 1,0 s | 5 à 8 |
| 100 × 40, diagonale de la variance | 6,9 à 7,2 s | 0,84 à 0,85 s | 8 |
| 100 × 40, sélection complète (75 évaluations) | 10,8 à 11,4 s | 4,6 à 5,2 s | 2 |
| 60 × 30 (p = 1 800), `predict()` vers 70 × 35 | 3,1 à 3,2 s | 0,79 s | 4 |
| 60 × 30, ajustement final | 0,45 à 0,55 s | 0,17 à 0,19 s | 3 |
| 60 × 30, sélection complète (177 évaluations) | 3,7 à 4,4 s | 3,6 à 3,8 s | 1 |

Le nombre d'évaluations du critère est le même dans les deux versions : la sélection elle-même n'a pas changé. Une
sélection par REML ne calcule pas la variance à chaque évaluation, d'où le gain plus faible sur la sélection complète ;
une sélection par AIC, BIC ou GCV la calcule à chaque fois. Sur ces grilles simulées, le log-taux est exactement
linéaire en âge : la 2.1 signale à juste titre le λ de cette dimension au bord de l'intervalle. Sortie brute :
`dev/bench/timings-2.0.0.9000.txt`.

### Écarts numériques

Les 34 ajustements de la 2.0.0 figés avant toute modification (`tests/testthat/fixtures/ref-2.0.0.rds`) sont
recalculés par la 2.1 (`dev/ref-2.0.0/compare_ref.R`, sortie dans `gaps-2.0.0.9000.txt`).

| Ajustements | λ | `y_hat` | `std_y_hat`, `edf` |
|---|--:|--:|--:|
| λ fixé, 1D et 2D (2 ajustements) | identique | identique | ≤ 2,1e-15 |
| 2D, λ sélectionné par Nelder-Mead (5) | identique, mêmes évaluations | identique | ≤ 5,3e-15 |
| 1D régression, λ sélectionné (12) | ≤ 2,3e-7 | ≤ 1,1e-9 | ≤ 3,8e-8 |
| 1D vraisemblance, λ sélectionné (12) | ≤ 4,3e-5 | ≤ 9,8e-7 | ≤ 1,7e-5 |
| Prédictions : 2D (2), 1D après sélection (1) | identique ; ≤ 8e-10 | ≤ 5,7e-12 | ≤ 2,9e-10 |

::: ecart
**Le seuil de 1e-10 n'est pas tenu pour un λ sélectionné en 1D**, en régression (jusqu'à 2,3e-7 sur λ) et surtout
en vraisemblance. Il l'est à λ fixé et en 2D. En 1D, la
sélection n'est définie qu'à sa propre tolérance près : celle de Brent (1,2e-4 sur log λ) et l'arrêt du PIRLS, dont
le critère hérite un bruit relatif d'environ 1e-8. Le moindre changement d'arrondi, nouveau noyau ou autre BLAS, peut
faire prendre à Brent un autre chemin dans cette indétermination. Calibration (`dev/ref-2.0.0/calib_noise.R`) : la
2.0.0 elle-même, avec des expositions multipliées par 1 + 1e-14, déplace ses sélections 1D en vraisemblance jusqu'à
4,1e-5 sur λ, 2e-7 sur `y_hat` et 3e-6 sur `std_y_hat`. Les écarts de la 2.1 sont du même ordre, jusqu'à cinq fois
plus sur `y_hat`. Aucune version ne peut tenir 1e-10 sur ces cas.
:::

## La suite de tests

| Fichier | Ce qu'il vérifie | Attentes |
|---|---|--:|
| `test-1D.R`, `test-2D.R` | noyaux bande contre l'algèbre dense (suite de la 2.0.0, inchangée) | 76 |
| `test-kernels.R` | Takahashi, variance d'extrapolation 2D, critère REML contre leurs formules denses | 13 |
| `test-bugs.R` | un bloc par bug, rouge sur la 2.0.0 | 84 |
| `test-inputs.R` | refus des entrées invalides | 16 |
| `test-convergence.R` | champ `convergence` et avertissements | 20 |
| `test-methods.R` | `print`, chaque vue de `plot`, `output_to_df`, `vcov` | 26 |
| `test-reference.R` | les 34 ajustements de la 2.0.0 | 145 |
| **Total** | | **380** |

Durée : 22 s pour l'ensemble des fichiers sur la machine chargée, 31 s de bout en bout par le gate
`dev/run_tests_staged.R` (compilation optimisée, puis 75 s quand le C++ est à recompiler). Le gate est rouge dès
qu'un test échoue, plante ou est sauté.

Face aux BLAS : les échecs de WH au contrôle MKL du CRAN en 2023 et 2024 venaient de tests qui comparaient un calcul
à lui-même avec une tolérance serrée (chantier « Erreurs MKL » d'ebmcd). Aucun test de la 2.1 ne compare deux
résultats flottants à l'identique : les noyaux se comparent à l'algèbre dense à 1e-8 ou 1e-10, et les λ sélectionnés
à l'étalon dans l'indétermination calibrée ci-dessus. `R CMD check` : 0 erreur, 0 avertissement, 1 note (le fichier `.git` du worktree, exclu depuis par le `.Rbuildignore`), vignette non construite ; son code et celui du README tournent sans aucun avertissement.

## Ce qui reste ouvert

- **Avant la soumission** : rejouer les contrôles additionnels du CRAN (MKL, OpenBLAS, sanitizers, valgrind) par le
  skill `cran-tests`, que tu tapes toi-même, puis la promotion de `dev` vers `master` avec ta validation.
- `output_to_df()` calcule les colonnes `wt`, `ec` et `d` puis ne les rend pas. Les rendre changerait la sortie :
  laissé tel quel.
- L'ajustement final d'une sélection 2D tourne dans l'orientation d'origine, à la bande $q_2 n_1$, même quand la
  sélection a tourné à la bande plus étroite de la table transposée. Gain résiduel estimé à 0,3 s sur 100×40 après
  Takahashi, pour un changement de l'orientation de `object$R` lue par `vcov()` et `predict()` : non fait.

## Ce qui reste à décider

::: question
**Quel seuil pour un λ sélectionné ?** Proposition : garder 1e-10 à λ fixé (tenu : ≤ 2e-14), et juger un λ
sélectionné contre la sensibilité propre de la 2.0.0, mesurée par `calib_noise.R`, comme le fait aujourd'hui
`test-reference.R` (tolérance 1e-3 sur λ, 1e-5 sur `y_hat`, 1e-4 sur `std_y_hat` et `edf`).
:::
