# WH — doctrine du dépôt

Paquet R de lissage de Whittaker-Henderson en 1 et 2 dimensions, pour la graduation de tables
actuarielles utilisées en assurance vie. Méthodes du papier Biessy (2025), arXiv:2306.06932.
**Publié sur le CRAN** : version 2.0.0, soumise le 2025-06-19, portée par `master`. La version de
développement `2.0.0.9000`, future 2.1, vit sur `dev` depuis le 2026-10-26 ; ses modifications sont
présentées dans `documentation/wh-2.1/note-wh-2-1.html` (source `.md` à côté).

WH est le **sous-ensemble léger assumé** de `WHfull` (`C:\Packages R\WHfull`) : sélection des
paramètres de lissage par emboîtement `outer` avec un optimiseur sans dérivée, et algèbre **bande**
seule — ni Fellner-Schall, ni Newton, ni réduction de rang, ni GLAM.

## Carte du code

- `R/main.R` : façade `WH()` (validation, dispatch 1D / 2D), `WH_outer` (sélection de λ par
  `optim` — Brent en 1D, Nelder-Mead en 2D), `WH_fixed_lambda`, les ajusteurs internes
  `fit_pen_reg`, `fit_WH_fixed_lambda_reg` et `fit_WH_fixed_lambda_ml` (PIRLS), les méthodes S3
  `predict` / `vcov` / `print` / `plot` en versions `WH_1d` et `WH_2d`, et `output_to_df`.
- `R/tools.R` : matrices du modèle et opérateurs de différences (`get_model_matrices`, `q_diff`),
  critères de sélection (`compute_REML`, `compute_K`) et diagnostics (`get_diagnosis`,
  `compute_res_deviance`), initialisation de λ (`find_starting_lambda`),
  log-pseudo-déterminant exact de la pénalité 1D (`log_pdet_P`, forme close en coefficients
  binomiaux), grille de prédiction (`get_pred_grid`), constantes nommées (`wh_log_lambda_interval`,
  `wh_edge_tol`, `wh_max_iter`), et les passerelles vers la représentation bande compacte
  (`combine_P_compact`, `get_y_new_compact`, `banded_to_compact`, `compact_to_sym`).
- `src/` : noyaux C++ — `compact.cpp` (résolutions en stockage bande compact, inverse sélective de
  Takahashi `band_inverse_compact_cpp` et sa diagonale `diag_V_compact_cpp`, variance
  d'extrapolation 2D `diag_V_pred_compact_cpp`), `lapack.cpp` (appels LAPACK bande : Cholesky,
  résolutions, valeurs propres) ; enveloppes générées dans `R/RcppExports.R`. Toute l'algèbre reste
  en stockage bande, mémoire O(p·q) : exigence de l'auteur (2026-10-26), seul `vcov()` est dense.
- `data/` : jeux `portfolio_mort` (1D) et `portfolio_LTC` (2D). `vignettes/WH.Rmd` ;
  `README.Rmd`, rendu en `README.md`.
- `tests/` : `testthat` édition 3 — `test-1D.R`, `test-2D.R` et `test-kernels.R` confrontent le
  chemin bande à l'algèbre dense équivalente écrite en base R dans le test même ; `test-bugs.R` (un
  bloc par bug de la 2.0.0, rouge sur elle), `test-inputs.R`, `test-convergence.R`,
  `test-methods.R` ; `test-reference.R` compare 34 ajustements à ceux de la 2.0.0
  (`fixtures/ref-2.0.0.rds`, recalculés par `helper-reference.R`) ; `spelling.R`, sauté sur le
  CRAN. Aucun test ne compare deux flottants à l'identique : c'est ce qui faisait échouer WH sous
  MKL au CRAN en 2023 et 2024.
- `dev/` (hors du tarball) : gate `run_tests_staged.R`, base `land-base.txt`, étalon de la 2.0.0
  (`ref-2.0.0/` : `make_ref.R`, `compare_ref.R`, `calib_noise.R` et leurs sorties), temps
  (`bench/bench.R`).
- `.github/workflows/` : `R-CMD-check.yaml` (5 plateformes, sur `master`) et `rhub.yaml` (manuel).
- `.doc-projet` : pointeur vers la vitrine HTML du package, qui vit hors du dépôt.

## Commandes

Le gate est `dev/run_tests_staged.R`, depuis la racine du paquet :

```bash
"/c/Program Files/R/R-4.5.3/bin/Rscript.exe" dev/run_tests_staged.R > dev/.last_test_run.log 2>&1
```

Il compile le C++ en `-O2` (`pkgbuild::compile_dll(debug = FALSE)` : `load_all` seul compile en
`-O0`, dix fois plus lent sur les tests 2D), charge le paquet par `load_all`, lance chaque fichier
de test et sort en code 1 dès qu'un test échoue, plante ou est sauté. Mesuré le 2026-10-26 :
380 attentes, 18 à 22 s de tests, 31 s de bout en bout, 75 s quand le C++ est à recompiler.
`land.sh` le détecte et le lance de lui-même. Un lancement en tâche de fond passe par le verrou
machine du skill `gate-lock-r` ; ⛔ jamais par `Rscript -e`, qui échoue sous Git Bash.

Une version de travail se mesure installée dans une **bibliothèque privée**
(`R CMD INSTALL --library=<dossier>`), ⛔ jamais par-dessus le WH 2.0.0 de la bibliothèque
utilisateur : `dev/bench/bench.R` et `dev/ref-2.0.0/compare_ref.R` prennent ce dossier en argument.

## Intégration

La base d'intégration est **`dev`**, déclarée par `dev/land-base.txt` (décision de l'auteur du
2026-10-26). `master`, branche par défaut que GitHub déclare, porte **la version publiée au CRAN**
(2.0.0) et ne reçoit le cycle 2.1 qu'à la promotion.

- **Naissance** : un worktree créé par l'application naît sur `origin/master`, la branche par
  défaut que **GitHub** déclare, ⛔ ni sur la référence locale `origin/HEAD`, ⛔ ni sur le HEAD du
  checkout principal — donc **hors de la base**. Premier geste dans un worktree neuf et propre :
  `git fetch origin dev && git merge --ff-only origin/dev` ; travail déjà committé :
  `git merge origin/dev`. ⛔ Pas de `git reset --hard`. Constat obligatoire :
  `git rev-list --count HEAD..origin/dev` = 0. Une branche qui ne porte pas `dev/land-base.txt`
  fait refuser le lander.
- **Atterrissage sur `dev`** : par `land.sh` (skill `gate-lock-r`), qui lit la base dans
  `dev/land-base.txt` et lance le gate ; ⛔ jamais un `--base` tapé à la main, ⛔ jamais un
  `git push` direct vers `dev`, que le garde refuse. Mesuré le 2026-10-26 : 25 commits atterris
  par cette voie, gate vert en 18 s.
- **Promotion `dev` → `master`** : c'est une soumission au CRAN. Elle exige la **validation express
  de l'auteur**, après les contrôles additionnels du CRAN (skill `cran-tests`, tapé par l'auteur) ;
  ⛔ ne rien promouvoir de sa propre initiative.
- ⚠️ **Tout ce que `R CMD build` empaquette est livrable** : `R/`, `src/`, `man/`, `tests/`,
  `data/`, `vignettes/`, `NAMESPACE`, `DESCRIPTION`, `NEWS.md`, `README.md`. Il atterrit sur `dev`,
  jamais directement sur `master`. Le piège est que `NEWS.md` et `README.md` sont des `.md` : une
  voie doc-only vers `master` les y pousserait **sans gate**, faisant diverger en silence la
  branche publiée du tarball soumis. Seuls les fichiers listés au `.Rbuildignore` (dont `dev/`,
  `documentation/` et ce `CLAUDE.md`) sont hors livrable.

## Décisions de l'auteur

- **2026-09-22**, verbatim : « Après réflexion, la règle, ici ainsi que dans tous les dépôts qui
  ont une suite de test (WH, WHfull et ebmcd) doit exister dans le Claude.md. Interdiction de la
  déléguer dans la fiche tech, et la mémoire court-terme n'est pas non plus le bon endroit après
  mûre réflexion. » La règle visée est le § Intégration ci-dessus. Cet arbitrage **révoque** celui
  du 2026-09-19, qui la logeait hors de ce fichier.
- **2026-10-26** — ouverture du cycle 2.1. Demande, verbatim : « crée une session dans le repo WH qui
  va créer une version de développement qui deviendra plus tard la 2.1 pour améliorer le package sur
  plusieurs aspects : aligner les développements numériques sur nos avancées les plus récentes
  (ebmcd, WHfull a normalement déjà partiellement été amélioré sur ce point), corriger les bugs
  existants, mettre en place une suite de tests digne de ce nom (mais légère et rapide). Toutes les
  modifications réalisées devront m'être présentées via une note ». Plan validé le même jour, verbatim :
  « Oui comme proposé. Tu me confirmes que ta méthodologie garde le stockage par bandes et
  l'empreinte mémoire limitée qui en résulte ? » — base `dev`, avertissements de convergence émis
  quel que soit `verbose`, `lambda_start_ratio` documenté, stockage bande et mémoire O(p·q)
  maintenus.
- **2026-06-04** (PR #11) — la partition WH / WHfull est un **choix de périmètre**, pas un manque :
  la confrontation du code au papier n'a laissé aucune divergence dure. Ce qui en découle, et qui
  est l'état du dépôt : WH reste figé sur sa version CRAN, et toute évolution de son source
  livrable, documentation roxygen comprise, relève d'un cycle dev délibéré.
