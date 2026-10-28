# WH — doctrine du dépôt

Paquet R de lissage de Whittaker-Henderson en 1 et 2 dimensions, pour la graduation de tables
actuarielles utilisées en assurance vie. Méthodes du papier Biessy (2025), arXiv:2306.06932.
**Publié sur le CRAN** : version 2.0.0, soumise le 2025-06-19.

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
  `compute_res_deviance`), initialisation de λ (`find_starting_lambda`), et les passerelles vers
  la représentation bande compacte (`combine_P_compact`, `get_y_new_compact`,
  `get_diag_V_pred_compact`, `banded_to_compact`, `compact_to_sym`).
- `src/` : noyaux C++ — `compact.cpp` (Cholesky et résolutions en stockage bande compact),
  `lapack.cpp` (appels LAPACK bande) ; enveloppes générées dans `R/RcppExports.R`.
- `data/` : jeux `portfolio_mort` (1D) et `portfolio_LTC` (2D). `vignettes/WH.Rmd` ;
  `README.Rmd`, rendu en `README.md`.
- `tests/` : `testthat` édition 3 — `test-1D.R` et `test-2D.R`, qui confrontent le chemin bande à
  l'algèbre dense équivalente écrite en base R dans le test même ; `spelling.R`, sauté sur le CRAN.
- `.github/workflows/` : `R-CMD-check.yaml` (5 plateformes, sur `master`) et `rhub.yaml` (manuel).
- `.doc-projet` : pointeur vers la vitrine HTML du package, qui vit hors du dépôt.

## Commandes

Le dépôt n'a **aucun script de gate** : `dev/run_tests_staged.R` n'existe pas, et c'est un fait,
pas un oubli — WH est figé sur sa version CRAN, sa suite ne sert qu'à un cycle de développement
délibéré. Deux conséquences mécaniques :

- `land-to-master.sh` **refuse** tout atterrissage sans `--test-cmd` explicite, le repli
  `--no-test` étant lui-même refusé ici (§ Intégration) : un futur cycle dev commence donc par
  se donner son gate ;
- aucune commande de suite n'est figée ici. L'écrire reviendrait à choisir ce gate sans l'avoir
  exercé. La suite est une suite `testthat` ordinaire, sur le paquet compilé et chargé, à lancer
  sous le verrou machine du skill `gate-lock-r` — et ⛔ jamais par `Rscript -e`, qui échoue sous
  Git Bash : passer par un fichier `.R`.

## Intégration

Le dépôt n'a **pas de base d'intégration déclarée** — `dev/land-base.txt` n'existe pas — et la
branche par défaut que GitHub déclare est **`master`**. Ici, donc, atterrir et promouvoir sont
**le même geste**, et `master` est la branche publiée au CRAN.

- **Naissance** : un worktree créé par l'application naît sur `origin/master`, la branche par
  défaut que **GitHub** déclare, ⛔ ni sur la référence locale `origin/HEAD`, ⛔ ni sur le HEAD du
  checkout principal. Il peut donc être en retard. Premier geste dans un worktree neuf et propre :
  `git fetch origin master && git merge --ff-only origin/master` ; travail déjà committé :
  `git merge origin/master`. ⛔ Pas de `git reset --hard`. Constat obligatoire :
  `git rev-list --count HEAD..origin/master` = 0.
- **Atterrissage — la chaîne `land.sh` ne passe pas ici, et c'est mécanique** (mesuré le
  2026-09-22) : `land-to-master.sh` exige une commande de tests **avant même** de regarder si le
  diff est documentaire, et le repli `--no-test` est refusé vers une branche protégée, ce
  qu'est `master` faute de base déclarée. Sans script de gate dans le dépôt, tout atterrissage
  par la chaîne échoue en précondition. Donc :
  - **diff 100 % documentaire** (`.md`, `.html`, `.gitignore`) → `git merge origin/master` puis
    `git push origin HEAD:master`. Le garde du verrou machine exempte ces diffs par construction,
    et aucune suite de tests ne les lit : c'est la voie normale, pas un contournement.
  - **tout le reste** → se donner d'abord un gate, le passer en `--test-cmd` à `land.sh`, et
    ⛔ jamais un `--base` tapé à la main.
- ⚠️ **Tout ce que `R CMD build` empaquette est livrable** : `R/`, `src/`, `man/`, `tests/`,
  `data/`, `vignettes/`, `NAMESPACE`, `DESCRIPTION`, `NEWS.md`, `README.md`. Y toucher ouvre un
  cycle CRAN — montée de version et resoumission — et exige la **validation express de
  l'utilisateur** : ⛔ ne rien lander de tel de sa propre initiative. Le piège est que `NEWS.md` et
  `README.md` sont des `.md` : la voie doc-only du lander les pousserait sur `master` **sans
  gate**, faisant diverger en silence la branche publiée du tarball soumis. Seuls les fichiers
  listés au `.Rbuildignore` sont hors livrable et relèvent du doc-only.
- Ce `CLAUDE.md` n'est **pas** inscrit au `.Rbuildignore` : `R CMD check` le signale donc comme
  fichier non standard à la racine du paquet — une note, pas une erreur. L'y inscrire touche un
  fichier non documentaire et relève du même cycle dev, seul en mesure de payer son gate.

## Décisions de l'auteur

- **2026-09-22**, verbatim : « Après réflexion, la règle, ici ainsi que dans tous les dépôts qui
  ont une suite de test (WH, WHfull et ebmcd) doit exister dans le Claude.md. Interdiction de la
  déléguer dans la fiche tech, et la mémoire court-terme n'est pas non plus le bon endroit après
  mûre réflexion. » La règle visée est le § Intégration ci-dessus. Cet arbitrage **révoque** celui
  du 2026-09-19, qui la logeait hors de ce fichier.
- **2026-06-04** (PR #11) — la partition WH / WHfull est un **choix de périmètre**, pas un manque :
  la confrontation du code au papier n'a laissé aucune divergence dure. Ce qui en découle, et qui
  est l'état du dépôt : WH reste figé sur sa version CRAN, et toute évolution de son source
  livrable, documentation roxygen comprise, relève d'un cycle dev délibéré.
