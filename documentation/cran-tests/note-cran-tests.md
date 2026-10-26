---
title: "Rejouer les contrôles du CRAN avant une soumission"
subtitle: "Les contrôles de BLAS du CRAN (MKL, OpenBLAS) et l'absence des paquets suggérés se rejouent fidèlement sous Windows, sur cette machine ; les sanitizers, valgrind, noLD, les compilateurs de Linux et le Mac exigent Linux ou macOS, et R-hub les fournit déjà configurés, sans rien installer ici."
surtitre: "WH — note technique"
date: "25 octobre 2026"
pied: "WH — octobre 2026"
---

## La question

Le CRAN contrôle chaque paquet sur 13 environnements réguliers (*check flavors*) et sur une trentaine de configurations
additionnelles (*additional issues*) : autres bibliothèques d'algèbre linéaire, sanitizers, valgrind, R sans
`long double`, Mac M1… C'est l'une d'elles, MKL, qui a fait échouer WH en 2023 et en 2024. La note répond à trois
questions : lesquelles de ces configurations peuvent faire échouer un paquet R avec du C++ comme WH, lesquelles se
rejouent sur cette machine (Windows 11, WSL Ubuntu 20.04, i7-1165G7 AVX-512) et à quel prix, et ce que donne WH 2.0.0,
la version sur le CRAN, sous chacune.

## Ce que fait le CRAN

Sources lues le 25/10/2026 : la page des environnements de contrôle (`cran.r-project.org/web/checks/check_flavors.html`),
celle des types de contrôles additionnels (`check_issue_kinds.html`) et les README de B. Ripley
(`stats.ox.ac.uk/pub/bdr/<type>/README.txt`).

### Les environnements réguliers

| Environnement | Système | Compilateur | Processeur |
|---|---|---|---|
| r-devel-linux-x86_64-debian-clang | Debian testing | clang/flang 23.1.0 | 2 × Xeon E5-2690 |
| r-devel-linux-x86_64-debian-gcc | Debian testing | GCC 16.2.0 | 2 × Xeon E5-2690 |
| r-devel-linux-x86_64-fedora-clang | Fedora 44 | clang/flang 23.1.0, libc++ | 2 × Xeon E5-2690 v4 |
| r-devel-linux-x86_64-fedora-gcc | Fedora 44 | GCC 16.2.1 | 2 × Xeon E5-2690 v4 |
| r-patched-linux, r-release-linux | Debian testing | GCC 15.3.0 | 2 × Xeon E5-2690 |
| r-devel-windows, r-release-windows | Windows Server 2022 | GCC 14.3.0 | 2 × AMD EPYC 7443 |
| r-oldrel-windows | Windows Server 2022 | GCC 14.3.0 | 2 × Xeon E5-2680 v4 |
| r-release-macos-arm64 | macOS 15.7.1 | Apple Clang 1700 | Apple M3 Ultra |
| r-oldrel-macos-arm64 | macOS 13.4 | Apple Clang 1400 | Apple M1 |
| r-release / r-oldrel-macos-x86_64 | macOS 13.3.1 | Apple Clang 1403 | Core i5 6 cœurs |

Les deux environnements Debian (gcc et clang) lient WH à la LAPACK du système et au BLAS de référence de R
(`-llapack -lRblas` dans les journaux d'installation de WH 2.0.0) ; fedora-gcc utilise la LAPACK du système (3.12.0),
fedora-clang la LAPACK interne de R (`--without-lapack`), comme le contrôle M1mac ; Windows utilise la LAPACK et le
BLAS internes de R. La configuration des environnements macOS réguliers n'a pas été lue.

### Les contrôles additionnels qui peuvent faire échouer WH

WH est du R et du C++ (Rcpp) qui appelle LAPACK en mode bande (`dpbtrf`, `dtbtrs`, `dsbev`), sans OpenMP, sans
`\donttest`, et dont les tests comparent des résultats numériques à une tolérance. Les types retenus :

| Type | Ce que fait le CRAN | Ce qu'il peut révéler dans WH |
|---|---|---|
| MKL | Fedora 44, R-devel gcc, oneMKL 2026.1.0 **en série** (`-lmkl_gf_lp64 -lmkl_core -lmkl_sequential`), BLAS et LAPACK ; 2023.2.0 avant 2026-07 | résultats qui changent avec les noyaux de MKL : les échecs de 2023 et 2024 |
| OpenBLAS | Fedora 44, OpenBLAS 0.3.29 en série à la place de `libRblas.so`, LAPACK du système (3.12.0) gardée | écart numérique au BLAS de référence |
| noSuggests | `_R_CHECK_DEPENDS_ONLY_=true` : les paquets de `Suggests` (testthat, knitr…) absents | test ou vignette qui suppose un paquet suggéré |
| gcc-ASAN, gcc-UBSAN | Fedora 44, gcc 16.1, `-fsanitize=address,undefined,bounds-strict`, R construit avec | lecture ou écriture hors bornes, comportement indéfini dans le C++ |
| clang-ASAN, clang-UBSAN | Fedora 44, clang 23 et libc++, `-fsanitize=address,undefined` | idem, avec la bibliothèque C++ de LLVM |
| valgrind | Fedora 44, valgrind 3.27.1, R instrumenté (niveau 2), `R CMD check --use-valgrind` | mémoire non initialisée, accès invalides |
| noLD | R-devel configuré `--disable-long-double` | sommes et moyennes de R calculées sans précision étendue, comme sur le Mac ARM |
| M1mac | R-devel sur Mac arm64, macOS 27, LAPACK interne | arithmétique et compilateur du Mac ARM |
| LTO | Fedora, `--enable-lto=R`, `INSTALL --use-LTO` | prototypes incohérents entre unités de compilation (avertissement `-Wlto-type-mismatch`) |
| clang23, gcc | compilateurs les plus récents (en-têtes de libc++ allégés) | en-tête C++ manquant |
| rchk | analyse statique du code natif (PROTECT) | protection oubliée ; peu pertinent avec Rcpp, que l'outil analyse mal |

Écartés, avec leur raison :

| Type | Raison |
|---|---|
| ATLAS | « no longer run regularly » depuis le 2026-09-01 (README des BLAS) |
| BLIS | serial BLIS 2.0 de Fedora ; aucune copie ici ni conteneur R-hub : [non reproduit]{.pastille .ko}, risque couvert par les trois autres BLAS |
| Intel | page « no longer updated » (compilateurs oneAPI 2024.2.1) |
| Strict, noRemap, C23, 0len | devenus le comportement par défaut de R-devel, ou archivés : tout contrôle R-devel les couvre |
| BLAS | appels Fortran sans `FCONE` : WH les a, et R ≥ 4.2 l'impose à la compilation |
| noOMP, donttest | WH n'a ni OpenMP ni `\donttest` |
| musl, linux-arm64, rcnst, rlibro, vnu | sans lien avec le code de WH ; vnu demande Java, absent ici |

## Où rejouer chaque contrôle : Windows, WSL ou R-hub

Trois lieux sont possibles. **Windows** : les copies privées de R 4.5.3 sous `C:\Users\Guillaume Biessy\R-prive\`,
posées par le chantier MKL d'ebmcd, et Rtools45. **WSL** : la distribution Ubuntu 20.04 de la machine. **R-hub** :
le service de contrôle de R-hub, qui tourne sur GitHub Actions dans le dépôt du paquet ; son flux
`.github/workflows/rhub.yaml` est déjà dans WH, le dépôt est public, l'exécution est gratuite.

Trois constats, mesurés le 25/10, tranchent la plupart des cas :

- **Rtools45 ne sait pas compiler avec les sanitizers.** `g++ -fsanitize=address` et `-fsanitize=undefined` échouent à
  l'édition de liens (`cannot find -lasan`, `-lubsan`) ; seul `-fsanitize-undefined-trap-on-error` compile, qui arrête
  le programme sans dire où. Valgrind n'existe pas sous Windows.
- **Le WSL ne compile rien en l'état** : ni `gcc`, ni `g++`, ni `make` ; R 4.3.1 sans testthat ni Rcpp. Ubuntu 20.04
  n'a que gcc 9 et valgrind 3.15, loin des gcc 16, clang 23 et valgrind 3.27 du CRAN.
- **Les noyaux de MKL en série se comportent de même sous Windows et sous Linux** (chantier MKL d'ebmcd,
  `dev/mkl/note-mkl.md`) : le défaut propre à Linux, la couche de fils Intel, ne touche pas la configuration en série
  du CRAN.

::: definition
**Substitut.** La configuration de cette machine ou de R-hub qui rejoue le mieux un contrôle du CRAN : même
bibliothèque à la même version, mêmes options, et le chemin d'instructions du processeur du CRAN quand la bibliothèque
en dépend. Sa **fidélité** dit ce qu'il ne reproduit pas.
:::

| Contrôle du CRAN | Substitut retenu | Lieu | Ce qu'il ne reproduit pas |
|---|---|---|---|
| r-release-windows | R 4.5.3, BLAS et LAPACK de référence, Rtools45 (GCC 14.2) | Windows | GCC 14.3 du CRAN ; processeur AMD |
| r-devel-windows | R-devel Windows | Windows | [à télécharger]{.pastille .hyp} : environ 90 Mo |
| MKL | copie `R-4.5.3-mkl` (oneMKL 2026.1, la version du CRAN), `MKL_THREADING_LAYER=SEQUENTIAL`, `MKL_ENABLE_INSTRUCTIONS=AVX2` | Windows | R 4.5.3 au lieu de R-devel ; Fedora et gcc |
| MKL avant 2026-07 | copie `R-4.5.3-mkl2023` (2023.2), série, AVX2 ; AVX pour le processeur de 2023-2024 | Windows | idem |
| OpenBLAS | copie `R-4.5.3-openblas-cran` : OpenBLAS 0.3.29 à la place de `Rblas.dll` seule, LAPACK de R gardée, `OPENBLAS_NUM_THREADS=1`, `OPENBLAS_CORETYPE=Haswell` | Windows | LAPACK 3.12.1 de R au lieu de la 3.12.0 de Fedora ; OpenBLAS compilé par mingw |
| noSuggests | `_R_CHECK_DEPENDS_ONLY_=true` sur la copie de référence | Windows | rien de propre à Linux |
| gcc-ASAN, gcc-UBSAN | conteneur `gcc-asan` : Fedora 44, R-devel, drapeaux du CRAN, `libasan` préchargée | R-hub | R lui-même non instrumenté |
| clang-ASAN, clang-UBSAN | conteneurs `clang-asan`, `clang-ubsan` : Ubuntu 26.04, clang et libc++ | R-hub | Ubuntu au lieu de Fedora ; `alignment` et `float-divide-by-zero` exclus |
| valgrind | conteneur `valgrind` : Fedora 44, R-devel construit pour valgrind | R-hub | niveau d'instrumentation de R [non vérifié]{.pastille .hyp} |
| noLD | conteneur `nold` : Ubuntu 26.04, R-devel `--disable-long-double` | R-hub | Ubuntu au lieu de Fedora |
| M1mac, macOS arm64 | plate-forme `macos-arm64` : Mac M1 de GitHub, R-devel | R-hub | macOS et réglages de GitHub, pas ceux de B. Ripley |
| LTO | conteneur `lto` | R-hub | — |
| clang23, gcc 16, Linux | conteneurs `clang23`, `gcc16`, `ubuntu-clang`, `ubuntu-gcc12` | R-hub | Ubuntu au lieu de Debian testing |
| rchk | conteneur `rchk` | R-hub | — |

Pour MKL, R-hub est **moins** fidèle que Windows : son conteneur `mkl` pose oneMKL 2025.3 (pas 2026.1) par `mkl_rt`
(`containers/mkl/Dockerfile` de R-hub, lu le 25/10), et les machines de GitHub peuvent être des AMD EPYC, sur lesquels MKL
ne prend pas les noyaux d'un Xeon [hypothèse]{.pastille .hyp}. Windows permet en outre d'imposer le chemin AVX2 du
processeur du CRAN.

### Le WSL, chiffré

Le WSL ne serait le bon lieu que pour rejouer **hors ligne** ce que R-hub fournit. Il demanderait, avant tout contrôle :
les outils de compilation et valgrind d'Ubuntu (`build-essential`, `gfortran`, `r-base-dev`, `valgrind`, de l'ordre de
300 Mo [estimation]{.pastille .hyp}), les paquets R Rcpp et testthat compilés depuis leurs sources, puis, pour noLD et
un valgrind instrumenté, deux constructions de R-devel depuis ses sources (40 Mo de sources, une vingtaine de minutes
chacune [estimation]{.pastille .hyp}). Le résultat resterait sur gcc 9 et valgrind 3.15, Ubuntu 20.04 n'étant plus
maintenue. R-hub donne les mêmes contrôles, à jour, sans rien poser ici.

## Ce que donne WH 2.0.0

Archive construite depuis `origin/master`, dont les fichiers du paquet sont ceux de la version 2.0.0 sur le CRAN.
Chaque ligne est un `R CMD check --as-cran` complet (installation, exemples, 76 tests, vignette, manuel), le
25/10/2026. La bannière est celle que la bibliothèque imprime elle-même au chargement, dans un processus lancé comme
ceux du contrôle.

| Substitut | Bibliothèque chargée | Statut | Tests |
|---|---|---|---|
| référence (r-release-windows) | BLAS et LAPACK de R | 1 NOTE | 76 / 76 |
| noSuggests | BLAS et LAPACK de R | 1 NOTE | 76 / 76 |
| OpenBLAS | OpenBLAS 0.3.29, noyaux Haswell, un fil | 1 NOTE | 76 / 76 |
| MKL | oneMKL 2026 Update 1, AVX-512, série | 1 NOTE | 76 / 76 |
| MKL, processeur du CRAN | oneMKL 2026 Update 1, AVX2, série | 1 NOTE | 76 / 76 |
| MKL avant 2026-07 | oneMKL 2023.0 Update 2, AVX2, série | 1 NOTE | 76 / 76 |
| MKL de 2023-2024 | oneMKL 2023.0 Update 2, AVX, série | 1 NOTE | 76 / 76 |

La NOTE est partout la même : `CLAUDE.md`, fichier non standard à la racine. Il est apparu dans le dépôt après la
soumission et n'est pas dans l'archive du CRAN ; son inscription au `.Rbuildignore` est un point déjà ouvert de la
doctrine du dépôt. Dans chaque contrôle, WH importe `Rlapack.dll` et rien d'autre.

**WH 2.0.0 passe toutes les configurations de BLAS du CRAN, y compris celle de 2023-2024 qui faisait échouer
WH 1.0.6.**

### Sur R-hub

Lancé le 25/10/2026 avec ton accord, sur `master` (exécution 40191631313 du flux `rhub.yaml`, onglet Actions du
dépôt WH). R-hub fait échouer une tâche si un fichier `.Rout` du contrôle (exemples, tests) contient
« runtime error: » ou « SUMMARY: UndefinedBehaviorSanitizer » ; la ligne de compilation de `compact.cpp` porte bien
les options du sanitizer ou de LTO.

| Plate-forme R-hub | Compilation de WH | Statut | Tests | Sanitizer |
|---|---|---|---|---|
| `gcc-asan` | `-fsanitize=address,undefined,bounds-strict`, `-flto=10` | 1 NOTE | 76 / 76 | aucun signalement |
| `clang-asan` | `-fsanitize=address,undefined`, libc++ | OK | 76 / 76 | aucun signalement |
| `clang-ubsan` | `-fsanitize=undefined`, libc++ | OK | 76 / 76 | aucun signalement |
| `nold` | R-devel sans `long double` | 1 NOTE | 76 / 76 | — |
| `macos-arm64` | Mac M1 de GitHub, R-devel | 1 NOTE | 76 / 76 | — |
| `lto` | `-flto`, aucun avertissement `-Wlto-type-mismatch` | 2 NOTE | 76 / 76 | — |
| `clang23` | clang 23 | 1 NOTE | 76 / 76 | — |
| `gcc16` | gcc 16 | 1 NOTE | 76 / 76 | — |
| `valgrind` | R-devel sous valgrind | sans verdict : annulé après 6 h 01 (25/10, 12 h 06 → 18 h 07 UTC), la durée maximale d'une tâche GitHub | — | — |

La NOTE commune est celle du `CLAUDE.md`. La seconde NOTE de `lto` dit que les tests ont pris 3,3 fois plus de temps
de calcul que de temps écoulé : du calcul parallèle, que WH ne fait pas. Le conteneur charge un OpenBLAS à plusieurs
fils (`openblas-pthread` 0.3.32 ; FlexiBLAS OpenMP pour `gcc-asan`, d'après le `sessionInfo()` des journaux), là où les
environnements réguliers du CRAN utilisent le BLAS de référence. C'est un effet du conteneur, et une raison de plus
pour rejouer les contrôles de BLAS sous Windows plutôt que sur R-hub.

**WH 2.0.0 passe aussi les sanitizers, noLD, le Mac ARM, LTO et les compilateurs les plus récents.**

## La recette

Les scripts sont dans `documentation/cran-tests/` (hors de l'archive du paquet, `^documentation$` étant au
`.Rbuildignore`). Une session les lance sous le verrou machine (skill `gate-lock-r`), depuis Git Bash :

```
bash documentation/cran-tests/run_windows.sh <dépôt> <dossier de sortie> [<référence git>] [<calibration | none>]
```

- `configs.sh` : un substitut par configuration (copie privée de R et variables d'environnement), sources en tête.
- `check_config.sh` : un `R CMD check --as-cran` sous un substitut ; écrit une ligne de bilan (statut, ligne de
  testthat, bannière de la bibliothèque chargée, DLL importées par le paquet compilé).
- `run_windows.sh` : construit l'archive depuis la référence git (par défaut `origin/master`, jamais l'arbre de
  travail), rejoue la calibration, puis toutes les configurations.
- `blas_probe.R` : fait annoncer à la bibliothèque sa version, son chemin d'instructions et sa couche de fils
  (`MKL_VERBOSE=1`, `OPENBLAS_VERBOSE=2`).
- `calibrer_nosuggests.sh` : prouve que `_R_CHECK_DEPENDS_ONLY_` prend effet, sur un paquet témoin.

Les calibrations précèdent tout verdict, parce qu'un contrôle qui passe ne prouve rien tant que l'instrument n'a pas
montré qu'il voit un défaut connu :

| Calibration | Attendu | Mesuré le 25/10 |
|---|---|---|
| WH 1.0.6 (R pur, échec MKL de 2023) sous MKL 2023.2, série, AVX | des tests en échec | FAIL 5, PASS 179 |
| WH 1.0.6 sous MKL 2026.1, série, AVX-512 | aucun échec | FAIL 0, PASS 184 |
| témoin dont l'exemple appelle un paquet de `Suggests`, sans puis avec `_R_CHECK_DEPENDS_ONLY_` | OK, puis erreur | `Status: OK`, puis « there is no package called 'spelling' » |

Le CRAN avait vu 7 échecs en 2023 ; le chantier MKL d'ebmcd en mesurait 5 à 9 selon les tirages d'alignement.

::: ecart titre="Une copie de R à BLAS remplacé peut faire charger une autre bibliothèque que celle qu'on teste"
Le premier passage du banc a rendu WH 40 fois plus lent sous OpenBLAS (exemples en 234 s, contre 13 s ailleurs).
Cause mesurée : lié contre la copie OpenBLAS, l'éditeur de liens de Rtools prend `dpbtrf`, `dtbtrs` et `dsbev` dans
la DLL OpenBLAS posée en `Rblas.dll` (elle exporte aussi LAPACK), et les importe sous son nom interne,
`libopenblas.dll`. Windows cherche alors ce nom sur le PATH, et y trouve celle d'Inkscape : le paquet tournait sur un
autre OpenBLAS que celui du contrôle. Correctif, dans `check_config.sh` : lier le paquet aux noms `Rblas.dll` et
`Rlapack.dll` de la copie de référence (`R_MAKEVARS_USER`) ; au chargement, ces noms désignent les DLL déjà chargées
par la copie testée, comme `-llapack -lRblas` au CRAN. Le banc vérifie désormais les DLL importées par le paquet.

Le R de production (`C:\Program Files\R\R-4.5.3`) a OpenBLAS à la fois en `Rblas.dll` et en `Rlapack.dll` : tout
paquet qui appelle LAPACK, compilé par lui, importera `libopenblas.dll` et chargera celui d'Inkscape. Aucun paquet
installé n'est touché aujourd'hui (aucun n'importe ce nom, contrôle du 25/10), mais le prochain compilé le sera.
:::

## Le skill « CRAN-tests » proposé

Une proposition, à valider : la création passera ensuite par le skill `helper-skill-creator`.

**Quand il s'arme** : avant une soumission au CRAN d'un paquet maison (WH, WHfull, ebmcd, Linkplots…), après un
signalement du CRAN dans un contrôle additionnel, et après une montée de version de R (les copies privées sont à
refaire).

**Deux étages.**

| Étage | Ce qu'il rejoue | Ce qu'il demande |
|---|---|---|
| Windows, sur cette machine | référence, noSuggests, OpenBLAS, MKL 2026.1 (AVX-512 et AVX2), MKL 2023.2 (AVX2 et AVX) ; R-devel Windows si posé | rien : lancé sous le verrou, environ 25 min, verdict en tableau |
| R-hub, sur GitHub | `gcc-asan`, `clang-asan`, `clang-ubsan`, `valgrind`, `nold`, `macos-arm64`, `lto`, `clang23`, `gcc16` ; `rchk` pour un paquet en C | ton accord à chaque soumission : les calculs tournent sur GitHub et leurs journaux sont visibles dans l'onglet Actions d'un dépôt public |

**Scripts**, repris de `documentation/cran-tests/` et généralisés à tout paquet :

- `configs.sh`, `check_config.sh`, `run_local.sh` (l'actuel `run_windows.sh`, avec la liste des configurations en
  option) et `blas_probe.R` ;
- les calibrations : `calibrer_mkl.sh` (WH 1.0.6 doit échouer sous MKL 2023.2 AVX), `calibrer_nosuggests.sh` ;
- `rhub_lancer.sh` : `gh workflow run rhub.yaml -f config=<liste>` sur la branche visée, puis une seule surveillance
  jusqu'à la fin des tâches ; `rhub_bilan.sh` : le statut de chaque tâche et l'extrait d'erreur ;
- `poser_copies.sh` : refait les copies privées après une montée de R (repris de `mkl_prive.sh` et
  `install_prive.sh` d'ebmcd), copie OpenBLAS « BLAS seul » comprise.

**Règles qu'il porte** :

- l'archive se construit depuis la référence git soumise, jamais depuis l'arbre de travail ;
- aucun verdict sans les deux calibrations vertes ; un paquet compilé qui importe autre chose que `Rblas.dll` et
  `Rlapack.dll` invalide son contrôle ;
- les durées ne se jugent pas (la machine est partagée), ni sur R-hub le rapport temps de calcul / temps écoulé, que
  gonfle l'OpenBLAS à plusieurs fils de ses conteneurs ;
- un succès R-hub se lit dans le journal de chaque tâche (API `actions/jobs/<id>/logs`, servi dès la fin de la tâche,
  sans attendre toute l'exécution) : statut du contrôle, ligne de testthat, options de compilation ;
- les sources du CRAN se relisent à chaque usage (README des BLAS : versions de MKL et d'OpenBLAS) : un changement de
  version demande un téléchargement, donc ton accord ;
- un paquet sans flux R-hub reçoit `.github/workflows/rhub.yaml` (hors archive, `^\.github$` au `.Rbuildignore`) ;
  pour un dépôt privé, R-hub consomme le quota gratuit de GitHub Actions du compte.

## Ce qui reste à décider

::: question titre="Valider la conception du skill « CRAN-tests »"
Périmètre, deux étages, scripts et règles de la section précédente. Une fois validée, le skill se crée par
`helper-skill-creator`, et reprend les scripts de `documentation/cran-tests/`.
:::

Deux points ne demandent pas de décision maintenant :

- **R-devel sous Windows** (environ 90 Mo [estimation]{.pastille .hyp}) n'apporte qu'une version de R plus récente
  sur des contrôles déjà verts ; le skill le proposera à la première soumission qui en aura besoin.
- **Le R de production sous OpenBLAS** (encadré de la recette) relève du dépôt ebmcd, où vit `poser_openblas.ps1` :
  une session à lancer là-bas.

Les scripts restent sur la branche `claude/cran-configs-ccce4a` : ce ne sont pas des fichiers documentaires, et
la doctrine de WH interdit de les pousser sur `master` sans suite de tests. Le skill les reprendra ; seule cette note
va sur `master`.
