# PRE-ETAPE DE PROMOTION : NAMESPACE, man/ et DESCRIPTION suivent-ils le source ?
#
# Decision de l'auteur (2026-10-31) : « on ne lance plus document que lors de la
# promotion ». Entre deux promotions, le controle statique
# tests/testthat/test-roxygen-exports.R garde les exports ; a la promotion, ce
# script regenere la doc roxygen sur une COPIE de l'arbre (l'arbre du lander
# n'est jamais sali) et refuse si elle differe de la doc committee. Remede :
# `devtools::document()`, puis un commit a part `docs(roxygen): ...`.
#
# Sur copie, le code se charge par `roxygen2::load_source` : il source `R/`
# sans compiler `src/`, si bien que le cout reste celui de roxygen seul, meme
# pour un paquet compile. Calibre le 2026-10-31 : 0 ecart contre la doc produite
# par document() sur Linkplots, ebmcd, WHfull et WH.
#
# Appel : `source("dev/verifier_roxygen.R"); verifier_roxygen(".")` -- stop()
# si ecart, avec la liste des fichiers. Kill-probe : `Rscript dev/verifier_roxygen.R
# --killprobe` retire un @export sur une copie et exige le refus.
# COPIE IDENTIQUE dans les depots qui promeuvent (Linkplots, ebmcd, WHfull, WH).
#
# QUAND : `roxygen_a_verifier()` dit si le gate courant doit verifier. Le lander
# (`land-to-master.sh`) exporte a la commande de gate `LAND_GATE_BASE` (base visee)
# et `LAND_DECLARED_BASE` (base de dev/land-base.txt, vide si aucune). On verifie a
# la PROMOTION (les deux posees et differentes), dans un depot sans base declaree,
# et quand la suite est lancee a la main (`LAND_GATE_BASE` absente) ; on saute sur
# un atterrissage vers la base declaree. Pourquoi : WH lance sa suite a chaque
# atterrissage sur `dev`, ou la verification refusait des `dev` (constate le
# 2026-10-31), contre la decision de l'auteur.

roxygen_a_verifier <- function(gate = Sys.getenv("LAND_GATE_BASE", NA),
                               declaree = Sys.getenv("LAND_DECLARED_BASE", "")) {
  if (is.na(gate) || !nzchar(gate) || !nzchar(declaree)) return(TRUE)
  !identical(gate, declaree)
}

verifier_roxygen <- function(racine = ".", copie_hook = NULL) {
  racine <- normalizePath(racine, winslash = "/")
  if (!requireNamespace("roxygen2", quietly = TRUE))
    stop("roxygen2 absent : la doc ne peut pas etre verifiee avant la promotion.")
  tmp <- tempfile("verif-roxygen-")
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  copie <- file.path(tmp, basename(racine))
  dir.create(copie)
  for (x in c("DESCRIPTION", "NAMESPACE", "R", "man", "inst", "data"))
    if (file.exists(file.path(racine, x)))
      file.copy(file.path(racine, x), copie, recursive = TRUE)
  if (is.function(copie_hook)) copie_hook(copie)

  # load_source ne charge pas `data/` : un jeu documente par son nom (`"immo"`)
  # doit exister dans l'environnement que roxygen inspecte.
  charger <- function(chemin) {
    env <- roxygen2::load_source(chemin)
    for (f in list.files(file.path(chemin, "data"), pattern = "\\.(rda|RData)$",
                         full.names = TRUE))
      load(f, envir = env)
    env
  }
  suppressMessages(roxygen2::roxygenise(copie, load_code = charger))

  lire <- function(f) if (file.exists(f)) gsub("\r", "", readLines(f, warn = FALSE)) else NULL
  rd <- function(d) if (dir.exists(d)) list.files(d, pattern = "\\.Rd$") else character(0)
  noms <- c("DESCRIPTION", "NAMESPACE",
            file.path("man", union(rd(file.path(racine, "man")), rd(file.path(copie, "man")))))
  ecarts <- noms[!vapply(noms, function(n)
    identical(lire(file.path(racine, n)), lire(file.path(copie, n))), NA)]
  if (length(ecarts))
    stop(sprintf(paste0("Doc roxygen en retard sur le source (%d fichier(s) : %s). ",
                        "Lancer devtools::document(), committer a part ",
                        "(docs(roxygen): ...), puis relancer la promotion."),
                 length(ecarts), paste(ecarts, collapse = ", ")), call. = FALSE)
  invisible(TRUE)
}

if (sys.nframe() == 0L && identical(commandArgs(TRUE), "--killprobe")) {
  # Arbre sain : vert. Un @export retire sur la copie : refus nommant NAMESPACE.
  t0 <- proc.time()[["elapsed"]]
  verifier_roxygen(".")
  cat(sprintf("temoin : vert (%.1f s)\n", proc.time()[["elapsed"]] - t0))
  retire_export <- function(copie) {
    for (f in list.files(file.path(copie, "R"), full.names = TRUE)) {
      l <- readLines(f, warn = FALSE)
      i <- grep("^#' @export\\s*$", l)
      if (length(i)) { writeLines(l[-i[1]], f); return(invisible()) }
    }
    stop("aucun @export a retirer")
  }
  msg <- tryCatch({ verifier_roxygen(".", retire_export); "VERT" },
                  error = function(e) conditionMessage(e))
  cat("mutant :", msg, "\n")
  if (!grepl("NAMESPACE", msg, fixed = TRUE)) stop("KILL-PROBE MORT : le retrait d'un @export n'a pas ete vu")
  cat("KILL-PROBE OK\n")
}
