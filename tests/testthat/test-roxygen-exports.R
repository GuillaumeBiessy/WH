# CONTROLE STATIQUE DES EXPORTS ROXYGEN -- sans document(), en R de base.
#
# POURQUOI. Une fonction inseree ENTRE un bloc roxygen et la fonction qu'il
# documente capte le bloc ET son `@export` : au document() suivant, l'export
# public disparait du NAMESPACE, en silence -- load_all() ignore le NAMESPACE,
# aucun test ne le voit. Variante soeur : inseree AU MILIEU du bloc, elle en
# capte les `@param` et laisse le NAMESPACE intact. Decision de l'auteur
# (2026-10-31) : document() ne tourne plus qu'a la promotion ; ce controle,
# qui lit les blocs de `R/*.R` contre le NAMESPACE COMMITTE, tourne avec les
# tests, cibles ou complets.
#
# ROUGE SEULEMENT SUR LE PIEGE :
#   (a) un nom du NAMESPACE committe -- `export()`, et `S3method(g, c)` lu
#       `g.c` -- qu'aucun `@export` / `@exportS3Method` ne porte plus dans le
#       source : export perdu (ou retire sans document(), ce qui casse aussi
#       l'installation) ;
#   (b) un `@export` pose sur un objet dont le nom commence par un point ;
#   (c) un `@param` absent des arguments des fonctions de son sujet (union
#       sur un `@rdname` partage) : bloc capte par une fonction inseree.
# Un `@export` NEUF, absent du NAMESPACE, reste vert : le NAMESPACE se
# regenere a la promotion. Un sujet dont une fonction n'est pas un
# `function(...)` litteral (Vectorize, memoise...) n'est pas juge en (c).
#
# COPIE IDENTIQUE dans les sept depots R maison (Linkplots, ebmcd, WHfull, WH,
# arrays, seriatim, seriatimLite) : rien de propre a un paquet ici. Saute sous
# CRAN, et quand les sources `R/` sont absentes (paquet installe, R CMD check).
# Les known-answers du piege sont dans ce fichier : ils prouvent que chaque
# rouge mord encore.

# Les lignes roxygen d'un bloc -> liste nommee de tags (nom, valeur).
.rox_tags <- function(lignes) {
  txt <- sub("^\\s*#' ?", "", lignes)
  debut <- grepl("^@[[:alnum:]_]+", txt)
  if (!any(debut)) return(list())
  grp <- cumsum(debut)
  lapply(split(txt[grp > 0L], grp[grp > 0L]), function(t) {
    nom <- sub("^@([[:alnum:]_]+).*$", "\\1", t[1])
    val <- c(sub("^@[[:alnum:]_]+\\s*", "", t[1]), t[-1])
    list(nom = nom, val = trimws(paste(val, collapse = "\n")))
  })
}

# L'objet qu'un bloc documente : nom, et arguments si c'est un `function`
# litteral (NULL sinon).
.rox_objet <- function(e) {
  if (is.call(e) && is.symbol(e[[1]]) &&
      as.character(e[[1]]) %in% c("<-", "=", "<<-")) {
    lhs <- e[[2]]
    nom <- if (is.symbol(lhs) || is.character(lhs)) as.character(lhs) else NA_character_
    rhs <- e[[3]]
    args <- NULL
    if (is.call(rhs) && identical(rhs[[1]], as.symbol("function")))
      args <- if (is.null(rhs[[2]])) character(0) else names(rhs[[2]])
    return(list(nom = nom, args = args, affectation = TRUE))
  }
  nom <- NA_character_
  if (is.character(e) && length(e) == 1L && !identical(e, "_PACKAGE")) nom <- e
  list(nom = nom, args = NULL, affectation = FALSE)
}

# Tous les blocs roxygen d'un fichier : comme roxygen2, les lignes `#'` entre
# la fin d'une expression de premier niveau et le debut de la suivante forment
# le bloc de cette derniere.
.rox_blocs <- function(fichier) {
  lignes <- readLines(fichier, warn = FALSE, encoding = "UTF-8")
  exprs <- parse(text = lignes, keep.source = TRUE)
  refs <- attr(exprs, "srcref")
  fin_prec <- 0L
  blocs <- list()
  for (i in seq_along(exprs)) {
    deb <- refs[[i]][1]
    zone <- if (deb - 1L > fin_prec) seq.int(fin_prec + 1L, deb - 1L) else integer(0)
    rox <- zone[grepl("^\\s*#'", lignes[zone])]
    fin_prec <- refs[[i]][3]
    if (!length(rox)) next
    obj <- .rox_objet(exprs[[i]])
    blocs[[length(blocs) + 1L]] <- c(obj, list(
      tags = .rox_tags(lignes[rox]),
      lieu = sprintf("%s:%d", basename(fichier), deb)))
  }
  blocs
}

# Noms qu'un bloc exporte, au sens du NAMESPACE (`S3method(g, c)` = `g.c`).
.rox_exports_bloc <- function(b) {
  noms <- vapply(b$tags, `[[`, "", "nom")
  vals <- lapply(b$tags, function(t) strsplit(t$val, "\\s+")[[1]])
  nom_tag <- unlist(vals[noms == "name"])[1]
  defaut <- if (!is.na(b$nom)) b$nom else if (length(nom_tag)) nom_tag else character(0)
  out <- character(0)
  for (k in seq_along(noms)) {
    v <- vals[[k]][nzchar(vals[[k]])]
    if (noms[k] == "export") {
      out <- c(out, if (length(v)) v else defaut)
    } else if (noms[k] == "exportS3Method") {
      if (identical(v, "NULL")) next
      out <- c(out, if (length(v) >= 2L) paste0(sub("^.*::", "", v[1]), ".", v[2]) else defaut)
    }
  }
  if (length(out) && "method" %in% noms) {
    m <- vals[[which(noms == "method")[1]]]
    if (length(m) >= 2L) out <- c(out, paste0(m[1], ".", m[2]))
  }
  unique(out)
}

# Noms exportes par le NAMESPACE committe (`S3method(g, c[, m])` = `g.c`, ou `m`).
.rox_namespace <- function(fichier) {
  nom_arg <- function(a) if (is.character(a)) a else if (is.symbol(a))
    as.character(a) else sub("^.*::", "", paste(deparse(a), collapse = ""))
  out <- character(0)
  for (e in as.list(parse(fichier, keep.source = FALSE))) {
    if (!is.call(e) || !is.symbol(e[[1]])) next
    f <- as.character(e[[1]])
    a <- as.list(e)[-1]
    if (f == "export") out <- c(out, vapply(a, nom_arg, ""))
    if (f == "S3method")
      out <- c(out, if (length(a) >= 3L) nom_arg(a[[3]])
                    else paste0(nom_arg(a[[1]]), ".", nom_arg(a[[2]])))
  }
  unique(out)
}

# Le controle : vecteur de violations (vide = vert).
.rox_controle <- function(racine) {
  fichiers <- list.files(file.path(racine, "R"), pattern = "\\.[Rr]$", full.names = TRUE)
  blocs <- unlist(lapply(fichiers, .rox_blocs), recursive = FALSE)
  viol <- character(0)

  exp_src <- unique(unlist(lapply(blocs, .rox_exports_bloc)))
  perdus <- setdiff(.rox_namespace(file.path(racine, "NAMESPACE")), exp_src)
  if (length(perdus))
    viol <- c(viol, sprintf("(a) export perdu : `%s` est dans le NAMESPACE committe, aucun @export ne le porte dans R/", perdus))

  for (b in blocs) {
    pts <- grep("^\\.", .rox_exports_bloc(b), value = TRUE)
    if (length(pts))
      viol <- c(viol, sprintf("(b) @export sur un objet pointe : `%s` (%s)", pts, b$lieu))
  }

  sujet <- vapply(blocs, function(b) {
    noms <- vapply(b$tags, `[[`, "", "nom")
    for (t in c("rdname", "describeIn", "name")) {
      k <- which(noms == t)
      if (length(k)) return(strsplit(b$tags[[k[1]]]$val, "\\s+")[[1]][1])
    }
    if (is.na(b$nom)) b$lieu else b$nom
  }, "")
  for (s in unique(sujet)) {
    grp <- blocs[sujet == s]
    if (any(vapply(grp, function(b) b$affectation && is.null(b$args), NA))) next
    args <- unique(unlist(lapply(grp, `[[`, "args")))
    if (!any(vapply(grp, function(b) !is.null(b$args), NA))) next
    for (b in grp) {
      p <- unlist(lapply(b$tags, function(t) if (t$nom == "param")
        trimws(strsplit(strsplit(t$val, "\\s+")[[1]][1], ",")[[1]])))
      orph <- setdiff(p[nzchar(p)], c(args, "..."))
      if (length(orph))
        viol <- c(viol, sprintf("(c) @param sans argument : %s, sujet `%s` (%s) -- bloc capte par une fonction inseree, ou @param perime ?",
                                paste0("`", orph, "`", collapse = ", "), s, b$lieu))
    }
  }
  viol
}

test_that("R/ ne perd aucun export du NAMESPACE committe, et aucun bloc n'est capte", {
  skip_on_cran()
  racine <- test_path("..", "..")
  skip_if_not(dir.exists(file.path(racine, "R")) && file.exists(file.path(racine, "NAMESPACE")),
              "sources R/ absentes (paquet installe)")
  v <- .rox_controle(racine)
  if (length(v)) {
    fail(paste(c("Controle roxygen statique :", v,
                 "Remede : inserer une fonction AVANT toute la suite contigue de #', jamais entre un bloc et sa fonction ;",
                 "un export retire expres se constate par document(), a la promotion."), collapse = "\n"))
  } else succeed()
})

# --- Known-answers : le piege reproduit sur un mini-paquet ---------------------
.rox_fixture <- function(a_R) {
  d <- tempfile("rox-exports-")
  dir.create(file.path(d, "R"), recursive = TRUE)
  writeLines(a_R, file.path(d, "R", "a.R"))
  writeLines(c(
    "#' Methode S3",
    "#' @export",
    "print.monobjet <- function(x, ...) invisible(x)",
    "",
    "#' @exportS3Method pkg::autoplot",
    "autoplot.monobjet <- function(object, ...) NULL",
    "",
    "#' Une constante",
    "#' @export",
    "CONST_A <- 3",
    "",
    "#' Un groupe",
    "#' @param u un u",
    "#' @param v un v, argument de g2 seulement",
    "#' @rdname groupe",
    "#' @export",
    "g1 <- function(u) u",
    "",
    "#' @rdname groupe",
    "#' @export",
    "g2 <- function(v, ...) v",
    "",
    "#' Aide interne",
    "#' @param w un w",
    "#' @noRd",
    ".interne <- function(w) w"), file.path(d, "R", "b.R"))
  writeLines(c("# Generated by roxygen2: do not edit by hand", "",
               "S3method(pkg::autoplot,monobjet)", "S3method(print,monobjet)",
               "export(CONST_A)", "export(g1)", "export(g2)", "export(publique)"),
             file.path(d, "NAMESPACE"))
  d
}
.rox_bloc_publique <- c("#' Publique", "#' @param x un x", "#' @param y,z deux autres")
.rox_fn_publique <- "publique <- function(x, y, z) x"

test_that("known-answer : le mini-paquet sain est vert, un @export neuf aussi", {
  skip_on_cran()
  d <- .rox_fixture(c(.rox_bloc_publique, "#' @export", .rox_fn_publique))
  on.exit(unlink(d, recursive = TRUE))
  expect_identical(.rox_controle(d), character(0))
  # Un @export neuf, absent du NAMESPACE : vert (le NAMESPACE suit a la promotion).
  cat("\n#' Neuve\n#' @export\nneuve <- function() 1\n", file = file.path(d, "R", "a.R"), append = TRUE)
  expect_identical(.rox_controle(d), character(0))
})

test_that("known-answer : une fonction inseree entre le bloc et sa fonction allume (a), (b) et (c)", {
  skip_on_cran()
  d <- .rox_fixture(c(.rox_bloc_publique, "#' @export", ".aide <- function(w) w", .rox_fn_publique))
  on.exit(unlink(d, recursive = TRUE))
  v <- .rox_controle(d)
  expect_length(v, 3L)
  expect_true(any(grepl("^\\(a\\) .*`publique`", v)))
  expect_true(any(grepl("^\\(b\\) .*`\\.aide`", v)))
  expect_true(any(grepl("^\\(c\\) .*`x`, `y`, `z`, sujet `\\.aide`", v)))
})

test_that("known-answer : une fonction inseree dans le bloc, avant son @export, allume (c) seul", {
  skip_on_cran()
  d <- .rox_fixture(c(.rox_bloc_publique, ".aide <- function(w) w", "#' @export", .rox_fn_publique))
  on.exit(unlink(d, recursive = TRUE))
  v <- .rox_controle(d)
  expect_length(v, 1L)
  expect_match(v, "^\\(c\\) .*sujet `\\.aide`")
})

test_that("known-answer : une methode S3 qui perd son @export allume (a)", {
  skip_on_cran()
  d <- .rox_fixture(c(.rox_bloc_publique, "#' @export", .rox_fn_publique))
  on.exit(unlink(d, recursive = TRUE))
  b <- readLines(file.path(d, "R", "b.R"))
  writeLines(b[-2L], file.path(d, "R", "b.R"))
  expect_identical(.rox_controle(d),
                   "(a) export perdu : `print.monobjet` est dans le NAMESPACE committe, aucun @export ne le porte dans R/")
})
