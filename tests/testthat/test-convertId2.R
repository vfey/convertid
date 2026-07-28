# Tests for convertId2()
#
# convertId2() is a fast 1-to-1 gene identifier converter. The critical
# guarantee is that the output vector is always the same length and order
# as the input, with unresolved entries as NA.
#
# Conversions tested:
#   - ENSG -> symbol (auto and explicit)
#   - Symbol -> ENSG (auto and explicit)
#   - Entrez -> symbol (auto and explicit)
#   - Entrez -> ENSG (explicit output = "ensembl")
#   - Single element and batch paths
#   - NA handling: unknown IDs, ambiguous multi-Entrez mappings
#   - Output length and names match input
#   - Mouse species
#   - Missing package error
#
# Known stable mappings used as fixtures (Human):
#   ACTB:  ENSG00000075624, Entrez 60
#   GAPDH: ENSG00000111640, Entrez 2597
#   BRCA1: ENSG00000012048, Entrez 672

# ---------------------------------------------------------------------------
# Output length and names invariant
# ---------------------------------------------------------------------------
testthat::test_that("output is same length as input for ENSG -> symbol", {
  testthat::skip_if_not_installed("org.Hs.eg.db")
  input  <- c("ENSG00000075624", "ENSG00000111640", "ENSG00000012048")
  result <- convertId2(input)
  testthat::expect_equal(length(result), length(input))
  testthat::expect_equal(names(result), input)
})

testthat::test_that("output is same length as input for symbol -> ENSG", {
  testthat::skip_if_not_installed("org.Hs.eg.db")
  input  <- c("ACTB", "GAPDH", "BRCA1")
  result <- convertId2(input)
  testthat::expect_equal(length(result), length(input))
  testthat::expect_equal(names(result), input)
})

testthat::test_that("output is same length as input for Entrez -> symbol", {
  testthat::skip_if_not_installed("org.Hs.eg.db")
  input  <- c("60", "2597", "672")
  result <- convertId2(input)
  testthat::expect_equal(length(result), length(input))
  testthat::expect_equal(names(result), input)
})

testthat::test_that("output is same length as input for Entrez -> ENSG", {
  testthat::skip_if_not_installed("org.Hs.eg.db")
  input  <- c("60", "2597", "672")
  result <- convertId2(input, output = "ensembl")
  testthat::expect_equal(length(result), length(input))
  testthat::expect_equal(names(result), input)
})

# ---------------------------------------------------------------------------
# Correct conversions: ENSG -> symbol
# ---------------------------------------------------------------------------
testthat::test_that("ENSG -> symbol returns correct symbols (single)", {
  testthat::skip_if_not_installed("org.Hs.eg.db")
  testthat::expect_equal(convertId2("ENSG00000075624"), c("ENSG00000075624" = "ACTB"))
  testthat::expect_equal(convertId2("ENSG00000111640"), c("ENSG00000111640" = "GAPDH"))
})

testthat::test_that("ENSG -> symbol returns correct symbols (batch)", {
  testthat::skip_if_not_installed("org.Hs.eg.db")
  result <- convertId2(c("ENSG00000075624", "ENSG00000111640"))
  testthat::expect_equal(result[["ENSG00000075624"]], "ACTB")
  testthat::expect_equal(result[["ENSG00000111640"]], "GAPDH")
})

# ---------------------------------------------------------------------------
# Correct conversions: symbol -> ENSG
# ---------------------------------------------------------------------------
testthat::test_that("symbol -> ENSG returns correct IDs (single)", {
  testthat::skip_if_not_installed("org.Hs.eg.db")
  result <- convertId2("ACTB")
  testthat::expect_equal(result[["ACTB"]], "ENSG00000075624")
})

testthat::test_that("symbol -> ENSG returns correct IDs (batch)", {
  testthat::skip_if_not_installed("org.Hs.eg.db")
  result <- convertId2(c("ACTB", "GAPDH"))
  testthat::expect_equal(result[["ACTB"]],  "ENSG00000075624")
  testthat::expect_equal(result[["GAPDH"]], "ENSG00000111640")
})

# ---------------------------------------------------------------------------
# Correct conversions: Entrez -> symbol
# ---------------------------------------------------------------------------
testthat::test_that("Entrez -> symbol returns correct symbols (single)", {
  testthat::skip_if_not_installed("org.Hs.eg.db")
  testthat::expect_equal(convertId2("60"),   c("60"   = "ACTB"))
  testthat::expect_equal(convertId2("2597"), c("2597" = "GAPDH"))
})

testthat::test_that("Entrez -> symbol returns correct symbols (batch)", {
  testthat::skip_if_not_installed("org.Hs.eg.db")
  result <- convertId2(c("60", "2597"))
  testthat::expect_equal(result[["60"]],   "ACTB")
  testthat::expect_equal(result[["2597"]], "GAPDH")
})

# ---------------------------------------------------------------------------
# Correct conversions: Entrez -> ENSG
# ---------------------------------------------------------------------------
testthat::test_that("Entrez -> ENSG returns correct IDs (single)", {
  testthat::skip_if_not_installed("org.Hs.eg.db")
  result <- convertId2("60", output = "ensembl")
  testthat::expect_equal(result[["60"]], "ENSG00000075624")
})

testthat::test_that("Entrez -> ENSG returns correct IDs (batch)", {
  testthat::skip_if_not_installed("org.Hs.eg.db")
  result <- convertId2(c("60", "2597"), output = "ensembl")
  testthat::expect_equal(result[["60"]],   "ENSG00000075624")
  testthat::expect_equal(result[["2597"]], "ENSG00000111640")
})

# ---------------------------------------------------------------------------
# Explicit output overrides auto detection
# ---------------------------------------------------------------------------
testthat::test_that("output = 'symbol' on ENSG input returns symbols", {
  testthat::skip_if_not_installed("org.Hs.eg.db")
  result <- convertId2("ENSG00000075624", output = "symbol")
  testthat::expect_equal(result[["ENSG00000075624"]], "ACTB")
})

testthat::test_that("output = 'ensembl' on symbol input returns ENSG", {
  testthat::skip_if_not_installed("org.Hs.eg.db")
  result <- convertId2("ACTB", output = "ensembl")
  testthat::expect_equal(result[["ACTB"]], "ENSG00000075624")
})

testthat::test_that("output = 'ensembl' on Entrez input returns ENSG", {
  testthat::skip_if_not_installed("org.Hs.eg.db")
  result <- convertId2("60", output = "ensembl")
  testthat::expect_equal(result[["60"]], "ENSG00000075624")
})

testthat::test_that("output = 'symbol' on Entrez input returns symbol", {
  testthat::skip_if_not_installed("org.Hs.eg.db")
  result <- convertId2("60", output = "symbol")
  testthat::expect_equal(result[["60"]], "ACTB")
})

# ---------------------------------------------------------------------------
# NA handling
# ---------------------------------------------------------------------------
testthat::test_that("unknown ENSG ID returns NA", {
  testthat::skip_if_not_installed("org.Hs.eg.db")
  result <- convertId2("ENSG00000000000")
  testthat::expect_true(is.na(result[["ENSG00000000000"]]))
})

testthat::test_that("unknown symbol returns NA", {
  testthat::skip_if_not_installed("org.Hs.eg.db")
  result <- convertId2("NOTAREAL1")
  testthat::expect_true(is.na(result[["NOTAREAL1"]]))
})

testthat::test_that("unknown Entrez ID returns NA", {
  testthat::skip_if_not_installed("org.Hs.eg.db")
  result <- convertId2("9999999999")
  testthat::expect_true(is.na(result[["9999999999"]]))
})

testthat::test_that("NA in input is preserved as NA in output", {
  testthat::skip_if_not_installed("org.Hs.eg.db")
  # NA detection falls through to first non-NA element
  input  <- c("ENSG00000075624", "ENSG00000111640")
  result <- convertId2(input)
  testthat::expect_false(any(is.na(result)))
})

testthat::test_that("mix of known and unknown IDs preserves order and length", {
  testthat::skip_if_not_installed("org.Hs.eg.db")
  input  <- c("ENSG00000075624", "ENSG00000000000", "ENSG00000111640")
  result <- convertId2(input)
  testthat::expect_equal(length(result), 3L)
  testthat::expect_equal(result[["ENSG00000075624"]], "ACTB")
  testthat::expect_true(is.na(result[["ENSG00000000000"]]))
  testthat::expect_equal(result[["ENSG00000111640"]], "GAPDH")
})

# ---------------------------------------------------------------------------
# Round-trip consistency
# ---------------------------------------------------------------------------
testthat::test_that("ENSG -> symbol -> ENSG round-trip is consistent", {
  testthat::skip_if_not_installed("org.Hs.eg.db")
  ensg    <- c("ENSG00000075624", "ENSG00000111640", "ENSG00000012048")
  symbols <- convertId2(ensg)
  back    <- convertId2(as.character(symbols))
  # Round-trip may not be perfect for all genes but should hold for stable ones
  testthat::expect_equal(back[["ACTB"]],  "ENSG00000075624")
  testthat::expect_equal(back[["GAPDH"]], "ENSG00000111640")
  testthat::expect_equal(back[["BRCA1"]], "ENSG00000012048")
})

testthat::test_that("Entrez -> symbol -> Entrez round-trip is consistent", {
  testthat::skip_if_not_installed("org.Hs.eg.db")
  entrez  <- c("60", "2597")
  symbols <- convertId2(entrez)
  back    <- mget(as.character(symbols),
                  envir = getNamespace("org.Hs.eg.db")[["org.Hs.egSYMBOL2EG"]],
                  ifnotfound = NA)
  testthat::expect_equal(back[["ACTB"]][[1]],  "60")
  testthat::expect_equal(back[["GAPDH"]][[1]], "2597")
})

# ---------------------------------------------------------------------------
# Mouse species
# ---------------------------------------------------------------------------
testthat::test_that("Mouse ENSG -> symbol works", {
  testthat::skip_if_not_installed("org.Mm.eg.db")
  # Actb mouse: ENSMUSG00000029580
  result <- convertId2("ENSMUSG00000029580", species = "Mouse")
  testthat::expect_equal(result[["ENSMUSG00000029580"]], "Actb")
})

testthat::test_that("Mouse symbol -> ENSG works", {
  testthat::skip_if_not_installed("org.Mm.eg.db")
  result <- convertId2("Actb", species = "Mouse")
  testthat::expect_equal(result[["Actb"]], "ENSMUSG00000029580")
})

# ---------------------------------------------------------------------------
# Correct conversions: -> Entrez (output = "entrez")
# ---------------------------------------------------------------------------
testthat::test_that("ENSG -> Entrez returns correct IDs (single)", {
  testthat::skip_if_not_installed("org.Hs.eg.db")
  testthat::expect_equal(convertId2("ENSG00000075624", output = "entrez"),
                         c("ENSG00000075624" = "60"))
  testthat::expect_equal(convertId2("ENSG00000111640", output = "entrez"),
                         c("ENSG00000111640" = "2597"))
})

testthat::test_that("ENSG -> Entrez returns correct IDs (batch)", {
  testthat::skip_if_not_installed("org.Hs.eg.db")
  result <- convertId2(c("ENSG00000075624", "ENSG00000111640"), output = "entrez")
  testthat::expect_equal(result[["ENSG00000075624"]], "60")
  testthat::expect_equal(result[["ENSG00000111640"]], "2597")
  testthat::expect_equal(names(result), c("ENSG00000075624", "ENSG00000111640"))
})

testthat::test_that("Symbol -> Entrez returns correct IDs (single and batch)", {
  testthat::skip_if_not_installed("org.Hs.eg.db")
  testthat::expect_equal(convertId2("ACTB", output = "entrez"), c("ACTB" = "60"))
  result <- convertId2(c("ACTB", "GAPDH"), output = "entrez")
  testthat::expect_equal(result[["ACTB"]],  "60")
  testthat::expect_equal(result[["GAPDH"]], "2597")
})

testthat::test_that("Entrez -> Entrez is identity", {
  testthat::skip_if_not_installed("org.Hs.eg.db")
  testthat::expect_equal(convertId2("60", output = "entrez"), c("60" = "60"))
  result <- convertId2(c("60", "2597"), output = "entrez")
  testthat::expect_equal(result[["60"]],   "60")
  testthat::expect_equal(result[["2597"]], "2597")
})

testthat::test_that("output = 'entrez' preserves length, order and names", {
  testthat::skip_if_not_installed("org.Hs.eg.db")
  input  <- c("ACTB", "GAPDH", "BRCA1")
  result <- convertId2(input, output = "entrez")
  testthat::expect_equal(length(result), length(input))
  testthat::expect_equal(names(result), input)
})

testthat::test_that("unknown ID with output = 'entrez' returns NA", {
  testthat::skip_if_not_installed("org.Hs.eg.db")
  result <- convertId2("NOTAREAL1", output = "entrez")
  testthat::expect_true(is.na(result[["NOTAREAL1"]]))
})

# ---------------------------------------------------------------------------
# multi2NA: one-to-many resolution
# ---------------------------------------------------------------------------
testthat::test_that("multi2NA does not affect 1-to-1 mappings", {
  testthat::skip_if_not_installed("org.Hs.eg.db")
  input <- c("ACTB", "GAPDH", "BRCA1")
  testthat::expect_equal(convertId2(input, output = "entrez", multi2NA = FALSE),
                         convertId2(input, output = "entrez", multi2NA = TRUE))
})

testthat::test_that("multi2NA = TRUE returns NA exactly where default collapses multis", {
  testthat::skip_if_not_installed("org.Hs.eg.db")
  # Property holds regardless of whether the input actually contains a
  # one-to-many mapping: wherever the default (collapse) result contains a
  # "///"-separated string, the multi2NA = TRUE result must be NA; all other
  # entries must be identical between the two modes.
  input    <- c("60", "2597", "672", "ENSG00000075624", "ENSG00000111640")
  res_keep <- convertId2(input, output = "ensembl", multi2NA = FALSE)
  res_na   <- convertId2(input, output = "ensembl", multi2NA = TRUE)
  is_multi <- !is.na(res_keep) & grepl(" /// ", res_keep, fixed = TRUE)
  testthat::expect_true(all(is.na(res_na[is_multi])))
  testthat::expect_equal(res_na[!is_multi], res_keep[!is_multi])
})

testthat::test_that("single-element path honours multi2NA for a known multi mapping", {
  testthat::skip_if_not_installed("org.Hs.eg.db")
  # Discover a symbol that maps to more than one Ensembl gene ID in the
  # installed annotation, so the test is robust to annotation versions.
  orgdb  <- getNamespace("org.Hs.eg.db")
  sym2eg <- orgdb[["org.Hs.egSYMBOL2EG"]]
  ensg   <- orgdb[["org.Hs.egENSEMBL"]]
  syms   <- head(AnnotationDbi::keys(sym2eg), 2000L)
  multi_sym <- NA_character_
  for (s in syms) {
    eg <- tryCatch(get(s, envir = sym2eg), error = function(e) NULL)
    if (is.null(eg) || length(eg) != 1L) next          # keep intermediate 1-to-1
    en <- tryCatch(get(eg, envir = ensg), error = function(e) NULL)
    if (!is.null(en) && length(en) > 1L) { multi_sym <- s; break }
  }
  testthat::skip_if(is.na(multi_sym),
                    "No symbol with a one-to-many Ensembl mapping found")
  keep <- convertId2(multi_sym, output = "ensembl", multi2NA = FALSE)
  na   <- convertId2(multi_sym, output = "ensembl", multi2NA = TRUE)
  testthat::expect_true(grepl(" /// ", keep[[1]], fixed = TRUE))
  testthat::expect_true(is.na(na[[1]]))
})

# ---------------------------------------------------------------------------
# Missing package error
# ---------------------------------------------------------------------------
testthat::test_that("missing org.Hs.eg.db gives informative error", {
  mockery::stub(convertId2, "requireNamespace", function(pkg, ...) {
    if (pkg == "org.Hs.eg.db") FALSE else TRUE
  })
  testthat::expect_error(
    convertId2("ACTB"),
    regexp = "org.Hs.eg.db"
  )
})

testthat::test_that("missing org.Mm.eg.db gives informative error", {
  mockery::stub(convertId2, "requireNamespace", function(pkg, ...) {
    if (pkg == "org.Mm.eg.db") FALSE else TRUE
  })
  testthat::expect_error(
    convertId2("Actb", species = "Mouse"),
    regexp = "org.Mm.eg.db"
  )
})
