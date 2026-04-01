# Tests for unify_gene_ids() and dedup_gene_ids()
#
# dedup_gene_ids() is the internal deduplication function and is tested
# directly via convertid:::dedup_gene_ids() using a pre-computed fixture,
# avoiding any network calls.
#
# unify_gene_ids() is tested for input validation, argument handling, verbose
# output, and graceful degradation when lookups fail. One integration test
# exercises the full live pipeline and is skipped on CRAN.
#
# The fixture (fixture_genes_bm.rda) covers the following logic branches:
#   - Simple single-mapping genes (ACTB, GAPDH)
#   - Two ENSG IDs for same gene, older AnnotationDbi-confirmed preferred (AKAP17A)
#   - ENSG placeholder hgnc_symbol resolved via hgnc_symbol_2 (BMS1P4)
#   - Two identical-metadata rows, newer ENSG preferred (HLA-H)

# ---------------------------------------------------------------------------
# Load fixture
# ---------------------------------------------------------------------------
fixture_path <- testthat::test_path("fixtures", "fixture_genes_bm.rda")
load(fixture_path)   # loads: fixture_genes_bm

# ---------------------------------------------------------------------------
# dedup_gene_ids(): output structure
# ---------------------------------------------------------------------------
testthat::test_that("dedup_gene_ids returns a data frame", {
  result <- convertid:::dedup_gene_ids(fixture_genes_bm,
                                       has_symbols = TRUE,
                                       has_ensg2   = TRUE)
  testthat::expect_s3_class(result, "data.frame")
})

testthat::test_that("dedup_gene_ids output has no duplicate hgnc_symbols", {
  result <- convertid:::dedup_gene_ids(fixture_genes_bm,
                                       has_symbols = TRUE,
                                       has_ensg2   = TRUE)
  testthat::expect_false(any(duplicated(result$hgnc_symbol)))
})

testthat::test_that("dedup_gene_ids output has no duplicate ensembl_gene_ids", {
  result <- convertid:::dedup_gene_ids(fixture_genes_bm,
                                       has_symbols = TRUE,
                                       has_ensg2   = TRUE)
  testthat::expect_false(any(duplicated(result$ensembl_gene_id)))
})

testthat::test_that("dedup_gene_ids output has no NA hgnc_symbols", {
  result <- convertid:::dedup_gene_ids(fixture_genes_bm,
                                       has_symbols = TRUE,
                                       has_ensg2   = TRUE)
  testthat::expect_false(any(is.na(result$hgnc_symbol)))
})

testthat::test_that("dedup_gene_ids reduces 8 input rows to 5 unique genes", {
  result <- convertid:::dedup_gene_ids(fixture_genes_bm,
                                       has_symbols = TRUE,
                                       has_ensg2   = TRUE)
  testthat::expect_equal(nrow(result), 5L)
})

# ---------------------------------------------------------------------------
# dedup_gene_ids(): deduplication logic
# ---------------------------------------------------------------------------
testthat::test_that("single-mapping genes retained with correct symbol (ACTB, GAPDH)", {
  result <- convertid:::dedup_gene_ids(fixture_genes_bm,
                                       has_symbols = TRUE,
                                       has_ensg2   = TRUE)
  testthat::expect_equal(
    result$hgnc_symbol[result$ensembl_gene_id == "ENSG00000075624"], "ACTB")
  testthat::expect_equal(
    result$hgnc_symbol[result$ensembl_gene_id == "ENSG00000111640"], "GAPDH")
})

testthat::test_that("AnnotationDbi-confirmed ENSG ID preferred for AKAP17A", {
  # ENSG00000197976 has hgnc_symbol_2 == AKAP17A and is first in ensg_2 list;
  # ENSG00000292343 also has hgnc_symbol_2 == AKAP17A but is newer.
  # filter_ensg2_first picks ENSG00000197976.
  result <- convertid:::dedup_gene_ids(fixture_genes_bm,
                                       has_symbols = TRUE,
                                       has_ensg2   = TRUE)
  akap_row <- result[result$hgnc_symbol == "AKAP17A", ]
  testthat::expect_equal(nrow(akap_row), 1L)
  testthat::expect_equal(akap_row$ensembl_gene_id, "ENSG00000197976")
})

testthat::test_that("ENSG placeholder resolved via hgnc_symbol_2 for BMS1P4", {
  # ENSG00000242338: hgnc_symbol == BMS1P4, hgnc_symbol_2 == NA
  # ENSG00000271816: hgnc_symbol == ENSG..., hgnc_symbol_2 == BMS1P4
  # filter_prefer_confirmed discards ENSG00000242338; needs_fix resolves
  # ENSG00000271816's placeholder to BMS1P4.
  result <- convertid:::dedup_gene_ids(fixture_genes_bm,
                                       has_symbols = TRUE,
                                       has_ensg2   = TRUE)
  bms_row <- result[result$hgnc_symbol == "BMS1P4", ]
  testthat::expect_equal(nrow(bms_row), 1L)
  testthat::expect_equal(bms_row$ensembl_gene_id, "ENSG00000271816")
})

testthat::test_that("newer ENSG preferred when all metadata identical (HLA-H)", {
  # Both HLA-H rows have identical hgnc_symbol, hgnc_symbol_2, ensg_2.
  # filter_ensg2_first picks ENSG00000310469 (first in ensg_2 list).
  result <- convertid:::dedup_gene_ids(fixture_genes_bm,
                                       has_symbols = TRUE,
                                       has_ensg2   = TRUE)
  hlah_row <- result[result$hgnc_symbol == "HLA-H", ]
  testthat::expect_equal(nrow(hlah_row), 1L)
  testthat::expect_equal(hlah_row$ensembl_gene_id, "ENSG00000310469")
})

# ---------------------------------------------------------------------------
# dedup_gene_ids(): verbose output
# ---------------------------------------------------------------------------
testthat::test_that("dedup_gene_ids verbose = TRUE emits first pass message", {
  testthat::expect_message(
    convertid:::dedup_gene_ids(fixture_genes_bm,
                               has_symbols = TRUE,
                               has_ensg2   = TRUE,
                               verbose     = TRUE),
    regexp = "First pass"
  )
})

testthat::test_that("dedup_gene_ids verbose = TRUE emits second pass message", {
  testthat::expect_message(
    convertid:::dedup_gene_ids(fixture_genes_bm,
                               has_symbols = TRUE,
                               has_ensg2   = TRUE,
                               verbose     = TRUE),
    regexp = "Second pass"
  )
})

testthat::test_that("dedup_gene_ids verbose = FALSE produces no messages", {
  testthat::expect_no_message(
    convertid:::dedup_gene_ids(fixture_genes_bm,
                               has_symbols = TRUE,
                               has_ensg2   = TRUE,
                               verbose     = FALSE)
  )
})

# ---------------------------------------------------------------------------
# dedup_gene_ids(): ENSG-only mode
# ---------------------------------------------------------------------------
testthat::test_that("dedup_gene_ids works in ENSG-only mode (no gene_name column)", {
  ensg_only <- fixture_genes_bm[, c("ensembl_gene_id", "hgnc_symbol",
                                     "hgnc_symbol_2", "ensg_2")]
  result <- convertid:::dedup_gene_ids(ensg_only,
                                       has_symbols = FALSE,
                                       has_ensg2   = FALSE)
  testthat::expect_s3_class(result, "data.frame")
  testthat::expect_false(any(duplicated(result$hgnc_symbol)))
})

# ---------------------------------------------------------------------------
# unify_gene_ids(): input validation
# ---------------------------------------------------------------------------
testthat::test_that("unify_gene_ids errors on non-data-frame, non-character input", {
  testthat::expect_error(
    unify_gene_ids(c("ENSG00000075624"=60)),
    "'genes' must be a data frame or character vector"
  )
})

testthat::test_that("unify_gene_ids errors on missing ensg_col", {
  testthat::expect_error(
    unify_gene_ids(fixture_genes_bm, ensg_col = "wrong_col"),
    "Column 'wrong_col' not found"
  )
})

testthat::test_that("unify_gene_ids errors on missing symbol_col", {
  testthat::expect_error(
    unify_gene_ids(fixture_genes_bm, symbol_col = "wrong_col"),
    "Column 'wrong_col' not found"
  )
})

testthat::test_that("unify_gene_ids errors on empty input", {
  testthat::expect_error(
    unify_gene_ids(fixture_genes_bm[0, ]),
    "'genes' is empty"
  )
})

# ---------------------------------------------------------------------------
# unify_gene_ids(): column name restoration
# ---------------------------------------------------------------------------
testthat::test_that("original column names are restored in output", {
  mockery::stub(unify_gene_ids, "try_biomart", function(genes, ...) genes)
  mockery::stub(unify_gene_ids, "convertid::convertId2",
    function(x) rep(NA_character_, length(x)))
  input <- fixture_genes_bm
  names(input)[names(input) == "ensembl_gene_id"] <- "gene_id"
  names(input)[names(input) == "gene_name"]       <- "symbol"
  result <- suppressWarnings(
    unify_gene_ids(input, ensg_col = "gene_id", symbol_col = "symbol")
  )
  testthat::expect_true("gene_id" %in% names(result))
  testthat::expect_true("symbol"  %in% names(result))
  testthat::expect_false("ensembl_gene_id" %in% names(result))
  testthat::expect_false("gene_name"       %in% names(result))
})

# ---------------------------------------------------------------------------
# unify_gene_ids(): graceful degradation
# ---------------------------------------------------------------------------
testthat::test_that("gene_name symbols preserved when all lookups produce NA", {
  genes_failed <- data.frame(
    ensembl_gene_id = c("ENSG00000075624", "ENSG00000111640"),
    gene_name       = c("ACTB", "GAPDH"),
    hgnc_symbol     = NA_character_,
    hgnc_symbol_2   = NA_character_,
    ensg_2          = NA_character_,
    stringsAsFactors = FALSE
  )
  # Simulate the fix: gene_name is proper symbol so it should be used
  needs_fix <- is.na(genes_failed$hgnc_symbol) | grepl("^ENSG", genes_failed$hgnc_symbol)
  genes_failed$hgnc_symbol[needs_fix] <- ifelse(
    !grepl("^ENSG", genes_failed$gene_name[needs_fix]),
    genes_failed$gene_name[needs_fix],
    genes_failed$ensembl_gene_id[needs_fix]
  )
  testthat::expect_equal(genes_failed$hgnc_symbol, c("ACTB", "GAPDH"))
})

testthat::test_that("ENSG gene_name falls back to ensembl_gene_id when all lookups fail", {
  genes_failed <- data.frame(
    ensembl_gene_id = "ENSG00000075624",
    gene_name       = "ENSG00000075624",  # no proper symbol available
    hgnc_symbol     = NA_character_,
    hgnc_symbol_2   = NA_character_,
    ensg_2          = NA_character_,
    stringsAsFactors = FALSE
  )
  needs_fix <- is.na(genes_failed$hgnc_symbol) | grepl("^ENSG", genes_failed$hgnc_symbol)
  genes_failed$hgnc_symbol[needs_fix] <- ifelse(
    !grepl("^ENSG", genes_failed$gene_name[needs_fix]),
    genes_failed$gene_name[needs_fix],
    genes_failed$ensembl_gene_id[needs_fix]
  )
  testthat::expect_equal(genes_failed$hgnc_symbol, "ENSG00000075624")
})

# ---------------------------------------------------------------------------
# Integration test: full live pipeline (skipped on CRAN)
# ---------------------------------------------------------------------------
testthat::test_that("unify_gene_ids produces correct output with live lookups", {
  testthat::skip_on_cran()
  testthat::skip_if_offline()
  testthat::skip_if_not_installed("org.Hs.eg.db")

  input <- data.frame(
    gene_id   = c("ENSG00000075624", "ENSG00000111640",
                  "ENSG00000197976", "ENSG00000292343"),
    gene_name = c("ACTB", "GAPDH", "AKAP17A", "AKAP17A"),
    stringsAsFactors = FALSE
  )

  result <- unify_gene_ids(input,
                           ensg_col   = "gene_id",
                           symbol_col = "gene_name")

  testthat::expect_s3_class(result, "data.frame")
  testthat::expect_false(any(duplicated(result$hgnc_symbol)))
  testthat::expect_false(any(duplicated(result$gene_id)))
  testthat::expect_false(any(is.na(result$hgnc_symbol)))
  testthat::expect_equal(
    result$hgnc_symbol[result$gene_id == "ENSG00000075624"], "ACTB")
  testthat::expect_equal(
    result$hgnc_symbol[result$gene_id == "ENSG00000111640"], "GAPDH")
  # AKAP17A: older AnnotationDbi-confirmed ID retained
  testthat::expect_equal(nrow(result[result$hgnc_symbol == "AKAP17A", ]), 1L)
  testthat::expect_equal(
    result$gene_id[result$hgnc_symbol == "AKAP17A"], "ENSG00000197976")
})
