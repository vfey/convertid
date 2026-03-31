# Tests for convert.bm()
#
# convert.bm() is a wrapper around get.bm()/getBM() that merges BioMart results
# with the input data frame. Key behaviours tested here:
#
#   - Input validation and argument handling
#   - IDs not found in BioMart are retained via all.y=TRUE merge (with NAs)
#     and then have their sym.col filled with the filter ID
#   - rm.dups removes rows with duplicated biom.filter values
#   - biom.filter is prepended to biom.attributes if missing
#   - row names are used when id = "row.names"
#   - verbose produces messages
#
# Network-dependent behaviour is tested in the integration test at the bottom.
# All unit tests mock get.bm() to avoid network calls.

# ---------------------------------------------------------------------------
# Helper: a minimal mock get.bm() result for a known set of IDs
# ---------------------------------------------------------------------------
mock_get.bm <- function(values, ...) {
  # Returns BioMart-like result for two IDs; third ID ("ENSG_MISSING") absent
  data.frame(
    ensembl_gene_id = c("ENSG00000075624", "ENSG00000111640"),
    hgnc_symbol     = c("ACTB", "GAPDH"),
    description     = c("actin beta", "glyceraldehyde-3-phosphate dehydrogenase"),
    stringsAsFactors = FALSE
  )
}

mock_get.bm_dup <- function(values, ...) {
  # Returns a result with a duplicated ensembl_gene_id (as BioMart can do
  # when multiple hgnc_symbols map to one ENSG)
  data.frame(
    ensembl_gene_id = c("ENSG00000075624", "ENSG00000075624", "ENSG00000111640"),
    hgnc_symbol     = c("ACTB", "ACTB_DUP", "GAPDH"),
    description     = c("actin beta", "actin beta dup", "GAPDH"),
    stringsAsFactors = FALSE
  )
}

# ---------------------------------------------------------------------------
# Input structure
# ---------------------------------------------------------------------------
testthat::test_that("convert.bm accepts a data frame with an ID column", {
  mockery::stub(convert.bm, "get.bm", mock_get.bm)
  dat <- data.frame(
    ID = c("ENSG00000075624", "ENSG00000111640"),
    stringsAsFactors = FALSE
  )
  result <- convert.bm(dat, id = "ID")
  testthat::expect_s3_class(result, "data.frame")
  testthat::expect_true("hgnc_symbol" %in% names(result))
})

testthat::test_that("convert.bm accepts row names via id = 'row.names'", {
  mockery::stub(convert.bm, "get.bm", mock_get.bm)
  dat <- data.frame(
    some_col = c(1L, 2L),
    row.names = c("ENSG00000075624", "ENSG00000111640")
  )
  result <- convert.bm(dat, id = "row.names")
  testthat::expect_s3_class(result, "data.frame")
  testthat::expect_true("hgnc_symbol" %in% names(result))
})

testthat::test_that("convert.bm prepends biom.filter to biom.attributes if absent", {
  captured_attrs <- NULL
  mock_capture <- function(values, biom.data.set, biom.mart, host,
                           biom.filter, biom.attributes, ...) {
    captured_attrs <<- biom.attributes
    mock_get.bm(values)
  }
  mockery::stub(convert.bm, "get.bm", mock_capture)
  dat <- data.frame(ID = "ENSG00000075624", stringsAsFactors = FALSE)
  # biom.attributes does NOT include the filter
  convert.bm(dat, id = "ID",
              biom.filter     = "ensembl_gene_id",
              biom.attributes = c("hgnc_symbol", "description"))
  testthat::expect_true("ensembl_gene_id" %in% captured_attrs)
  testthat::expect_equal(captured_attrs[1], "ensembl_gene_id")
})

# ---------------------------------------------------------------------------
# Missing IDs: retained with sym.col filled by biom.filter value
# ---------------------------------------------------------------------------
testthat::test_that("IDs absent from BioMart result are retained with ENSG as symbol", {
  mockery::stub(convert.bm, "get.bm", mock_get.bm)
  dat <- data.frame(
    ID = c("ENSG00000075624", "ENSG00000111640", "ENSG_MISSING"),
    stringsAsFactors = FALSE
  )
  result <- convert.bm(dat, id = "ID")
  # All three IDs retained
  testthat::expect_equal(nrow(result), 3L)
  # Missing ID gets its ENSG as hgnc_symbol
  missing_row <- result[result$ensembl_gene_id == "ENSG_MISSING", ]
  testthat::expect_equal(missing_row$hgnc_symbol, "ENSG_MISSING")
})

testthat::test_that("empty hgnc_symbol is replaced by ensembl_gene_id", {
  mock_empty_sym <- function(values, ...) {
    data.frame(
      ensembl_gene_id = "ENSG00000075624",
      hgnc_symbol     = "",
      description     = "some gene",
      stringsAsFactors = FALSE
    )
  }
  mockery::stub(convert.bm, "get.bm", mock_empty_sym)
  dat <- data.frame(ID = "ENSG00000075624", stringsAsFactors = FALSE)
  result <- convert.bm(dat, id = "ID")
  testthat::expect_equal(result$hgnc_symbol, "ENSG00000075624")
})

# ---------------------------------------------------------------------------
# rm.dups
# ---------------------------------------------------------------------------
testthat::test_that("rm.dups = FALSE retains duplicated biom.filter rows", {
  mockery::stub(convert.bm, "get.bm", mock_get.bm_dup)
  dat <- data.frame(
    ID = c("ENSG00000075624", "ENSG00000111640"),
    stringsAsFactors = FALSE
  )
  result <- convert.bm(dat, id = "ID", rm.dups = FALSE)
  testthat::expect_equal(nrow(result), 3L)
})

testthat::test_that("rm.dups = TRUE removes duplicated biom.filter rows", {
  mockery::stub(convert.bm, "get.bm", mock_get.bm_dup)
  dat <- data.frame(
    ID = c("ENSG00000075624", "ENSG00000111640"),
    stringsAsFactors = FALSE
  )
  result <- convert.bm(dat, id = "ID", rm.dups = TRUE)
  testthat::expect_equal(nrow(result), 2L)
  testthat::expect_false(any(duplicated(result$ensembl_gene_id)))
})

# ---------------------------------------------------------------------------
# verbose
# ---------------------------------------------------------------------------
testthat::test_that("verbose = TRUE emits messages when replacements occur", {
  mockery::stub(convert.bm, "get.bm", mock_get.bm)
  dat <- data.frame(
    ID = c("ENSG00000075624", "ENSG00000111640", "ENSG_MISSING"),
    stringsAsFactors = FALSE
  )
  testthat::expect_message(
    convert.bm(dat, id = "ID", verbose = TRUE),
    regexp = "missing"
  )
})

testthat::test_that("verbose = TRUE emits message when rm.dups removes rows", {
  mockery::stub(convert.bm, "get.bm", mock_get.bm_dup)
  dat <- data.frame(
    ID = c("ENSG00000075624", "ENSG00000111640"),
    stringsAsFactors = FALSE
  )
  testthat::expect_message(
    convert.bm(dat, id = "ID", rm.dups = TRUE, verbose = TRUE),
    regexp = "duplicated"
  )
})

testthat::test_that("verbose = FALSE produces no messages", {
  mockery::stub(convert.bm, "get.bm", mock_get.bm)
  dat <- data.frame(
    ID = c("ENSG00000075624", "ENSG00000111640", "ENSG_MISSING"),
    stringsAsFactors = FALSE
  )
  testthat::expect_no_message(
    convert.bm(dat, id = "ID", verbose = FALSE)
  )
})

# ---------------------------------------------------------------------------
# Output structure
# ---------------------------------------------------------------------------
testthat::test_that("result contains all input columns plus BioMart attributes", {
  mockery::stub(convert.bm, "get.bm", mock_get.bm)
  dat <- data.frame(
    ID       = c("ENSG00000075624", "ENSG00000111640"),
    extra_col = c("a", "b"),
    stringsAsFactors = FALSE
  )
  result <- convert.bm(dat, id = "ID")
  testthat::expect_true("extra_col"       %in% names(result))
  testthat::expect_true("hgnc_symbol"     %in% names(result))
  testthat::expect_true("description"     %in% names(result))
  testthat::expect_true("ensembl_gene_id" %in% names(result))
})

testthat::test_that("result row count equals input row count", {
  mockery::stub(convert.bm, "get.bm", mock_get.bm)
  dat <- data.frame(
    ID = c("ENSG00000075624", "ENSG00000111640", "ENSG_MISSING"),
    stringsAsFactors = FALSE
  )
  result <- convert.bm(dat, id = "ID")
  testthat::expect_equal(nrow(result), nrow(dat))
})

# ---------------------------------------------------------------------------
# Integration test: real BioMart call (skipped on CRAN)
# ---------------------------------------------------------------------------
testthat::test_that("convert.bm returns correct symbols with live BioMart", {
  testthat::skip_on_cran()
  testthat::skip_if_offline()

  dat <- data.frame(
    ID = c("ENSG00000075624", "ENSG00000111640"),
    stringsAsFactors = FALSE
  )
  result <- convert.bm(dat, id = "ID")

  testthat::expect_s3_class(result, "data.frame")
  testthat::expect_equal(nrow(result), 2L)
  testthat::expect_equal(
    result$hgnc_symbol[result$ensembl_gene_id == "ENSG00000075624"], "ACTB")
  testthat::expect_equal(
    result$hgnc_symbol[result$ensembl_gene_id == "ENSG00000111640"], "GAPDH")
})
