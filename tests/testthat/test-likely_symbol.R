# Tests for likely_symbol()
#
# Key behaviours tested:
#   - Character vector and data frame input
#   - output = "likely", "symbols", "all" return correct columns
#   - Symbol resolution: current, alias, previous, unrecognised
#   - alias_sym and prev_sym flags
#   - index_threshold: correct path selection, identical results from both paths
#   - Session cache: reused on second call, bypassed when hgnc supplied
#   - refresh = TRUE forces re-download and updates cache timestamp
#   - Stale cache (> 3 days) emits a warning message
#   - Unsupported organism returns input unchanged
#   - verbose = TRUE / FALSE
#
# All unit tests supply a minimal hgnc fixture or manipulate the package cache
# directly to avoid any network calls.

# ---------------------------------------------------------------------------
# Minimal HGNC fixture
#   ACTB  - current symbol, no aliases, no previous
#   KYAT1 - alias CCBL1, previous symbol KAAT1
#   BRCA1 - alias RNF53
# ---------------------------------------------------------------------------
hgnc_fixture <- data.frame(
  symbol        = c("ACTB",  "KYAT1", "BRCA1"),
  alias_symbol  = c("",      "CCBL1", "RNF53"),
  prev_symbol   = c("",      "KAAT1", ""),
  stringsAsFactors = FALSE
)

# Helper: clear the package cache before tests that depend on cache state
clear_cache <- function() {
  if (exists("hgnc", envir = .likely_symbol_cache, inherits = FALSE))
    rm("hgnc", "downloaded_at", envir = .likely_symbol_cache)
}

# ---------------------------------------------------------------------------
# Input types
# ---------------------------------------------------------------------------
testthat::test_that("likely_symbol accepts a character vector", {
  result <- likely_symbol("ACTB", hgnc = hgnc_fixture, verbose = FALSE)
  testthat::expect_s3_class(result, "data.frame")
  testthat::expect_true("input_symbol"  %in% names(result))
  testthat::expect_true("likely_symbol" %in% names(result))
})

testthat::test_that("likely_symbol accepts a data frame with Material, org, sp columns", {
  dat <- data.frame(Material = 1L, org = "human", sp = "ACTB",
                    stringsAsFactors = FALSE)
  result <- likely_symbol(dat, hgnc = hgnc_fixture, verbose = FALSE)
  testthat::expect_s3_class(result, "data.frame")
})

testthat::test_that("likely_symbol errors on data frame missing required columns", {
  dat <- data.frame(wrong = "ACTB", stringsAsFactors = FALSE)
  testthat::expect_error(
    likely_symbol(dat, hgnc = hgnc_fixture, verbose = FALSE),
    "Need the following columns"
  )
})

# ---------------------------------------------------------------------------
# output options
# ---------------------------------------------------------------------------
testthat::test_that("output = 'likely' returns likely_symbol and input_symbol only", {
  result <- likely_symbol("ACTB", hgnc = hgnc_fixture,
                          output = "likely", verbose = FALSE)
  testthat::expect_equal(sort(names(result)), sort(c("likely_symbol", "input_symbol")))
})

testthat::test_that("output = 'symbols' returns correct columns", {
  result <- likely_symbol("ACTB", hgnc = hgnc_fixture,
                          output = "symbols", verbose = FALSE)
  testthat::expect_equal(sort(names(result)),
    sort(c("current_symbols", "likely_symbol", "input_symbol", "all_symbols")))
})

testthat::test_that("output = 'all' returns correct columns", {
  result <- likely_symbol("ACTB", hgnc = hgnc_fixture,
                          output = "all", verbose = FALSE)
  testthat::expect_equal(sort(names(result)),
    sort(c("current_symbols", "likely_symbol", "previous_symbol",
           "input_symbol", "all_symbols")))
})

# ---------------------------------------------------------------------------
# Symbol resolution logic
# ---------------------------------------------------------------------------
testthat::test_that("current symbol resolves to itself", {
  result <- likely_symbol("ACTB", hgnc = hgnc_fixture,
                          output = "likely", verbose = FALSE)
  testthat::expect_equal(result$input_symbol,  "ACTB")
  testthat::expect_equal(result$likely_symbol, "ACTB")
})

testthat::test_that("alias resolves to current symbol", {
  result <- likely_symbol("CCBL1", hgnc = hgnc_fixture,
                          output = "likely", verbose = FALSE)
  testthat::expect_equal(result$input_symbol,  "CCBL1")
  testthat::expect_equal(result$likely_symbol, "KYAT1")
})

testthat::test_that("previous symbol resolves to current symbol", {
  result <- likely_symbol("KAAT1", hgnc = hgnc_fixture,
                          output = "likely", verbose = FALSE)
  testthat::expect_equal(result$input_symbol,  "KAAT1")
  testthat::expect_equal(result$likely_symbol, "KYAT1")
})

testthat::test_that("unrecognised symbol returns input unchanged", {
  result <- likely_symbol("NOTAREAL1", hgnc = hgnc_fixture,
                          output = "likely", verbose = FALSE)
  testthat::expect_equal(result$input_symbol,  "NOTAREAL1")
  testthat::expect_equal(result$likely_symbol, "NOTAREAL1")
})

testthat::test_that("multiple symbols are all resolved correctly", {
  result <- likely_symbol(c("ACTB", "CCBL1", "KAAT1"),
                          hgnc = hgnc_fixture, output = "likely", verbose = FALSE)
  testthat::expect_equal(nrow(result), 3L)
  testthat::expect_equal(result$likely_symbol[result$input_symbol == "CCBL1"], "KYAT1")
  testthat::expect_equal(result$likely_symbol[result$input_symbol == "KAAT1"], "KYAT1")
  testthat::expect_equal(result$likely_symbol[result$input_symbol == "ACTB"],  "ACTB")
})

# ---------------------------------------------------------------------------
# alias_sym and prev_sym flags
# ---------------------------------------------------------------------------
testthat::test_that("alias_sym = FALSE prevents alias lookup", {
  result <- likely_symbol("CCBL1", hgnc = hgnc_fixture,
                          alias_sym = FALSE, prev_sym = FALSE,
                          output = "likely", verbose = FALSE)
  testthat::expect_equal(result$likely_symbol, "CCBL1")
})

testthat::test_that("prev_sym = FALSE prevents previous symbol lookup", {
  result <- likely_symbol("KAAT1", hgnc = hgnc_fixture,
                          alias_sym = FALSE, prev_sym = FALSE,
                          output = "likely", verbose = FALSE)
  testthat::expect_equal(result$likely_symbol, "KAAT1")
})

# ---------------------------------------------------------------------------
# index_threshold
# ---------------------------------------------------------------------------
testthat::test_that("row-scan used below index_threshold (no index message)", {
  testthat::expect_no_message(
    likely_symbol("ACTB", hgnc = hgnc_fixture,
                  index_threshold = 10L, verbose = TRUE),
    message = "Building alias symbol index"
  )
})

testthat::test_that("index built at or above index_threshold", {
  testthat::expect_message(
    likely_symbol("ACTB", hgnc = hgnc_fixture,
                  index_threshold = 1L, verbose = TRUE),
    regexp = "Building alias symbol index"
  )
})

testthat::test_that("both paths produce identical results for alias lookup", {
  r_scan  <- likely_symbol("CCBL1", hgnc = hgnc_fixture,
                            index_threshold = 99L, output = "likely", verbose = FALSE)
  r_index <- likely_symbol("CCBL1", hgnc = hgnc_fixture,
                            index_threshold = 1L,  output = "likely", verbose = FALSE)
  testthat::expect_equal(r_scan, r_index)
})

testthat::test_that("both paths produce identical results for previous symbol lookup", {
  r_scan  <- likely_symbol("KAAT1", hgnc = hgnc_fixture,
                            index_threshold = 99L, output = "likely", verbose = FALSE)
  r_index <- likely_symbol("KAAT1", hgnc = hgnc_fixture,
                            index_threshold = 1L,  output = "likely", verbose = FALSE)
  testthat::expect_equal(r_scan, r_index)
})

testthat::test_that("both paths produce identical results for batch input", {
  syms    <- c("ACTB", "CCBL1", "KAAT1", "NOTAREAL1")
  r_scan  <- likely_symbol(syms, hgnc = hgnc_fixture,
                            index_threshold = 99L, output = "likely", verbose = FALSE)
  r_index <- likely_symbol(syms, hgnc = hgnc_fixture,
                            index_threshold = 1L,  output = "likely", verbose = FALSE)
  r_scan  <- r_scan[order(r_scan$input_symbol), ]
  r_index <- r_index[order(r_index$input_symbol), ]
  testthat::expect_equal(r_scan, r_index)
})

# ---------------------------------------------------------------------------
# Session cache
# ---------------------------------------------------------------------------
testthat::test_that("hgnc table is cached after first call", {
  clear_cache()
  mockery::stub(likely_symbol, "read.delim", function(...) hgnc_fixture)
  likely_symbol("ACTB", verbose = FALSE)
  testthat::expect_true(
    exists("hgnc", envir = .likely_symbol_cache, inherits = FALSE))
  testthat::expect_true(
    exists("downloaded_at", envir = .likely_symbol_cache, inherits = FALSE))
})

testthat::test_that("second call reuses cache without downloading", {
  clear_cache()
  # Seed the cache
  .likely_symbol_cache$hgnc          <- hgnc_fixture
  .likely_symbol_cache$downloaded_at <- Sys.time()
  # read.delim should never be called on second use
  mockery::stub(likely_symbol, "read.delim", function(...) {
    stop("Should not download when cache is available")
  })
  testthat::expect_no_error(
    likely_symbol("ACTB", verbose = FALSE)
  )
})

testthat::test_that("cache reuse emits message when verbose = TRUE", {
  clear_cache()
  .likely_symbol_cache$hgnc          <- hgnc_fixture
  .likely_symbol_cache$downloaded_at <- Sys.time()
  testthat::expect_message(
    likely_symbol("ACTB", verbose = TRUE),
    regexp = "Reusing cached HGNC table"
  )
})

testthat::test_that("stale cache (> 3 days) emits staleness warning", {
  clear_cache()
  .likely_symbol_cache$hgnc          <- hgnc_fixture
  .likely_symbol_cache$downloaded_at <- Sys.time() - (4 * 24 * 3600)
  testthat::expect_message(
    likely_symbol("ACTB", verbose = TRUE),
    regexp = "may be outdated"
  )
})

testthat::test_that("supplying hgnc bypasses cache entirely", {
  clear_cache()
  # Seed cache with a different table that would give wrong results
  wrong_hgnc <- data.frame(symbol = "WRONG", alias_symbol = "ACTB",
                            prev_symbol = "", stringsAsFactors = FALSE)
  .likely_symbol_cache$hgnc          <- wrong_hgnc
  .likely_symbol_cache$downloaded_at <- Sys.time()
  # Supplying hgnc directly should ignore the cache
  result <- likely_symbol("ACTB", hgnc = hgnc_fixture,
                           output = "likely", verbose = FALSE)
  testthat::expect_equal(result$likely_symbol, "ACTB")
})

testthat::test_that("supplying hgnc emits 'Using supplied' message when verbose", {
  testthat::expect_message(
    likely_symbol("ACTB", hgnc = hgnc_fixture, verbose = TRUE),
    regexp = "Using supplied HGNC table"
  )
})

# ---------------------------------------------------------------------------
# refresh = TRUE
# ---------------------------------------------------------------------------
testthat::test_that("refresh = TRUE triggers re-download even when cache exists", {
  clear_cache()
  .likely_symbol_cache$hgnc          <- hgnc_fixture
  .likely_symbol_cache$downloaded_at <- Sys.time()
  download_called <- FALSE
  mockery::stub(likely_symbol, "read.delim", function(...) {
    download_called <<- TRUE
    hgnc_fixture
  })
  likely_symbol("ACTB", refresh = TRUE, verbose = FALSE)
  testthat::expect_true(download_called)
})

testthat::test_that("refresh = TRUE updates the cache timestamp", {
  clear_cache()
  old_time <- Sys.time() - (2 * 24 * 3600)
  .likely_symbol_cache$hgnc          <- hgnc_fixture
  .likely_symbol_cache$downloaded_at <- old_time
  mockery::stub(likely_symbol, "read.delim", function(...) hgnc_fixture)
  likely_symbol("ACTB", refresh = TRUE, verbose = FALSE)
  testthat::expect_gt(
    as.numeric(.likely_symbol_cache$downloaded_at),
    as.numeric(old_time)
  )
})

testthat::test_that("refresh = TRUE emits refreshing message when verbose", {
  clear_cache()
  .likely_symbol_cache$hgnc          <- hgnc_fixture
  .likely_symbol_cache$downloaded_at <- Sys.time()
  mockery::stub(likely_symbol, "read.delim", function(...) hgnc_fixture)
  testthat::expect_message(
    likely_symbol("ACTB", refresh = TRUE, verbose = TRUE),
    regexp = "Refreshing HGNC table"
  )
})

# ---------------------------------------------------------------------------
# Organism handling
# ---------------------------------------------------------------------------
testthat::test_that("unsupported organism returns input unchanged", {
  result <- likely_symbol("Actb", orgnsm = "rat",
                          output = "likely", verbose = FALSE)
  testthat::expect_equal(result$input_symbol, "Actb")
})

testthat::test_that("NULL orgnsm is treated as empty string without error", {
  testthat::expect_no_error(
    likely_symbol("ACTB", orgnsm = NULL, hgnc = hgnc_fixture, verbose = FALSE)
  )
})

# ---------------------------------------------------------------------------
# verbose
# ---------------------------------------------------------------------------
testthat::test_that("verbose = TRUE emits symbol alias message", {
  testthat::expect_message(
    likely_symbol("ACTB", hgnc = hgnc_fixture, output = "likely", verbose = TRUE),
    regexp = "Getting symbol aliases"
  )
})

testthat::test_that("verbose = FALSE suppresses all messages", {
  testthat::expect_no_message(
    likely_symbol("ACTB", hgnc = hgnc_fixture, output = "likely", verbose = FALSE)
  )
})

# ---------------------------------------------------------------------------
# Integration tests: real HGNC download (skipped on CRAN)
# ---------------------------------------------------------------------------
testthat::test_that("likely_symbol resolves known alias with live HGNC download", {
  testthat::skip_on_cran()
  testthat::skip_if_offline()
  clear_cache()
  result <- likely_symbol("CCBL1", output = "likely", verbose = FALSE)
  testthat::expect_s3_class(result, "data.frame")
  testthat::expect_equal(result$input_symbol,  "CCBL1")
  testthat::expect_equal(result$likely_symbol, "KYAT1")
})

testthat::test_that("second live call reuses cache without re-downloading", {
  testthat::skip_on_cran()
  testthat::skip_if_offline()
  # First call populates cache
  likely_symbol("CCBL1", output = "likely", verbose = FALSE)
  # Second call should reuse cache — stub read.delim to confirm no download
  mockery::stub(likely_symbol, "read.delim", function(...) {
    stop("Should not download on second call")
  })
  testthat::expect_no_error(
    likely_symbol("KAAT1", output = "likely", verbose = FALSE)
  )
})

testthat::test_that("index and row-scan paths agree with live HGNC download", {
  testthat::skip_on_cran()
  testthat::skip_if_offline()
  syms    <- c("CCBL1", "KAAT1", "ACTB", "NOTAREAL1")
  r_scan  <- likely_symbol(syms, index_threshold = 99L,
                            output = "likely", verbose = FALSE)
  r_index <- likely_symbol(syms, index_threshold = 1L,
                            output = "likely", verbose = FALSE)
  r_scan  <- r_scan[order(r_scan$input_symbol), ]
  r_index <- r_index[order(r_index$input_symbol), ]
  testthat::expect_equal(r_scan, r_index)
})
