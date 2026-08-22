# Test helper: skip a test when BioMart is not reachable.
#
# A port-level connectivity check (testthat::skip_if_offline()) is not
# sufficient for BioMart: the host can be reachable while the service itself
# is down or the mart listing fails. This helper therefore attempts an actual
# useEnsembl() session and skips when that errors.
#
# Kept in tests/testthat/ rather than R/ so that the package code does not
# depend on testthat (a Suggests package) and does not open network
# connections outside of tests.

skip_if_biomart_unavailable <- function() {
  available <- tryCatch(
    {
      biomaRt::useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl")
      TRUE
    },
    error = function(e) FALSE
  )
  testthat::skip_if_not(available, "BioMart is not available")
}
