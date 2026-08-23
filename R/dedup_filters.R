# =============================================================================
# Row-selection filters used by .dedup_gene_ids()
#
# Each filter takes one group of candidate rows (a data frame sharing a
# gene_name or an hgnc_symbol) and returns the subset that survives its
# criterion. .apply_filters() runs them in order and stops at the first one
# that narrows the group to exactly one row, so the order of the chain encodes
# the reconciliation priority; see .dedup_gene_ids() for the chains themselves.
#
# The filters are package-level rather than closures inside .dedup_gene_ids()
# so that each can be tested and stubbed on its own.
# =============================================================================

#' Unexported functions
#' Apply a chain of row filters until one resolves to a single row
#' @description \command{.apply_filters()} calls each filter in \code{filters}
#' on \code{x} in turn and returns the first result that contains exactly one
#' row. The chain therefore encodes a priority order: earlier filters win.
#' @param x (\code{data.frame}). A group of candidate rows.
#' @param filters (\code{list}). Filter functions, each taking \code{x} and
#' returning a subset of its rows.
#' @param force_single (\code{logical}). What to do when no filter resolves to a
#' single row: \code{TRUE} falls back to the first row of \code{x},
#' \code{FALSE} returns \code{x} unchanged. Defaults to \code{FALSE}.
#' @return A \code{data.frame}: the first single-row filter result, or
#' \code{x[1L, ]} when \code{force_single = TRUE}, or \code{x} unchanged.
#' @keywords internal
.apply_filters <- function(x, filters, force_single = FALSE) {
  for (f in filters) {
    candidate <- f(x)
    if (nrow(candidate) == 1L) return(candidate)
  }
  if (force_single) x[1L, ] else x
}

#' Unexported functions
#' Test whether an Ensembl gene ID occurs in a '///'-separated list
#' @description \command{.ensg_in_list()} checks, element by element, whether
#' \code{ensg_id} appears among the IDs held in the corresponding
#' \code{ensg_2} entry, which \command{convertId2()} returns as a
#' \code{" /// "}-separated string for one-to-many mappings.
#' @param ensg_id (\code{character}). Ensembl gene IDs.
#' @param ensg_2 (\code{character}). Matching \code{///}-separated ID lists.
#' @return A \code{logical} vector, \code{FALSE} where \code{ensg_2} is
#' \code{NA}.
#' @keywords internal
.ensg_in_list <- function(ensg_id, ensg_2) {
  !is.na(ensg_2) & mapply(function(id, e2) {
    id %in% trimws(strsplit(e2, "///")[[1L]])
  }, ensg_id, ensg_2)
}

#' Unexported functions
#' Prefer rows that AnnotationDbi confirmed
#' @description \command{.filter_prefer_confirmed()} discards rows with a
#' missing \code{hgnc_symbol_2} whenever a sibling row in the same group does
#' carry an AnnotationDbi symbol. Without this pre-filter an unconfirmed BioMart
#' row could be chosen over a confirmed one merely because its
#' \code{hgnc_symbol} happens to equal \code{gene_name}. It fires both when the
#' confirmed sibling matches \code{gene_name} and when it belongs to a different
#' symbol entirely.
#' @param x (\code{data.frame}). A group of candidate rows.
#' @return A \code{data.frame} with the unconfirmed rows dropped, or \code{x}
#' unchanged when no row is confirmed.
#' @keywords internal
.filter_prefer_confirmed <- function(x) {
  # Case 1: confirmed sibling has hgnc_symbol_2 == gene_name
  if (any(!is.na(x$hgnc_symbol_2) &
          x$hgnc_symbol_2 == x$gene_name &
          !grepl("^ENSG", x$hgnc_symbol_2)))
    return(x[!is.na(x$hgnc_symbol_2), ])
  # Case 2: confirmed sibling belongs to a different symbol entirely
  # (e.g. ENSG00000202198 has hgnc_symbol_2 == RN7SK, not 7SK)
  # Discard unconfirmed NA siblings, keeping only the confirmed row(s)
  if (any(!is.na(x$hgnc_symbol_2) &
          !grepl("^ENSG", x$hgnc_symbol_2) &
          x$hgnc_symbol_2 != x$gene_name))
    return(x[!is.na(x$hgnc_symbol_2), ])
  x
}

#' Unexported functions
#' Keep rows whose BioMart symbol matches the gene name
#' @description \command{.filter_symbol_matches_name()} keeps rows where
#' \code{hgnc_symbol} equals \code{gene_name} and is not a raw ENSG
#' placeholder. Requires a \code{gene_name} column.
#' @param x (\code{data.frame}). A group of candidate rows.
#' @return A \code{data.frame} with the matching rows.
#' @keywords internal
.filter_symbol_matches_name <- function(x)
  x[!is.na(x$hgnc_symbol) &
      x$hgnc_symbol == x$gene_name &
      !grepl("^ENSG", x$hgnc_symbol), ]

#' Unexported functions
#' Keep rows whose AnnotationDbi symbol matches the gene name
#' @description \command{.filter_symbol2_matches_name()} keeps rows where
#' \code{hgnc_symbol_2} equals \code{gene_name} and is not a raw ENSG
#' placeholder. Requires a \code{gene_name} column.
#' @param x (\code{data.frame}). A group of candidate rows.
#' @return A \code{data.frame} with the matching rows.
#' @keywords internal
.filter_symbol2_matches_name <- function(x)
  x[!is.na(x$hgnc_symbol_2) &
      x$hgnc_symbol_2 == x$gene_name &
      !grepl("^ENSG", x$hgnc_symbol_2), ]

#' Unexported functions
#' Keep rows whose symbol equals the gene name
#' @description \command{.filter_symbol_matches_gene_name()} is the plain
#' equality form of \command{.filter_symbol_matches_name()}, without the ENSG
#' placeholder test. It is used in the second deduplication pass, where the
#' placeholders have already been resolved.
#' @param x (\code{data.frame}). A group of candidate rows.
#' @return A \code{data.frame} with the matching rows.
#' @keywords internal
.filter_symbol_matches_gene_name <- function(x)
  x[x$hgnc_symbol == x$gene_name, ]

#' Unexported functions
#' Keep rows confirmed by the first AnnotationDbi Ensembl ID
#' @description \command{.filter_ensg2_first()} keeps rows whose
#' \code{ensembl_gene_id} equals the first entry of the \code{///}-separated
#' \code{ensg_2} list. The ordering of that list is not a strong preference
#' signal on its own, which is why this filter sits after the source-agreement
#' filters in the chains.
#' @param x (\code{data.frame}). A group of candidate rows.
#' @return A \code{data.frame} with the confirmed rows.
#' @keywords internal
.filter_ensg2_first <- function(x) {
  first_ensg <- trimws(sapply(strsplit(as.character(x$ensg_2), "///"), `[`, 1L))
  x[!is.na(x$ensg_2) & x$ensembl_gene_id == first_ensg, ]
}

#' Unexported functions
#' Drop rows whose Ensembl ID appears in the AnnotationDbi list
#' @description \command{.filter_ensg_not_in_list()} keeps rows whose
#' \code{ensembl_gene_id} is absent from their \code{ensg_2} list, preferring
#' the more canonical ID when several map to one symbol.
#' @param x (\code{data.frame}). A group of candidate rows.
#' @return A \code{data.frame} with the remaining rows.
#' @seealso \code{\link{.ensg_in_list}}
#' @keywords internal
.filter_ensg_not_in_list <- function(x)
  x[!.ensg_in_list(x$ensembl_gene_id, x$ensg_2), ]

#' Unexported functions
#' Keep rows where both annotation sources agree
#' @description \command{.filter_symbols_agree()} keeps rows where the BioMart
#' symbol and the AnnotationDbi symbol are identical, i.e. cross-source
#' confirmation. It needs neither \code{gene_name} nor \code{ensg_2} and so
#' appears in every filter chain.
#' @param x (\code{data.frame}). A group of candidate rows.
#' @return A \code{data.frame} with the agreeing rows.
#' @keywords internal
.filter_symbols_agree <- function(x)
  x[!is.na(x$hgnc_symbol_2) & x$hgnc_symbol == x$hgnc_symbol_2, ]

#' Unexported functions
#' Drop rows whose symbol is still a raw Ensembl ID
#' @description \command{.filter_drop_ensg_symbol()} deprioritises rows whose
#' \code{hgnc_symbol} is an unresolved ENSG placeholder rather than a symbol.
#' @param x (\code{data.frame}). A group of candidate rows.
#' @return A \code{data.frame} with the placeholder rows removed.
#' @keywords internal
.filter_drop_ensg_symbol <- function(x)
  x[!grepl("^ENSG", x$hgnc_symbol), ]

#' Unexported functions
#' Keep rows carrying an AnnotationDbi symbol
#' @description \command{.filter_has_symbol2()} keeps rows with a non-missing
#' \code{hgnc_symbol_2}, used as a late tiebreaker.
#' @param x (\code{data.frame}). A group of candidate rows.
#' @return A \code{data.frame} with the rows that have an AnnotationDbi symbol.
#' @keywords internal
.filter_has_symbol2 <- function(x)
  x[!is.na(x$hgnc_symbol_2), ]

#' Unexported functions
#' Keep rows carrying an AnnotationDbi Ensembl ID
#' @description \command{.filter_has_ensg2()} keeps rows with a non-missing
#' \code{ensg_2}, used as the last tiebreaker before falling back to the first
#' row.
#' @param x (\code{data.frame}). A group of candidate rows.
#' @return A \code{data.frame} with the rows that have an AnnotationDbi ID.
#' @keywords internal
.filter_has_ensg2 <- function(x)
  x[!is.na(x$ensg_2), ]

#' Unexported functions
#' Take the first row when nothing can disambiguate the group
#' @description \command{.filter_last_resort()} returns the first row, but only
#' when every disambiguation field (\code{ensg_2} and \code{hgnc_symbol_2}) is
#' \code{NA} across the whole group. Otherwise the group is returned unchanged
#' so that later filters still get their chance.
#' @param x (\code{data.frame}). A group of candidate rows.
#' @return A \code{data.frame}: \code{x[1L, ]} when the group carries no
#' disambiguating information, otherwise \code{x}.
#' @keywords internal
.filter_last_resort <- function(x) {
  if (all(is.na(x$ensg_2)) && all(is.na(x$hgnc_symbol_2))) x[1L, ] else x
}
