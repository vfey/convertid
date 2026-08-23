#' Unexported functions
#' Deduplicate reconciled gene ID rows
#' @description \command{.dedup_gene_ids()} takes a data frame in which the
#' BioMart and AnnotationDbi lookups have already been carried out, i.e. one
#' that carries populated \code{hgnc_symbol}, \code{hgnc_symbol_2} and
#' \code{ensg_2} columns, and reduces it to one row per gene with a unified
#' \code{hgnc_symbol}.
#'
#' Two passes are run. The first groups by \code{gene_name} (or by
#' \code{ensembl_gene_id} in ENSG-only mode) and resolves several Ensembl IDs
#' mapping to one gene; the second groups by \code{hgnc_symbol} and resolves
#' several gene names mapping to one symbol, always down to a single row. Each
#' pass runs a chain of filters via \command{.apply_filters()}; which chain is
#' used depends on whether gene symbols and AnnotationDbi Ensembl IDs are
#' available.
#' @param genes (\code{data.frame}). Columns \code{ensembl_gene_id},
#' \code{hgnc_symbol}, \code{hgnc_symbol_2}, \code{ensg_2}, and optionally
#' \code{gene_name}.
#' @param has_symbols (\code{logical}). Is a \code{gene_name} column present?
#' @param has_ensg2 (\code{logical}). Does \code{ensg_2} hold any non-\code{NA}
#' value?
#' @param verbose (\code{logical}). Should progress and summary messages be
#' written to the console? Defaults to \code{FALSE}.
#' @return A \code{data.frame} with one row per gene and unified
#' \code{hgnc_symbol} values.
#' @seealso \code{\link{unify_gene_ids}}, \code{\link{.apply_filters}}
#' @keywords internal
.dedup_gene_ids <- function(genes, has_symbols, has_ensg2, verbose = FALSE) {

  # ---------------------------------------------------------------------------
  # Build filter chains based on available data. The order of each chain is the
  # reconciliation priority; see the individual filters in dedup_filters.R.
  # ---------------------------------------------------------------------------
  if (has_symbols && has_ensg2) {
    # Full filter chain: gene symbols + ensg_2 available
    dedup_filters_by_name <- list(
      .filter_prefer_confirmed,          # 0.  Discard NA hgnc_symbol_2 when confirmed alternative exists
      .filter_symbol_matches_name,       # 1.  BioMart hgnc matches gene_name
      .filter_symbol2_matches_name,      # 2.  AnnotationDbi hgnc matches gene_name
      .filter_symbols_agree,             # 3.  Both sources agree on symbol
      .filter_ensg2_first,               # 4.  ENSG confirmed by first entry in ensg_2
      .filter_drop_ensg_symbol,          # 5.  Drop raw ENSG placeholders
      .filter_last_resort                # 6.  All NA: take first row
    )
    dedup_filters_by_symbol <- list(
      .filter_symbols_agree,             # 1. Both sources agree
      .filter_symbol_matches_gene_name,  # 2. hgnc_symbol matches gene_name
      .filter_ensg2_first,               # 3. ENSG confirmed by first entry in ensg_2
      .filter_ensg_not_in_list,          # 4. Prefer more canonical ENSG
      .filter_has_symbol2,               # 5. Has AnnotationDbi symbol
      .filter_has_ensg2                  # 6. Has AnnotationDbi ENSG
    )
  } else if (has_symbols) {
    # Gene symbols available but ensg_2 lookup failed or unavailable
    dedup_filters_by_name <- list(
      .filter_prefer_confirmed,          # 0.  Discard NA hgnc_symbol_2 when confirmed alternative exists
      .filter_symbol_matches_name,       # 1.  BioMart hgnc matches gene_name
      .filter_symbol2_matches_name,      # 2.  AnnotationDbi hgnc matches gene_name
      .filter_symbols_agree,             # 3.  Both sources agree on symbol
      .filter_drop_ensg_symbol,          # 4.  Drop raw ENSG placeholders
      .filter_last_resort                # 5.  All NA: take first row
    )
    dedup_filters_by_symbol <- list(
      .filter_symbols_agree,             # 1. Both sources agree
      .filter_symbol_matches_gene_name,  # 2. hgnc_symbol matches gene_name
      .filter_has_symbol2,               # 3. Has AnnotationDbi symbol
      .filter_has_ensg2                  # 4. Has AnnotationDbi ENSG
    )
  } else {
    # ENSG-only: no gene symbols, no ensg_2
    dedup_filters_by_name <- list(
      .filter_symbols_agree,             # 1. Both sources agree
      .filter_drop_ensg_symbol,          # 2. Drop raw ENSG placeholders
      .filter_has_symbol2,               # 3. Has AnnotationDbi symbol
      .filter_last_resort                # 4. All NA: take first row
    )
    dedup_filters_by_symbol <- list(
      .filter_symbols_agree,             # 1. Both sources agree
      .filter_ensg_not_in_list,          # 2. Prefer more canonical ENSG
      .filter_has_symbol2,               # 3. Has AnnotationDbi symbol
      .filter_has_ensg2                  # 4. Has AnnotationDbi ENSG
    )
  }

  # ---------------------------------------------------------------------------
  # First pass: deduplicate by gene_name (if available) or ensembl_gene_id
  # ---------------------------------------------------------------------------
  ddply_col <- if (has_symbols) "gene_name" else "ensembl_gene_id"
  if (verbose) message(sprintf("First pass: deduplicating by %s...", ddply_col))

  genes2_intermediate <- plyr::ddply(genes, ddply_col, function(x) {

    if (nrow(x) > 1L) {
      # Remove duplicate ensembl_gene_ids upfront, preferring hgnc_symbol == gene_name
      if (any(duplicated(x$ensembl_gene_id))) {
        if (has_symbols) {
          x <- x[order(x$ensembl_gene_id, x$hgnc_symbol != x$gene_name), ]
        } else {
          x <- x[order(x$ensembl_gene_id), ]
        }
        x <- x[!duplicated(x$ensembl_gene_id), ]
      }
      x <- .apply_filters(x, dedup_filters_by_name, force_single = FALSE)
    }

    # Fix hgnc_symbol if it is a raw ENSG placeholder. When gene_name is
    # available the replacement is applied only if gene_name is itself a proper
    # symbol (not an ENSG ID) and differs from hgnc_symbol. In ENSG-only mode
    # there is no gene_name column, so every placeholder qualifies -- note that
    # x$gene_name would be NULL there, which would silently collapse the
    # condition to length zero and disable the fix altogether.
    is_placeholder <- grepl("^ENSG", x$hgnc_symbol)
    if (has_symbols) {
      needs_fix <- is_placeholder &
                   !is.na(x$gene_name) &
                   !grepl("^ENSG", x$gene_name) &
                   x$gene_name != x$hgnc_symbol
    } else {
      needs_fix <- is_placeholder
    }
    # Guard against NA from any of the comparisons above: NA indices are not
    # allowed in subscripted assignment and would abort the whole call.
    needs_fix[is.na(needs_fix)] <- FALSE

    from_symbol2  <- needs_fix & !is.na(x$hgnc_symbol_2)
    from_fallback <- needs_fix &  is.na(x$hgnc_symbol_2)
    x$hgnc_symbol[from_symbol2] <- x$hgnc_symbol_2[from_symbol2]
    x$hgnc_symbol[from_fallback] <-
      if (has_symbols) x$gene_name[from_fallback]
      else             x$ensembl_gene_id[from_fallback]
    x
  })

  if (verbose) {
    n_inter <- nrow(genes2_intermediate)
    message(sprintf("First pass: %d -> %d rows (%d removed).",
      nrow(genes), n_inter, nrow(genes) - n_inter))
  }

  # ---------------------------------------------------------------------------
  # Second pass: deduplicate by hgnc_symbol, always resolving to 1 row
  # ---------------------------------------------------------------------------
  if (verbose) message("Second pass: deduplicating by hgnc_symbol...")

  genes2 <- plyr::ddply(genes2_intermediate, "hgnc_symbol", function(x) {
    if (nrow(x) == 1L) return(x)
    .apply_filters(x, dedup_filters_by_symbol, force_single = TRUE)
  })

  if (verbose) {
    n_final <- nrow(genes2)
    message(sprintf("Second pass: %d -> %d rows (%d removed).",
      nrow(genes2_intermediate), n_final,
      nrow(genes2_intermediate) - n_final))
    message(sprintf("Done. Final output: %d rows.", n_final))
    message(sprintf("  Duplicate hgnc_symbols:     %d", sum(duplicated(genes2$hgnc_symbol))))
    message(sprintf("  Duplicate ensembl_gene_ids: %d", sum(duplicated(genes2$ensembl_gene_id))))
    message(sprintf("  NA hgnc_symbols:            %d", sum(is.na(genes2$hgnc_symbol))))
  }

  genes2
}


# =============================================================================
#' Unify gene IDs from BioMart and AnnotationDbi lookups
#'
#' Takes a data frame with Ensembl gene IDs (and optionally gene symbols) and
#' returns a deduplicated data frame with unified HGNC symbols, using a
#' priority-based reconciliation of BioMart and AnnotationDbi results.
#'
#' @param genes A data frame with at minimum an Ensembl gene ID column or a
#'   character vector of Ensembl gene IDs.
#' @param ensg_col Name of the column containing Ensembl gene IDs.
#'   Default: \code{"ensembl_gene_id"}.
#' @param symbol_col Name of the column containing gene symbols, or \code{NULL}
#'   if absent (ENSG-only mode). Default: \code{NULL}.
#' @param host BioMart host URL. Default: \code{"https://www.ensembl.org"}.
#' @param biomart_fallback Character vector of fallback BioMart host URLs to try
#'   if the primary host fails. Set to \code{NULL} to disable fallback.
#' @param keep_intermediates Logical; if \code{TRUE}, the intermediate lookup
#'   columns \code{hgnc_symbol_2} and \code{ensg_2} are retained in the output.
#'   Useful for debugging. Default: \code{FALSE}.
#' @param verbose Logical; if \code{TRUE}, print progress and summary messages.
#'   Default: \code{FALSE}.
#'
#' @return A deduplicated data frame with unified HGNC symbols in the
#'   \code{hgnc_symbol} column. The Ensembl ID and gene symbol columns are
#'   returned under the names they were given in, i.e. those passed as
#'   \code{ensg_col} and \code{symbol_col}. The intermediate lookup columns
#'   \code{hgnc_symbol_2} and \code{ensg_2} are dropped unless
#'   \code{keep_intermediates = TRUE}.
#'
#' @details
#' Requires the Bioconductor packages \pkg{org.Hs.eg.db} and \pkg{AnnotationDbi}.
#' These are not hard dependencies but will be checked at runtime with an
#' informative error if missing.
#'
#' \strong{Deduplication passes}
#'
#' The function performs two sequential deduplication passes via the internal
#' \code{.dedup_gene_ids()} function:
#' \enumerate{
#'   \item Deduplicate by \code{gene_name} (if available) or \code{ensembl_gene_id},
#'     resolving multiple ENSG IDs mapping to the same gene name.
#'   \item Deduplicate by \code{hgnc_symbol}, resolving cases where multiple
#'     gene names resolve to the same symbol.
#' }
#'
#' \strong{Symbol assignment priority}
#'
#' The guiding principle is that \strong{AnnotationDbi confirmation outranks
#' BioMart ordering}. AnnotationDbi (\pkg{org.Hs.eg.db}) reflects a stable,
#' versioned annotation database, while BioMart returns the current Ensembl
#' release which may be ahead of annotations used to build real-world count
#' matrices. Preferring AnnotationDbi-confirmed IDs therefore maximises
#' compatibility with count matrices from sequencing providers whose pipelines
#' are not frequently updated.
#'
#' Within each group of rows sharing a \code{gene_name}, the following priority
#' order is applied until a single row is selected:
#' \enumerate{
#'   \item \strong{Pre-filter:} If any row has \code{hgnc_symbol_2 == gene_name}
#'     (AnnotationDbi confirms the symbol), rows with \code{hgnc_symbol_2 == NA}
#'     are discarded first. This ensures that an AnnotationDbi-confirmed row is
#'     never passed over in favour of an unconfirmed one merely because the
#'     latter happens to have \code{hgnc_symbol == gene_name} from BioMart.
#'   \item \strong{BioMart symbol match:} Rows where \code{hgnc_symbol == gene_name}
#'     (and is not a raw ENSG placeholder).
#'   \item \strong{AnnotationDbi symbol match:} Rows where
#'     \code{hgnc_symbol_2 == gene_name} (and is not a raw ENSG placeholder).
#'   \item \strong{Both sources agree:} Rows where
#'     \code{hgnc_symbol == hgnc_symbol_2}, indicating cross-source confirmation.
#'   \item \strong{BioMart ENSG confirmation:} Rows whose \code{ensembl_gene_id}
#'     matches the first entry in the \code{ensg_2} \code{///}-separated list
#'     returned by AnnotationDbi. Note that \code{ensg_2} list ordering is not
#'     considered a reliable preference signal on its own; this filter is
#'     intentionally placed after source-agreement filters.
#'   \item \strong{Drop ENSG placeholders:} Rows where \code{hgnc_symbol} is
#'     still a raw ENSG ID are deprioritised.
#'   \item \strong{Last resort:} When all disambiguation fields
#'     (\code{hgnc_symbol_2}, \code{ensg_2}) are \code{NA} across the entire
#'     group, the first row is taken. When rows are otherwise identical in all
#'     metadata, the newer ENSG ID (as returned by BioMart) is preferred as the
#'     more current annotation.
#' }
#'
#' The second pass (by \code{hgnc_symbol}) applies the same principle but
#' additionally prefers rows whose \code{hgnc_symbol} matches \code{gene_name},
#' and uses AnnotationDbi ENSG confirmation as a tiebreaker before falling back
#' to \code{x[1, ]}.
#'
#' \strong{ENSG placeholder resolution}
#'
#' After the filter chain, any remaining rows where \code{hgnc_symbol} is a raw
#' ENSG placeholder are fixed: if \code{hgnc_symbol_2} is available it is used;
#' otherwise \code{gene_name} is used (or \code{ensembl_gene_id} in ENSG-only
#' mode). This allows rows with ENSG placeholders from BioMart to be correctly
#' resolved in the second pass via their \code{hgnc_symbol_2} value.
#'
#' \strong{BioMart fallback}
#'
#' BioMart queries are attempted with graceful fallback through mirror hosts.
#' If all hosts fail the function proceeds with AnnotationDbi results only.
#' If both BioMart and AnnotationDbi fail entirely, the input is returned with
#' ENSG IDs used as \code{hgnc_symbol} values.
#'
#' @examples
#' \dontrun{
#' # Example input: two-column data frame with Ensembl IDs and gene symbols,
#' # as typically produced by a sequencing provider's count matrix annotation
#' my_genes <- data.frame(
#'   gene_id   = c("ENSG00000000003", "ENSG00000000419", "ENSG00000000460",
#'                 "ENSG00000012048", "ENSG00000075624", "ENSG00000111640",
#'                 "ENSG00000141510", "ENSG00000146648"),
#'   gene_name = c("TSPAN6", "DPM1", "FIRRM",
#'                 "BRCA1",  "ACTB",  "GAPDH",
#'                 "TP53",   "EGFR"),
#'   stringsAsFactors = FALSE
#' )
#'
#' # With gene symbols (full mode)
#' result <- unify_gene_ids(my_genes,
#'                          ensg_col   = "gene_id",
#'                          symbol_col = "gene_name",
#'                          verbose    = TRUE)
#'
#' # ENSG-only (e.g. from count matrix row names, no symbol column available)
#' ensg_only <- data.frame(
#'   ensembl_gene_id  = my_genes$gene_id,
#'   stringsAsFactors = FALSE
#' )
#' result_ensg <- unify_gene_ids(ensg_only, verbose = TRUE)
#' }
#'
#' @export
# =============================================================================
unify_gene_ids <- function(genes,
                           ensg_col         = "ensembl_gene_id",
                           symbol_col       = NULL,
                           host             = "https://www.ensembl.org",
                           biomart_fallback = c("https://uswest.ensembl.org",
                                                "https://asia.ensembl.org",
                                                "https://useast.ensembl.org"),
                           keep_intermediates = FALSE,
                           verbose          = FALSE) {

  # ---------------------------------------------------------------------------
  # Input validation
  # ---------------------------------------------------------------------------
  if (!is.data.frame(genes)) {
    if (!is.character(genes))
      stop("'genes' must be a data frame or character vector.")
    else
      genes <- data.frame(
        ensembl_gene_id  = genes,
        stringsAsFactors = FALSE
      )
  }
  if (!ensg_col %in% names(genes))
    stop(sprintf("Column '%s' not found in 'genes'.", ensg_col))
  if (!is.null(symbol_col) && !symbol_col %in% names(genes))
    stop(sprintf("Column '%s' not found in 'genes'.", symbol_col))
  if (nrow(genes) == 0L)
    stop("'genes' is empty.")

  # ---------------------------------------------------------------------------
  # Check for required Bioconductor packages
  # ---------------------------------------------------------------------------
  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE))
    stop(
      "The Bioconductor package 'org.Hs.eg.db' is required for this function.\n",
      "Install it via: BiocManager::install('org.Hs.eg.db')",
      call. = FALSE
    )
  if (!requireNamespace("AnnotationDbi", quietly = TRUE))
    stop(
      "The Bioconductor package 'AnnotationDbi' is required for this function.\n",
      "Install it via: BiocManager::install('AnnotationDbi')",
      call. = FALSE
    )

  has_symbols <- !is.null(symbol_col)
  n_input     <- nrow(genes)

  if (verbose) message(sprintf(
    "Input: %d rows, %s.", n_input,
    if (has_symbols) sprintf("symbol column '%s' provided", symbol_col)
    else "no symbol column (ENSG-only mode)"))

  # Rename columns internally for consistent handling
  names(genes)[names(genes) == ensg_col] <- "ensembl_gene_id"
  if (has_symbols) names(genes)[names(genes) == symbol_col] <- "gene_name"

  # ---------------------------------------------------------------------------
  # Lookups
  # ---------------------------------------------------------------------------
  if (verbose) message("Querying BioMart...")
  genes <- .try_biomart(genes, host, biomart_fallback, verbose)
  biomart_ok <- !all(is.na(genes$hgnc_symbol))
  if (verbose) {
    if (biomart_ok)
      message(sprintf("BioMart: resolved %d/%d symbols.",
        sum(!is.na(genes$hgnc_symbol)), n_input))
    else
      message("BioMart: no symbols resolved, proceeding with AnnotationDbi only.")
  }

  if (verbose) message("Querying AnnotationDbi (ENSG -> symbol)...")
  genes$hgnc_symbol_2 <- tryCatch(
    convertId2(genes$ensembl_gene_id),
    error = function(e) {
      warning(sprintf(
        "AnnotationDbi ENSG lookup failed ('%s'). Proceeding without hgnc_symbol_2.",
        conditionMessage(e)))
      NA_character_
    }
  )
  if (verbose)
    message(sprintf("AnnotationDbi (ENSG -> symbol): resolved %d/%d symbols.",
      sum(!is.na(genes$hgnc_symbol_2)), n_input))

  if (has_symbols) {
    if (verbose) message("Querying AnnotationDbi (symbol -> ENSG)...")
    genes$ensg_2 <- tryCatch(
      convertId2(genes$gene_name),
      error = function(e) {
        warning(sprintf(
          "AnnotationDbi symbol lookup failed ('%s'). Proceeding without ensg_2.",
          conditionMessage(e)))
        NA_character_
      }
    )
    if (verbose)
      message(sprintf("AnnotationDbi (symbol -> ENSG): resolved %d/%d entries.",
        sum(!is.na(genes$ensg_2)), n_input))
  } else {
    genes$ensg_2 <- NA_character_
  }

  # Early exit if all lookups failed
  if (!biomart_ok && all(is.na(genes$hgnc_symbol_2))) {
    warning("Both BioMart and AnnotationDbi lookups failed. ",
            "Returning input with ENSG IDs as hgnc_symbol.")
    # Use gene_name if it's a proper symbol, otherwise fall back to ensembl_gene_id
    needs_fix <- is.na(genes$hgnc_symbol) | grepl("^ENSG", genes$hgnc_symbol)
    if (has_symbols) {
      genes$hgnc_symbol[needs_fix] <- ifelse(
        !grepl("^ENSG", genes$gene_name[needs_fix]),
        genes$gene_name[needs_fix],
        genes$ensembl_gene_id[needs_fix]
      )
    } else {
      genes$hgnc_symbol[needs_fix] <- genes$ensembl_gene_id[needs_fix]
    }
    if (!keep_intermediates) {
      drop_cols <- intersect(c("hgnc_symbol_2", "ensg_2"), names(genes))
      genes     <- genes[, !names(genes) %in% drop_cols, drop = FALSE]
    }
    names(genes)[names(genes) == "ensembl_gene_id"] <- ensg_col
    if (has_symbols) names(genes)[names(genes) == "gene_name"] <- symbol_col
    return(genes)
  }

  has_ensg2 <- has_symbols && !all(is.na(genes$ensg_2))

  # ---------------------------------------------------------------------------
  # Deduplicate via internal function
  # ---------------------------------------------------------------------------
  genes2 <- .dedup_gene_ids(genes,
                            has_symbols = has_symbols,
                            has_ensg2   = has_ensg2,
                            verbose     = verbose)

  # Drop intermediate lookup columns unless explicitly requested
  if (!keep_intermediates) {
    drop_cols <- intersect(c("hgnc_symbol_2", "ensg_2"), names(genes2))
    genes2    <- genes2[, !names(genes2) %in% drop_cols, drop = FALSE]
  }

  # Restore original column names
  names(genes2)[names(genes2) == "ensembl_gene_id"] <- ensg_col
  if (has_symbols) names(genes2)[names(genes2) == "gene_name"] <- symbol_col

  genes2
}
