# Retrieve Symbol Aliases and Previous symbols to determine a likely current symbol
#
# Author: vidal
###############################################################################

# Package-level environment used to cache the HGNC table across calls within
# the same R session. Avoids repeated downloads when likely_symbol() is called
# multiple times. The cache is reset automatically when the R session restarts.
.likely_symbol_cache <- new.env(hash = FALSE, parent = emptyenv())


#' Retrieve Symbol Aliases and Previous symbols to determine a likely current symbol
#' @description \command{likely_symbol()} downloads the latest version of the HGNC gene symbol database as a text
#'     file and query it to obtain symbol aliases, previous symbols and all symbols currently in use. (Optionally)
#'     assuming the input ID to be either an Alias or a Symbol or a Previous Symbol it performs multiple queries and
#'     compares the results of all possible combinations to determine a likely current Symbol.
#'     The downloaded HGNC table is cached for the duration of the R session to avoid repeated downloads.
#' @param syms (\code{character}). Vector of Gene Symbols to be tested.
#' @param alias_sym (\code{logical}). Should the input be assumed to be an Alias? Defaults to \code{TRUE}.
#' @param prev_sym (\code{logical}). Should the input be assumed to be a Previous Symbol? Defaults to \code{TRUE}.
#' @param orgnsm (\code{character}). The organism for which the Symbols are tested.
#' @param hgnc (\code{data.frame}). An optional data frame with the needed HGNC annotations. (Needs to match the
#'     format available at \code{hgnc_url}!) When supplied, bypasses both the cache and any download.
#' @param hgnc_url (\code{character}). URL where to download the HGNC annotation dataset. Defaults to
#'     \code{"https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt"}.
#' @param output (\code{character}). One of \code{"likely"}, \code{"symbols"} and \code{"all"}. Determines the scope
#'     of the output data frame. Defaults to \code{"likely"} which will return the input Symbol and the determined
#'     likely Symbol.
#' @param index_threshold (\code{integer}). Minimum number of unique input symbols above which inverted indices are
#'     pre-built for alias and previous symbol lookups, giving a substantial speedup for large inputs. Below this
#'     threshold the original row-scan is used, which is faster for very small inputs (e.g. a single symbol lookup)
#'     where the index-building overhead would dominate. Defaults to \code{10L}.
#' @param refresh (\code{logical}). Should the cached HGNC table be discarded and re-downloaded? Defaults to
#'     \code{FALSE}. Use \code{TRUE} to force a fresh download within the same R session, e.g. after a known
#'     HGNC update.
#' @param verbose (\code{logical}). Should messages be written to the console? Defaults to \code{TRUE}.
#' @details
#' The HGNC table is downloaded once per R session and cached in a package-level environment. Subsequent calls
#' reuse the cached table without any network access. If the cached table is more than 3 days old a warning message
#' is emitted recommending a refresh, since the HGNC database is updated monthly. To force a fresh download within
#' the same session use \code{refresh = TRUE} or start a new R session.
#'
#' When the number of unique input symbols is at or above \code{index_threshold}, inverted indices (hash tables)
#' are pre-built from the HGNC table so that each per-symbol lookup is O(1) rather than O(nrow(hgnc)), giving
#' roughly a 50-100x speedup for batch inputs. For small inputs the original row-scan is retained to avoid the
#' index-building overhead.
#' @return A \code{data.frame} with the following columns depending on the \code{output} setting.
#' \strong{\code{output="likely"}}:
#' \tabular{ll}{
#' \tab 'likely_symbol'\cr
#' \tab 'input_symbol'\cr
#' }
#' \strong{\code{output="symbols"}}:
#' \tabular{ll}{
#' \tab 'current_symbols'\cr
#' \tab 'likely_symbol'\cr
#' \tab 'input_symbol'\cr
#' \tab 'all_symbols'\cr
#' }
#' \strong{\code{output="all"}}:
#' \tabular{ll}{
#' \tab 'orig_input'\cr
#' \tab 'organism'\cr
#' \tab 'current_symbols'\cr
#' \tab 'likely_symbol'\cr
#' \tab 'input_symbol'\cr
#' \tab 'all_symbols'\cr
#' }
#' @note Only fully implemented for Human for now.
#' @examples
#' \dontrun{
#' # Single symbol lookup (uses row-scan, no index overhead)
#' likely_symbol("CCBL1")
#'
#' # Second call reuses cached HGNC table — no download
#' likely_symbol("KAAT1")
#'
#' # Force a fresh download within the same session
#' likely_symbol("CCBL1", refresh = TRUE)
#'
#' # Batch lookup (builds index for speed)
#' likely_symbol(c("ABCC4", "ACPP", "KIAA1524"))
#'
#' # Supply a pre-loaded table to bypass cache and download entirely
#' likely_symbol(c("ABCC4", "ACPP"), hgnc = my_hgnc_table)
#' }
#' @export
likely_symbol <-
  function(syms, alias_sym = TRUE, prev_sym = TRUE, orgnsm = "human", hgnc = NULL,
           hgnc_url = NULL, output = c("likely", "symbols", "all"),
           index_threshold = 10L, refresh = FALSE, verbose = TRUE)
  {
    if (is.null(orgnsm) || is.na(orgnsm))
      orgnsm <- ""

    if (!is.data.frame(syms)) {
      dq <- data.frame(Material = seq_along(syms), org = orgnsm,
                       sp = syms, stringsAsFactors = FALSE)
    } else {
      if (!all(grepl("Material|org|sp", names(syms))))
        stop("Need the following columns: 'Material', 'org', 'sp'")
      dq <- syms[syms$org %in% orgnsm, c("Material", "org", "sp")]
    }

    if (orgnsm == "human") {

      # ---------------------------------------------------------------------
      # HGNC table: use supplied table, cached table, or download fresh
      # ---------------------------------------------------------------------
      if (!is.null(hgnc)) {
        # User supplied table explicitly — use as-is, bypass cache entirely
        if (verbose) message("   ::: Using supplied HGNC table.", domain = NA)

      } else if (!refresh &&
                 exists("hgnc", envir = .likely_symbol_cache, inherits = FALSE)) {
        # Reuse cached table; warn if stale
        cached_age_days <- as.numeric(
          difftime(Sys.time(),
                   .likely_symbol_cache$downloaded_at,
                   units = "days")
        )
        if (cached_age_days > 3) {
          if (verbose)
            message(sprintf(
              paste0("   ::: Reusing cached HGNC table (downloaded %.0f days ago).\n",
                     "       Table may be outdated; use refresh = TRUE to re-download."),
              cached_age_days
            ), domain = NA)
        } else {
          if (verbose)
            message(sprintf(
              "   ::: Reusing cached HGNC table (downloaded %.1f hours ago).",
              cached_age_days * 24
            ), domain = NA)
        }
        hgnc <- .likely_symbol_cache$hgnc

      } else {
        # Download fresh (either first call or refresh = TRUE)
        if (refresh && verbose)
          message("   ::: Refreshing HGNC table...", domain = NA)
        if (is.null(hgnc_url))
          hgnc_url <- "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt"
        if (verbose)
          message(paste("   :::Downloading HGNC definition file from ",
                        sQuote(stringr::str_trunc(hgnc_url, 28, "right")), "> "),
                  domain = NA, appendLF = FALSE)
        hgnc <- read.delim(hgnc_url)
        if (verbose) message("done", domain = NA)
        # Store in session cache
        .likely_symbol_cache$hgnc          <- hgnc
        .likely_symbol_cache$downloaded_at <- Sys.time()
      }

      hg <- data.frame(
        symbol       = as.character(hgnc$symbol),
        alias_symbol = as.character(hgnc$alias_symbol),
        prev_symbol  = as.character(hgnc$prev_symbol),
        stringsAsFactors = FALSE
      )

      # ---------------------------------------------------------------------
      # Decide whether to build inverted indices based on input size
      # ---------------------------------------------------------------------
      n_unique  <- length(unique(dq$sp))
      use_index <- n_unique >= index_threshold

      if (use_index) {
        if (alias_sym) {
          if (verbose) message("   > Building alias symbol index...", domain = NA)
          alias_index <- new.env(hash = TRUE, parent = emptyenv())
          for (i in seq_len(nrow(hg))) {
            tokens <- trimws(strsplit(hg$alias_symbol[i], "\\|")[[1]])
            tokens <- tokens[nchar(tokens) > 0L]
            for (tok in tokens) {
              if (exists(tok, envir = alias_index, inherits = FALSE))
                alias_index[[tok]] <- c(alias_index[[tok]], i)
              else
                alias_index[[tok]] <- i
            }
          }
        }
        if (prev_sym) {
          if (verbose) message("   > Building previous symbol index...", domain = NA)
          prev_index <- new.env(hash = TRUE, parent = emptyenv())
          for (i in seq_len(nrow(hg))) {
            tokens <- trimws(strsplit(hg$prev_symbol[i], "\\|")[[1]])
            tokens <- tokens[nchar(tokens) > 0L]
            for (tok in tokens) {
              if (exists(tok, envir = prev_index, inherits = FALSE))
                prev_index[[tok]] <- c(prev_index[[tok]], i)
              else
                prev_index[[tok]] <- i
            }
          }
        }
      }
    }

    # ---------------------------------------------------------------------
    # Iterate through all supplied organisms and (unique) symbols
    # ---------------------------------------------------------------------

    dhg <- plyr::ddply(dq, "org", function(x) {
      if (unique(x$org) == "human") {
        if (verbose) message("   > Getting symbol aliases for 'Human'...", domain = NA)

        # ---------------------------------------------------------------------
        # Since we're testing human symbols they need to upper case
        # ---------------------------------------------------------------------

        x2        <- droplevels(x[!duplicated(toupper(x$sp)), ])
        unique.sp <- unique(toupper(x2$sp))

        hgd <- plyr::ldply(seq_along(unique.sp), function(y) {
          as  <- unique.sp[y]
          as1 <- strsplit(as, "\\|")[[1]]
          as2 <- paste(paste0("^", as1, "$"), collapse = "|")

          # ------------------------------------------------------------------
          # Alias lookup
          # ------------------------------------------------------------------
          if (alias_sym) {
            if (verbose)
              message(paste0("   >>> Querying symbol aliases for ", as,
                             " [", y, "/", length(unique.sp), "]..."),
                      domain = NA)
            if (use_index) {
              alias_rows <- unique(unlist(lapply(as1, function(s) {
                if (exists(s, envir = alias_index, inherits = FALSE))
                  alias_index[[s]]
                else
                  NULL
              })))
              if (length(alias_rows)) {
                hga1 <- do.call("rbind", lapply(alias_rows, function(z) {
                  data.frame(
                    symbol       = hg$symbol[z],
                    alias_symbol = trimws(strsplit(hg$alias_symbol[z], "\\|")[[1]])[1L],
                    stringsAsFactors = FALSE
                  )
                }))
              } else {
                hga1 <- NULL
              }
            } else {
              hga1 <- plyr::llply(seq_len(nrow(hg)), function(z) {
                hgz <- grep(as2, strsplit(hg$alias_symbol[z], "\\|")[[1]])
                if (length(hgz)) {
                  hgaz <- hg[z, ]
                  data.frame(symbol       = hgaz[[1]],
                             alias_symbol = strsplit(hgaz[[2]], "\\|")[[1]],
                             stringsAsFactors = FALSE)
                } else {
                  NULL
                }
              })
              if (length(unlist(hga1)))
                hga1 <- do.call("rbind", hga1)
              else
                hga1 <- NULL
            }
          } else {
            hga1 <- NULL
          }

          # ------------------------------------------------------------------
          # Current symbol lookup: vectorised grep (always fast)
          # ------------------------------------------------------------------
          hga2_idx <- grep(as2, hg$symbol)
          if (length(hga2_idx)) {
            hga2s <- ifelse(hg$symbol[hga2_idx]       == "", NA_character_, hg$symbol[hga2_idx])
            hga2a <- ifelse(hg$alias_symbol[hga2_idx] == "", NA_character_, hg$alias_symbol[hga2_idx])
            hga2p <- ifelse(hg$prev_symbol[hga2_idx]  == "", NA_character_, hg$prev_symbol[hga2_idx])
            hga2  <- data.frame(symbol       = hga2s,
                                alias_symbol = hga2a,
                                prev_symbol  = hga2p,
                                stringsAsFactors = FALSE)
          } else {
            hga2 <- hga2s <- hga2a <- hga2p <- NULL
          }

          # ------------------------------------------------------------------
          # Previous symbol lookup
          # ------------------------------------------------------------------
          if (prev_sym) {
            if (verbose)
              message(paste("   >>> Querying previous symbols for", as, "..."),
                      domain = NA)
            if (use_index) {
              prev_rows <- unique(unlist(lapply(as1, function(s) {
                if (exists(s, envir = prev_index, inherits = FALSE))
                  prev_index[[s]]
                else
                  NULL
              })))
              if (length(prev_rows)) {
                hga3 <- do.call("rbind", lapply(prev_rows, function(i) {
                  prev_toks <- trimws(strsplit(hg$prev_symbol[i], "\\|")[[1]])
                  prev_toks <- if (length(prev_toks)) prev_toks else ""
                  data.frame(symbol      = hg$symbol[i],
                             prev_symbol = prev_toks,
                             stringsAsFactors = FALSE)
                }))
                previous_sym <- hga3$prev_symbol
              } else {
                hga3         <- NULL
                previous_sym <- ""
              }
            } else {
              hga3 <- plyr::llply(seq_len(nrow(hg)), function(i) {
                hgi <- grep(as2, strsplit(hg$prev_symbol[i], "\\|")[[1]])
                if (length(hgi)) {
                  hgai      <- hg[i, ]
                  hgai_prev <- strsplit(hgai[[3]], "\\|")[[1]]
                  if (!length(hgai_prev)) hgai_prev <- ""
                  data.frame(symbol      = hgai[[1]],
                             prev_symbol = hgai_prev,
                             stringsAsFactors = FALSE)
                } else {
                  NULL
                }
              })
              if (length(unlist(hga3))) {
                hga3         <- do.call("rbind", hga3)
                previous_sym <- hga3$prev_symbol
              } else {
                hga3         <- NULL
                previous_sym <- ""
              }
            }
          } else {
            hga3         <- NULL
            previous_sym <- ""
          }

          # ------------------------------------------------------------------
          # Determine likely symbol (logic unchanged from original)
          # ------------------------------------------------------------------
          if (is.null(hga2)) {
            hga2_123 <- NULL
          } else {
            hga2_123 <- unlist(strsplit(unlist(hga2), "\\|"))
          }

          hga123 <- unique(c(as1, unlist(hga1), hga2_123, unlist(hga3)))
          hga123 <- hga123[nchar(hga123) > 0L]
          hga123 <- unlist(strsplit(hga123, "\\|"))

          current_sym <- paste(
            unique(c(hga1$symbol, hga2s, hga3$symbol, as1)),
            collapse = "|"
          )

          hga1_psbl <- hga2_psbl <- hga3_psbl <- NULL
          if (!is.null(hga1)) {
            hga1_match <- unlist(apply(hga1, 2, function(x) grep(as2, x)))
            if (length(hga1_match))
              hga1_psbl <- unique(hga1[hga1_match, "symbol"])
          }
          if (!is.null(hga2)) {
            hga2_match <- unlist(plyr::llply(
              strsplit(unlist(hga2), "\\|"), function(x) grep(as2, x)))
            if (length(hga2_match))
              hga2_psbl <- unique(hga2$symbol)
          }
          if (!is.null(hga3)) {
            hga3_match <- unlist(apply(hga3, 2, function(x) grep(as2, x)))
            if (length(hga3_match))
              hga3_psbl <- unique(hga3[hga3_match, "symbol"])
          }

          psbl_syms <- unique(c(hga1_psbl, hga2_psbl, hga3_psbl))
          psbl_alt  <- unique(c(hga2_psbl, hga3_psbl))

          if (length(psbl_syms) == 1L) {
            likely_sym <- psbl_syms
          } else if (length(psbl_alt) == 1L) {
            likely_sym <- psbl_alt
          } else {
            likely_sym <- as1
          }

          data.frame(
            current_sym  = current_sym,
            likely_sym   = likely_sym,
            previous_sym = paste(unique(previous_sym[nchar(previous_sym) > 0L]),
                                 collapse = "|"),
            all_symbols  = paste(sort(hga123), collapse = "|"),
            stringsAsFactors = FALSE
          )
        })

        x2$all_symbols     <- hgd$all_symbols
        x2$current_symbols <- hgd$current_sym
        x2$likely_symbol   <- hgd$likely_sym
        x2$previous_symbol <- hgd$previous_sym
        x <- merge(x, x2, by = "Material", all.x = TRUE)
        x <- plyr::ddply(x, "sp.x", function(spx) {
          spx$all_symbols[is.na(spx$all_symbols)] <-
            spx$all_symbols[!is.na(spx$all_symbols)][1]
          if (length(unique(spx$Material)) == 1L) spx[1, ] else spx
        })
        if (orgnsm == "human") x else x[, c("Material", "org.x", "sp.x", "all_symbols")]

      } else if (unique(x$org) == "mouse") {
        if (verbose) message("  Getting symbol aliases for 'Mouse'...", domain = NA)
        x2        <- droplevels(x[!duplicated(x$sp), ])
        unique.sp <- unique(x2$sp)
        mx2       <- convert.alias(unique.sp, species = "Mouse")
        mx2       <- paste(sort(mx2), collapse = "|")
        x2$all_symbols <- mx2
        x <- merge(x, x2, by = "Material", all.x = TRUE)
        plyr::ddply(x, "sp.x", function(spx) {
          spx$all_symbols[is.na(spx$all_symbols)] <-
            spx$all_symbols[!is.na(spx$all_symbols)][1]
          if (length(unique(spx$Material)) == 1L) spx[1, ] else spx
        })

      } else {
        if (verbose)
          message(paste("  Species", unique(x$org),
                        "not found. Returning unchanged input."), domain = NA)
        x$all_symbols <- x$sp
        names(x)      <- c("Material", "org.x", "sp.x", "all_symbols")
        plyr::ddply(x, "sp.x", function(spx) {
          spx$all_symbols[is.na(spx$all_symbols)] <-
            spx$all_symbols[!is.na(spx$all_symbols)][1]
          if (length(unique(spx$Material)) == 1L) spx[1, ] else spx
        })
      }
    })

    output <- match.arg(output)
    if (orgnsm == "human") {
      dhg <- dhg[, c("Material", "org.x", "current_symbols", "likely_symbol",
                     "previous_symbol", "sp.x", "all_symbols")]
      names(dhg) <- c("orig_input", "organism", "current_symbols", "likely_symbol",
                      "previous_symbol", "input_symbol", "all_symbols")
      cat("\n\n")
      if (output == "likely") {
        return(dhg[, c("likely_symbol", "input_symbol")])
      } else if (output == "symbols") {
        return(dhg[, c("current_symbols", "likely_symbol", "input_symbol", "all_symbols")])
      } else {
        return(dhg[, c("current_symbols", "likely_symbol", "previous_symbol",
                       "input_symbol", "all_symbols")])
      }
    } else {
      dhg        <- dhg[, c("Material", "org.x", "sp.x", "all_symbols")]
      names(dhg) <- c("orig_input", "organism", "input_symbol", "all_symbols")
      if (verbose) cat("\n\n")
      return(dhg)
    }
  }
