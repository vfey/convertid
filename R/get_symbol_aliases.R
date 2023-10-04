# TODO: Add comment
#
# Author: vidal
###############################################################################

#' Retrieve Symbol Aliases and Previous symbols to determine a likely current symbol
#' @description \command{likely_symbol()} downloads the latest version of the HGNC gene symbol database as a text
#'     file and query it to obtain symbol aliases, previous symbols and all symbols currently in use. (Optionally)
#'     assuming the input ID to be either an Alias or a Symbol or a Previous Symbol it performs multiple queries and
#'     compares the results of all possible combinations to determine a likely current Symbol.
#' @param syms (\code{character}). Vector of Gene Symbols to be tested.
#' @param alias_sym (\code{logical}). Should the input be assumed to be an Alias? Defaults to \code{TRUE}.
#' @param prev_sym (\code{logical}). Should the input be assumed to be a Previous Symbol? Defaults to \code{TRUE}.
#' @param orgnsm (\code{character}). The organism for which the Symbols are tested.
#' @param hgnc (\code{data.frame}). An optional data frame with the needed HGNC annotations. (Needs to match the format
#'     available at \code{hgnc_utl}!)
#' @param hgnc_url (\code{character}). URL where to download the HGNC annotation dataset. Defaults to
#'     \code{"ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/hgnc_complete_set.txt"}.
#' @param output (\code{character}). One of "likely", "symbols" and "all". Determines the scope of the output data frame.
#'     Defaults to \code{"likely"} which will return the inout Symbol and the determined likely Symbol.
#' @param verbose (\code{logical}). Should messages be written to the console? Defaults to \code{TRUE}.
#' @details Please note that the algorithm is very slow for large input vectors.
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
#' likely_symbol(c("ABCC4", "ACPP", "KIAA1524"))
#' }
#' @export
likely_symbol <-
  function(syms, alias_sym=TRUE, prev_sym=TRUE, orgnsm="human", hgnc=NULL,
           hgnc_url = NULL, output = c("likely", "symbols", "all"), verbose=TRUE)
  {
    if (is.null(orgnsm) || is.na(orgnsm))
      orgnsm <- ""
    if (!is.data.frame(syms)) {
      dq <- data.frame(Material=seq_along(syms), org=orgnsm, sp=toupper(syms), stringsAsFactors = FALSE)
    } else {
      if (!length(grep("Material&org&sp", names(syms))))
        stop("Need the following columns: 'Material', 'org', 'sp'")
      dq <- syms[syms$org %in% orgnsm, c("Material", "org", "sp")]
    }
    if (orgnsm=="human") {
      if (is.null(hgnc)) {
        if (is.null(hgnc_url)) {
          hgnc_url <- "ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/hgnc_complete_set.txt"
        }
        if (verbose) message(paste("   :::Downloading HGNC definition file from ", sQuote(stringr::str_trunc(hgnc_url, 28, "right")), ">>> "), domain=NA, appendLF = FALSE)
        hgnc <- read.delim(hgnc_url)
        if (verbose) message("done", domain=NA)
      }
      hg <- data.frame(symbol=as.character(hgnc$symbol), alias_symbol=as.character(hgnc$alias_symbol),
                       prev_symbol=as.character(hgnc$prev_symbol), stringsAsFactors=FALSE)
    }
    dhg <- plyr::ddply(dq, "org", function(x) {
      if (unique(x$org)=="human") {
        if (verbose) message(paste("  Getting symbol aliases for 'Human'..."), domain=NA)
        x2 <- droplevels(x[!duplicated(x$sp), ])
        unique.sp <- unique(x2$sp)
        hgd <- plyr::ldply(1:length(unique.sp), function(y) {
          as <- unique.sp[y]
          as1 <- strsplit(as, "\\|")[[1]]
          as2 <- paste(paste0("^", as1, "$"), collapse="|")
          if (alias_sym) {
            if (verbose) message(paste0("   >>> Querying symbol aliases for ", as, " [", y, "/", length(unique.sp), "]", "..."), domain=NA)
            hga1 <- plyr::llply(1:nrow(hg), function(z) {
              hgz <- grep(as2, strsplit(hg$alias_symbol[z], "\\|")[[1]])
              if (length(hgz)) {
                hgaz <- hg[z, ]
                data.frame(symbol=hgaz[[1]], alias_symbol=strsplit(hgaz[[2]], "\\|")[[1]], stringsAsFactors = FALSE)
              } else {
                NULL
              }
            })
            if (length(unlist(hga1))) {
              hga1 <- do.call("rbind", hga1)
            } else {
              hga1 <- NULL
            }
          } else {
            hga1 <- NULL
          }
          hga2 <- grep(as2, hg$symbol)
          if (length(hga2)) {
            hga2s <- ifelse(hg$symbol[hga2]=="", NA_character_, hg$symbol[hga2])
            hga2a <- ifelse(hg$alias_symbol[hga2]=="", NA_character_, hg$alias_symbol[hga2])
            hga2p <- ifelse(hg$prev_symbol[hga2]=="", NA_character_, hg$prev_symbol[hga2])
            hga2 <- data.frame(symbol=hga2s, alias_symbol=hga2a, prev_symbol=hga2p, stringsAsFactors = FALSE)
          } else {
            hga2 <- hga2s <- hga2a <- hga2p <- NULL
          }

          if (prev_sym) {
            if (verbose) message(paste("   >>> Querying previous symbols for", as, "..."), domain=NA)
            hga3 <- plyr::llply(1:nrow(hg), function(i) {
              hgi <- grep(as2, strsplit(hg$prev_symbol[i], "\\|")[[1]])
              if (length(hgi)) {
                hgai <- hg[i, ]
                hgai_prev <- strsplit(hgai[[3]], "\\|")[[1]]
                if (!length(hgai_prev))
                  hgai_prev <- ""
                data.frame(symbol=hgai[[1]], prev_symbol=hgai_prev, stringsAsFactors = FALSE)
              } else {
                NULL
              }
            })
            if (length(unlist(hga3))) {
              hga3 <- do.call("rbind", hga3)
            } else {
              hga3 <- NULL
            }
          } else {
            hga3 <- NULL
          }
          if (is.null(hga2)) {
            hga2_123 <- NULL
          } else {
            hga2_123 <- unlist(strsplit(unlist(hga2), "\\|"))
          }
          hga123 <- unique(c(as1, unlist(hga1), hga2_123, unlist(hga3)))
          hga123 <- hga123[hga123!=""]
          hga123 <- unlist(strsplit(hga123, "\\|"))
          current_sym <- paste(unique(c(hga1$symbol, hga2s, hga3$symbol, as1)), collapse="|")
          hga1_psbl <- hga2_psbl <- hga3_psbl <- NULL
          if (!is.null(hga1)) {
            hga1_match <- unlist(apply(hga1, 2, function(x) grep(as2, x)))
            if (length(hga1_match)) {
              hga1_psbl <- unique(c(hga1[hga1_match, "symbol"]))
            }
          }
          if (!is.null(hga2)) {
            hga2_match <- unlist(plyr::llply(strsplit(unlist(hga2), "\\|"), function(x) grep(as2, x)))
            if (length(hga2_match)) {
              hga2_psbl <- unique(c(hga2$"symbol"))
            }
          }
          if (!is.null(hga3)) {
            hga3_match <- unlist(apply(hga3, 2, function(x) grep(as2, x)))
            if (length(hga3_match)) {
              hga3_psbl <- unique(c(hga3[hga3_match, "symbol"]))
            }
          }
          psbl_syms <- unique(c(hga1_psbl, hga2_psbl, hga3_psbl))
          psbl_alt <- unique(c(hga2_psbl, hga3_psbl))
          if (length(psbl_syms)==1) {
            likely_sym <- psbl_syms
          } else if (length(psbl_alt)==1) {
            likely_sym <- psbl_alt
          } else {
            likely_sym <- as1
          }
          data.frame(current_sym=current_sym, likely_sym=likely_sym, all_symbols=paste(sort(hga123), collapse="|"), stringsAsFactors = FALSE)
        })
        x2$all_symbols <- hgd$all_symbols
        x2$current_symbols <- hgd$current_sym
        x2$likely_symbol <- hgd$likely_sym
        x <- merge(x, x2, by="Material", all.x=TRUE)
        x <- plyr::ddply(x, "sp.x", function(spx) {
          spx$all_symbols[is.na(spx$all_symbols)] <- spx$all_symbols[!is.na(spx$all_symbols)][1]
          if (length(unique(spx$Material))==1)
            spx[1, ]
          else
            spx
        })
        if (orgnsm=="human") {
          x
        } else {
          x[, c("Material", "org.x", "sp.x", "all_symbols")]
        }
      } else if (unique(x$org)=="mouse") {
        if (verbose) message(paste("  Getting symbol aliases for 'Mouse'..."), domain=NA)
        x2 <- droplevels(x[!duplicated(x$sp), ])
        unique.sp <- unique(x2$sp)
        mx2 <- convert.alias(unique.sp, species="Mouse")
        mx2 <- paste(sort(mx2), collapse="|")
        x2$all_symbols <- mx2
        x <- merge(x, x2, by="Material", all.x=TRUE)
        plyr::ddply(x, "sp.x", function(spx) {
          spx$all_symbols[is.na(spx$all_symbols)] <- spx$all_symbols[!is.na(spx$all_symbols)][1]
          if (length(unique(spx$Material))==1)
            spx[1, ]
          else
            spx
        })
      } else if (unique(x$org)=="rat") {
        if (verbose) message(paste("  Getting symbol aliases for 'Rat'... (not implemented)"), domain=NA)
        x$all_symbols <- x$sp
        names(x) <- c("Material", "org.x", "sp.x", "all_symbols")
        plyr::ddply(x, "sp.x", function(spx) {
          spx$all_symbols[is.na(spx$all_symbols)] <- spx$all_symbols[!is.na(spx$all_symbols)][1]
          if (length(unique(spx$Material))==1)
            spx[1, ]
          else
            spx
        })
      } else if (unique(x$org)=="pig") {
        if (verbose) message(paste("  Getting symbol aliases for 'Pig'... (not implemented)"), domain=NA)
        x$all_symbols <- x$sp
        names(x) <- c("Material", "org.x", "sp.x", "all_symbols")
        plyr::ddply(x, "sp.x", function(spx) {
          spx$all_symbols[is.na(spx$all_symbols)] <- spx$all_symbols[!is.na(spx$all_symbols)][1]
          if (length(unique(spx$Material))==1)
            spx[1, ]
          else
            spx
        })
      } else {
        if (verbose) message(paste("  Species", unique(x$org), "not found. Returning unchanged input."), domain=NA)
        x$all_symbols <- x$sp
        names(x) <- c("Material", "org.x", "sp.x", "all_symbols")
        plyr::ddply(x, "sp.x", function(spx) {
          spx$all_symbols[is.na(spx$all_symbols)] <- spx$all_symbols[!is.na(spx$all_symbols)][1]
          if (length(unique(spx$Material))==1)
            spx[1, ]
          else
            spx
        })
      }
    })
    output <- match.arg(output)
    if (orgnsm=="human") {
      dhg <- dhg[, c("Material", "org.x", "current_symbols", "likely_symbol", "sp.x", "all_symbols")]
      names(dhg) <- c("orig_input", "organism", "current_symbols", "likely_symbol", "input_symbol", "all_symbols")
      if (output == "likely") {
        return(dhg[, c("likely_symbol", "input_symbol")])
      } else if (output == "symbols") {
        return(dhg[, c("current_symbols", "likely_symbol", "input_symbol", "all_symbols")])
      } else {
      }
    } else {
      dhg <- dhg[, c("Material", "org.x", "sp.x", "all_symbols")]
      names(dhg) <- c("orig_input", "organism", "input_symbol", "all_symbols")
    }
    return(dhg)
  }
