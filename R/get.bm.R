#' Make a Query to Biomart.
#' @description \command{get.bm()} is a user-friendly wrapper for \command{getBM()} from the \emph{biomaRt} package with default
#'     settings for Human and Mouse.
#' It sets all needed variables and performs the query.
#' @param values \code{character} vector of ids to be converted.
#' @param biom.data.set \code{character} of length one. Biomart data set to use. Defaults to 'human' (internally translated to "hsapiens_gene_ensembl" if \code{biom.mart="ensembl"}).
#' @param biom.mart \code{character} vector. Biomart to use (uses the first element of the vector), defaults to "ensembl".
#' @param host \code{character} of length one. Host URL.
#' @param biom.filter \code{character} of length one. Name of biomart filter, i.e., type of query ids, defaults to "ensembl_gene_id".
#' @param biom.attributes \code{character} vector. Biomart attributes, i.e., type of desired result(s); make sure query id type is included!
#' @param biom.cache \code{character}. Path name giving the location of the cache \command{getBM()} uses if \code{use.cache=TRUE}. Defaults to the value in the \emph{BIOMART_CACHE} environment variable.
#' @param use.cache (\code{logical}). Should \command{getBM()} use the cache? Defaults to \code{TRUE} as in the \command{getBM()} function and is passed on to that.
#' @param biomart.fallback \code{character} vector. Fallback host URLs to try if the primary
#'   \code{host} fails. Set to \code{NULL} to disable fallback. Defaults to Ensembl mirror sites.
#' @param chunk.size \code{integer} of length one. Maximum number of IDs per BioMart query.
#'   Large ID lists are split into chunks of this size to avoid server timeouts.
#'   Set to \code{Inf} to disable chunking. Defaults to \code{500}.
#' @param verbose (\code{logical}). Should verbose output be written to the console? Defaults to \code{FALSE}.
#' @return  A data frame with the retrieved information.
#' @author Vidal Fey
#' @seealso \command{\link[biomaRt]{getBM}}
#' @examples
#' \dontrun{
#' val <- c("ENSG00000111199", "ENSG00000134121", "ENSG00000176102", "ENSG00000171611")
#' bm <- get.bm(val)
#' bm
#' }
#' @keywords utilities
#' @export
get.bm <-
  function(values,
           biom.data.set = c("human", "mouse"),
           biom.mart = c("ensembl", "mouse", "snp", "funcgen", "plants"),
           host = "https://www.ensembl.org",
           biom.filter = "ensembl_gene_id",
           biom.attributes = c("ensembl_gene_id",
                               "hgnc_symbol", "description"),
           biom.cache = rappdirs::user_cache_dir("biomaRt"),
           use.cache = TRUE,
           biomart.fallback = c("https://useast.ensembl.org",
                                "https://uswest.ensembl.org",
                                "https://asia.ensembl.org"),
           chunk.size = 500L,
           verbose = FALSE)
  {
    # -------------------------------------------------------------------------
    # Input validation: detect versioned Ensembl IDs
    # -------------------------------------------------------------------------
    if (!is.list(values)) {
      values <- as.character(values)
    }

    if (identical(biom.filter, "ensembl_gene_id")) {
      vals_check <- if (is.list(values)) unlist(values) else values
      versioned <- grep("\\.[0-9]+$", vals_check, value = TRUE)
      if (length(versioned) > 0L) {
        n_ver <- length(versioned)
        examples <- paste(utils::head(versioned, 3L), collapse = ", ")
        if (n_ver > 3L) examples <- paste0(examples, ", ...")
        stop(
          sprintf(
            "%d of %d input IDs appear to be versioned Ensembl IDs (e.g. %s).\n",
            n_ver, length(vals_check), examples
          ),
          "BioMart's 'ensembl_gene_id' filter expects unversioned IDs ",
          "(e.g. 'ENSG00000111199', not 'ENSG00000111199.5').\n",
          "Strip version suffixes first, e.g.: sub('\\\\.[0-9]+$', '', ids)",
          call. = FALSE
        )
      }
    }

    # -------------------------------------------------------------------------
    # Cache setup
    # -------------------------------------------------------------------------
    if (use.cache) {
      cache <- .setCacheLocation(cache.dir = biom.cache)
      if (verbose) message("  Using biomaRt cache directory ", sQuote(cache))
    }

    # -------------------------------------------------------------------------
    # Resolve mart and dataset
    # -------------------------------------------------------------------------
    biom <- match.arg(biom.mart)
    if (biom == "plants" && host == "https://www.ensembl.org") {
      if (verbose) message(sQuote("Plants"), " mart requested. Setting host to ", sQuote("https://plants.ensembl.org"), "...")
      host <- "https://plants.ensembl.org"
      # Plants mirrors don't exist; disable fallback
      biomart.fallback <- NULL
    }

    # -------------------------------------------------------------------------
    # Resolve dataset name
    # -------------------------------------------------------------------------
    if (any(biom.data.set %in% c("human", "mouse"))) {
      biom.data.set <- match.arg(biom.data.set)
    } else {
      # A data set named by the user. Reduce to length one here: a multi-element
      # value would otherwise reach an if() condition below and abort on
      # R >= 4.2 with "the condition has length > 1".
      biom.data.set <- as.character(biom.data.set)[1L]
    }
    biom.data.set.in <- biom.data.set
    biom.data.set <- switch(biom.data.set,
                            human = "hsapiens_gene_ensembl",
                            mouse = "mmusculus_gene_ensembl",
                            biom.data.set)
    if (verbose && !identical(biom.data.set, biom.data.set.in))
      message("Setting data set to ", sQuote(biom.data.set), "...")

    if (verbose) message("Input ID type is ", sQuote(biom.filter))

    # -------------------------------------------------------------------------
    # Query with host fall-back and return results if any
    # -------------------------------------------------------------------------
    .with_biomart_fallback(
      fn = function(h) {
        mart <- .connect_mart(h, biom, biom.data.set, use.cache, verbose)
        if (verbose) message("  Information requested: ",
                             paste(sQuote(setdiff(biom.attributes, biom.filter)), collapse = ", "), "...")
        .chunked_getBM(mart, values, biom.attributes, biom.filter,
                       use.cache, chunk.size, verbose)
      },
      host = host,
      fallback_hosts = biomart.fallback,
      verbose = verbose
    )
  }
