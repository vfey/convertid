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
    # Helper: connect to a mart on a given host
    # -------------------------------------------------------------------------
    .connect_mart <- function(h, biom, biom.data.set, use.cache, verbose) {
      if (verbose) message("Getting CURL SSL options for securely contacting host ", sQuote(h), "...")
      httr_config <- .get.httr_config(host = h, use.cache = use.cache)
      marts <- biomaRt::listMarts(host = h, http_config = httr_config)[["biomart"]]
      marts1 <- sub("mart", "", tolower(marts))
      marts1 <- unlist(lapply(strsplit(tolower(marts1), "_"), function(x) x[length(x)]))
      biom_name <- marts[grep(biom, marts1)]
      if (verbose) message("Using BioMart: ", sQuote(biom_name))
      biomaRt::useDataset(
        dataset = biom.data.set,
        mart = biomaRt::useMart(biomart = biom_name, host = h)
      )
    }

    # -------------------------------------------------------------------------
    # Resolve dataset name
    # -------------------------------------------------------------------------
    if (any(biom.data.set %in% c("human", "mouse"))) {
      biom.data.set <- match.arg(biom.data.set)
    }
    if (biom.data.set == "human") {
      biom.data.set <- "hsapiens_gene_ensembl"
      if (verbose) message("Setting data set to ", sQuote(biom.data.set), "...")
    } else if (biom.data.set == "mouse") {
      biom.data.set <- "mmusculus_gene_ensembl"
      if (verbose) message("Setting data set to ", sQuote(biom.data.set), "...")
    }

    if (verbose) message("Input ID type is ", sQuote(biom.filter))

    # -------------------------------------------------------------------------
    # Helper: chunked getBM query
    # -------------------------------------------------------------------------
    .chunked_getBM <- function(mart, values, biom.attributes, biom.filter,
                               use.cache, chunk.size, verbose) {
      vals <- if (is.list(values)) values else as.character(values)
      # Lists (e.g. for multi-value filters) are not chunked
      if (is.list(vals) || length(vals) <= chunk.size) {
        if (verbose) message("  Querying ", length(vals), " IDs in a single batch...")
        return(
          biomaRt::getBM(
            attributes = biom.attributes,
            filters    = biom.filter,
            values     = vals,
            mart       = mart,
            useCache   = use.cache
          )
        )
      }
      chunks <- split(vals, ceiling(seq_along(vals) / chunk.size))
      if (verbose) message(
        "  Splitting ", length(vals), " IDs into ", length(chunks),
        " chunks of up to ", chunk.size, " IDs..."
      )
      results <- vector("list", length(chunks))
      for (i in seq_along(chunks)) {
        if (verbose) message(
          "  Chunk ", i, "/", length(chunks),
          " (", length(chunks[[i]]), " IDs)..."
        )
        results[[i]] <- biomaRt::getBM(
          attributes = biom.attributes,
          filters    = biom.filter,
          values     = chunks[[i]],
          mart       = mart,
          useCache   = use.cache
        )
      }
      do.call(rbind, results)
    }

    # -------------------------------------------------------------------------
    # Query with host fallback
    # -------------------------------------------------------------------------
    all_hosts <- c(host, biomart.fallback)
    result <- NULL

    for (h in all_hosts) {
      if (verbose) message("Trying host: ", sQuote(h), "...")
      result <- tryCatch(
        {
          mart <- .connect_mart(h, biom, biom.data.set, use.cache, verbose)
          if (verbose) message("  Information requested: ", sQuote(setdiff(biom.attributes, biom.filter)), "...")
          .chunked_getBM(mart, values, biom.attributes, biom.filter,
                         use.cache, chunk.size, verbose)
        },
        error = function(e) {
          if (verbose) message("  Host ", sQuote(h), " failed: ", conditionMessage(e))
          NULL
        }
      )
      if (!is.null(result)) {
        if (verbose && h != host)
          message("Succeeded with fallback host ", sQuote(h), ".")
        return(result)
      }
    }

    # All hosts failed
    hosts_tried <- paste(sQuote(all_hosts), collapse = ", ")
    stop(
      "All BioMart hosts failed (tried: ", hosts_tried, ").\n",
      "Check your internet connection or try again later.",
      call. = FALSE
    )
  }
