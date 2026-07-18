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
           verbose = FALSE)
  {
    if (use.cache) {
      cache <- .setCacheLocation(cache.dir = biom.cache)
      if (verbose) message("  Using biomaRt cache directory ", sQuote(cache))
    }
    biom <- match.arg(biom.mart)
    if (biom=="plants" && host == "https://www.ensembl.org") {
      if (verbose) message(sQuote("Plants"), "mart requested. Setting host to ", sQuote("https://plants.ensembl.org"), "...")
      host <- "https://plants.ensembl.org"
    }
    if (verbose) message("Getting CURL SSL options for securely contacting host ", sQuote(host), "...")
    httr_config <- .get.httr_config(host = host, use.cache = use.cache)
    marts <- biomaRt::listMarts(host=host, http_config=httr_config)[["biomart"]]
    marts1 <- sub("mart", "", tolower(marts))
    marts1 <- unlist(lapply(strsplit(tolower(marts1), "_"), function(x) x[length(x)]))
    biom <- marts[grep(biom, marts1)]
    if (verbose) message("Using BioMart: ", sQuote(biom))
    if (any(biom.data.set %in% c("human", "mouse"))) {
      biom.data.set <- match.arg(biom.data.set)
    }
    if (biom.data.set=="human") {
      if (biom=="ENSEMBL_MART_ENSEMBL") {
        if (verbose) message("Setting data set to ", sQuote("hsapiens_gene_ensembl"), "...")
        biom.data.set <- "hsapiens_gene_ensembl"
      } else {
        stop("'biom.mart' needs to be 'ensembl' to use data set 'human'!")
      }
    }
    if (biom.data.set=="mouse") {
      if (biom=="ENSEMBL_MART_ENSEMBL") {
        if (verbose) message("Setting data set to ", sQuote("mmusculus_gene_ensembl"), "...")
        biom.data.set <- "mmusculus_gene_ensembl"
      } else {
        stop("'biom.mart' needs to be 'ensembl' to use data set 'mouse'!")
      }
    }

    if (verbose) message("Input ID type is ", sQuote(biom.filter))
    mart <- biomaRt::useDataset(dataset=biom.data.set, mart=biomaRt::useMart(biomart=biom, host=host))

    if (!is.list(values)) {
      values <- as.character(values)
    }

    if (verbose) message("  Information requested: ", sQuote(setdiff(biom.attributes, biom.filter)), "...")
    biomaRt::getBM(attributes=biom.attributes, filters=biom.filter, values=values, mart=mart, useCache = use.cache)
  }
