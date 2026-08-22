#' Unexported functions
#' Set the location for the \emph{\code{biomaRt}} cache
#' @description \command{.setBiomaRtCacheLocation()} attempts to set the cache location
#' used by the functions in the \emph{\code{biomaRt}} package and defined in the BIOMART_CACHE
#' environment variable.
#' If that variable is set and the defined location exists and is writable nothing is done.
#' If the system default cache location exists and is writable a sub-folder \code{app} is used (and created if necessary).
#' If the above don't work a new path is constructed from \code{cache.dir} and the \code{app} folder and an attempt is made to create that.
#' If all of the above fail the function attempts to create \code{file.path(tempdir(), app)}. If tat fails, too,
#' an exception is thrown.
#' @param cache.dir (\code{character}). Optional user-defined path used as cache parent directory. Defaults to \code{rappdirs::user_cache_dir()}. If a custom path is given that has to exist.
#' @param app (\code{character}). Optional application-specific cache sub-directory, i.e., the actually used cache location. Defaults to "biomaRt".
#' @return The value of the BIOMART_CACHE environment variable, i.e., the cache location.
#' @seealso \code{\link[rappdirs]{user_cache_dir}}, \code{\link[BiocFileCache]{BiocFileCache}}
#' @examples
#' \dontrun{.setCacheLocation()}
#' @keywords internal
.setCacheLocation <-
  function(cache.dir = rappdirs::user_cache_dir(), app = "biomaRt") {
    assertthat::assert_that(is.character(cache.dir), length(cache.dir) == 1L, !is.na(cache.dir),
              is.character(app), length(app) == 1L, !is.na(app))
    if (app == basename(cache.dir)) {
      cache.dir <- dirname(cache.dir)
    }
    system.cache <- rappdirs::user_cache_dir()
    current.cache <- Sys.getenv(x = "BIOMART_CACHE")
    if (.cache.writable(current.cache)) {
      cache.path <- current.cache
    } else if (.cache.writable(system.cache)) {
      cache.path <- rappdirs::user_cache_dir(app)
    } else {
      cache.path <- file.path(cache.dir, app)
    }
    path.ok <- .create.cache(cache.path)
    if (!path.ok) {
      cache.path <- file.path(tempdir(), app)
      message("Using temporary cache ", sQuote(cache.path), domain = NA)
      path.ok <- .create.cache(cache.path)
    }
    if (path.ok) {
      Sys.setenv(BIOMART_CACHE = cache.path)
    } else {
      stop("Could not create cache directory. Make sure the location in 'cache.dir' is writable.")
    }
    return(path.expand(Sys.getenv(x = "BIOMART_CACHE")))
  }

#' Unexported functions
#' Test if a path exists and is writable
#' @description \command{.cache.writable()} uses \command{file.access()} to test if a
#' given location exists and is writable by the user.
#' @param path (\code{character}). The path to be tested.
#' @return TRUE if both conditions are met, FALSE if not.
#' @seealso \code{\link[base]{file.access}}
#' @examples
#' \dontrun{.cache.writable(rappdirs::user_cache_dir())}
#' @keywords internal
.cache.writable <-
  function(path) {
    assertthat::assert_that(is.character(path), length(path) == 1L, !is.na(path))
    as.logical(file.access(path)+1) &&
      as.logical(file.access(path, 2)+1)
  }

#' Unexported functions
#' Create a file cache directory at a given location.
#' @description \command{.create.cache()} attempts to create a cache directory based on a given path name. Typically, such path
#' is specific to the package from within the function is called. The default settings refer to the file cache framework in the \emph{\code{biomaRt}} package.
#' @param cache.path (\code{character}). The path to use for the cached files.
#' @return TRUE if the location was successfully set up, FALSE if not.
#' @seealso \code{\link[rappdirs]{user_cache_dir}}
#' @examples
#' \dontrun{.create.cache(rappdirs::user_cache_dir("biomaRt"))}
#' @keywords internal
.create.cache <-
  function(cache.path = rappdirs::user_cache_dir("biomaRt")) {
    assertthat::assert_that(is.character(cache.path), length(cache.path) == 1L, !is.na(cache.path))
    if (!file.exists(cache.path)) {
      return(dir.create(cache.path, recursive = TRUE, showWarnings = FALSE))
    } else {
      .cache.writable(cache.path)
    }
  }

#' Unexported functions
#' Get httr configuration, i.e., current CURL options for data fetching functions.
#' @description \command{.get.httr_config()} retrieves the current CURL options and in particular tests and gets
#' the options used with \emph{"^https://.*ensembl.org"} URLs. The code was partly copied from \code{\link[biomaRt]{listMarts}()}.
#' @param httr_config (\code{list}). A R object of class \code{request} listing current CURL options. Missing (and meant to be created) by default.
#' @param host (\code{character}) Host URL.
#' @param use.cache (\code{logical}) Should \emph{\code{biomaRt}} functions use the cache? Defaults to \code{TRUE}.
#' @return A R object of class \code{request} listing current CURL options.
#' @seealso \code{\link[httr]{config}}
#' @examples
#' \dontrun{.get.httr_config()}
#' @keywords internal
.get.httr_config <-
  function(httr_config, host="https://www.ensembl.org", use.cache = TRUE) {
    if (grepl(pattern = "^https://.*ensembl.org", x = host) && missing(httr_config)) {
      httr_config <- .get.Ensembl_config(use.cache = use.cache)
    }
    if (missing(httr_config)) {
      httr_config <- httr::config()
    }
    return(httr_config)
  }

#' Unexported functions
#' Test and retrieve Ensembl-specific CURL SSL configuration.
#' @description \command{.get.Ensembl_config()} tests and gets CURL options used with \emph{"^https://.*ensembl.org"} URLs.
#' The function is a modified version of \code{.getEnsemblSSL} from the \emph{\code{biomaRt}} package.
#' @param use.cache (\code{logical}) Should \emph{\code{biomaRt}} functions use the cache? Defaults to \code{TRUE}.
#' @return A R object of class \code{request} listing current CURL options.
#' @seealso \code{\link[httr]{config}}
#' @examples
#' \dontrun{.get.Ensembl_config()}
#' @keywords internal
.get.Ensembl_config <-
  function (use.cache = TRUE) {
    if (use.cache) {
      cache <- .biomartCacheLocation()
      bfc <- BiocFileCache::BiocFileCache(cache, ask = FALSE)
      if (.checkInCache(bfc, hash = "ensembl-ssl-settings")) {
        ensembl_config <- .readFromCache(bfc, "ensembl-ssl-settings")
      } else {
        ensembl_config <- .checkEnsemblSSL()
        .addToCache(bfc, ensembl_config, hash = "ensembl-ssl-settings")
      }
    } else {
      ensembl_config <- .checkEnsemblSSL()
    }
    return(ensembl_config)
  }

#' Try BioMart query with fallback hosts
#'
#' Iterates through \code{host} and \code{fallback_hosts}, returning the first
#' successful result from \code{\link{convert.bm}}. If all hosts fail, emits a
#' warning and returns \code{genes} with \code{hgnc_symbol} set to
#' \code{NA_character_}.
#'
#' @param genes Data frame with an \code{ensembl_gene_id} column.
#' @param host Primary BioMart host URL.
#' @param fallback_hosts Character vector of fallback host URLs.
#' @param verbose Logical; emit progress messages.
#' @return Data frame with BioMart annotations, or input with NA symbols on failure.
#' @keywords internal
try_biomart <- function(genes, host, fallback_hosts, verbose) {
  tryCatch(
    convert.bm(genes, id = "ensembl_gene_id", host = host,
               biomart.fallback = fallback_hosts, verbose = verbose),
    error = function(e) {
      warning("All BioMart hosts failed. Proceeding without BioMart results.")
      genes$hgnc_symbol <- NA_character_
      genes
    }
  )
}

#' Try a function across BioMart hosts with fallback
#'
#' @param fn Function that takes a host URL and returns a result.
#' @param host Primary host URL.
#' @param fallback_hosts Character vector of fallback URLs, or NULL.
#' @param verbose Logical.
#' @return The first successful result from \code{fn}.
#' @keywords internal
.with_biomart_fallback <- function(fn, host, fallback_hosts, verbose) {
  all_hosts <- c(host, fallback_hosts)
  last_error <- NULL
  for (h in all_hosts) {
    if (verbose) message("Trying host: ", sQuote(h), "...")
    result <- tryCatch(
      fn(h),
      error = function(e) {
        if (verbose) message("  Host ", sQuote(h), " failed: ", conditionMessage(e))
        last_error <<- e
        NULL
      }
    )
    if (!is.null(result)) {
      if (verbose && h != host)
        message("Succeeded with fallback host ", sQuote(h), ".")
      return(result)
    }
  }
  hosts_tried <- paste(sQuote(all_hosts), collapse = ", ")
  stop(
    "All BioMart hosts failed (tried: ", hosts_tried, ").\n",
    "Last error: ", conditionMessage(last_error),
    call. = FALSE
  )
}

#' Connect to a mart on a given host
#'
#' @param h (\code{character}) Host URL.
#' @param biom \code{character} vector. Biomart to use (uses the first element of the vector), defaults to "ensembl".
#' @param biom.data.set \code{character} of length one. Biomart data set to use. Defaults to 'human' (internally translated to "hsapiens_gene_ensembl" if \code{biom.mart="ensembl"}).
#' @param use.cache (\code{logical}) Should \emph{\code{biomaRt}} functions use the cache? Defaults to \code{TRUE}.
#' @param verbose Logical.
#' @return A \code{mart} object for use in biomart queries.
#' @keywords internal
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

#' Make chunked getBM query
#'
#' @param mart Biomart \code{mart} object.
#' @param values \code{character} vector of ids to be converted.
#' @param biom.attributes \code{character} vector. Biomart attributes, i.e., type of desired result(s); make sure query id type is included!
#' @param biom.filter \code{character} of length one. Name of biomart filter, i.e., type of query ids, defaults to "ensembl_gene_id".
#' @param use.cache (\code{logical}). Should \command{getBM()} use the cache? Defaults to \code{TRUE} as in the \command{getBM()} function and is passed on to that.
#' @param chunk.size \code{integer} of length one. Maximum number of IDs per BioMart query.
#'   Large ID lists are split into chunks of this size to avoid server timeouts.
#'   Set to \code{Inf} to disable chunking. Defaults to \code{500}.
#' @param verbose (\code{logical}). Should verbose output be written to the console?
#' @return A data frame with the retrieved information.
#' @keywords internal
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

