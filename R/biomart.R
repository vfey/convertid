# These functions are copied from the biomaRt package to avoid ':::' operators.

## from 'caching.R':

#' Add values to cache
#' @param bfc Object of class BiocFileCache, created by a call to
#' BiocFileCache::BiocFileCache()
#' @param result character; name of the file written to chache
#' @param hash unique hash representing a query.
#' @return Invisibly, \code{TRUE} if the value was added to the cache and
#' \code{FALSE} if an entry with the same hash already existed.
#' @keywords internal
.addToCache <- function(bfc, result, hash) {

  if(!dir.exists(.biomartCacheLocation()))
    dir.create(.biomartCacheLocation())

  ## write our file to the biomart cache location directly
  tf <- tempfile(tmpdir = .biomartCacheLocation())
  saveRDS(result, file = tf)

  ## check once more that there isn't an entry with this hash
  ## if its free add our new file
  ## if there's a clash don't add anything and tidy up
  if(!.checkInCache(bfc, hash = hash)) {
    bfcadd(bfc, rname = hash, fpath = tf, action = "asis")
    res <- TRUE
  } else {
    file.remove(tf)
    res <- FALSE
  }
  return(invisible(res))
}

#' Read values from cache
#' @param bfc Object of class BiocFileCache, created by a call to
#' BiocFileCache::BiocFileCache()
#' @param hash unique hash representing a query.
#' @return The cached R object stored under \code{hash}.
#' @keywords internal
.readFromCache <- function(bfc, hash) {

  cache_hits <- bfcquery(bfc, hash, field = "rname")
  if(nrow(cache_hits) > 1) {
    stop("Multiple cache results found.",
         "\nPlease clear your cache by running biomartCacheClear()")
  } else {
    rid <- cache_hits$rid
    result <- readRDS( bfc[[ rid ]] )
    return(result)
  }
}

#' Check whether value in cache exists
#' @param bfc Object of class BiocFileCache, created by a call to
#' BiocFileCache::BiocFileCache()
#' @param hash unique hash representing a query.
#' @param verbose logical; should additional verbose output be printed? Not currently used.
#' @details This function returns TRUE if a record with the requested hash already
#' exists in the file cache, otherwise returns FALSE.
#' @return \code{TRUE} if a record with the requested hash exists in the file
#' cache, otherwise \code{FALSE}.
#' @keywords internal
.checkInCache <- function(bfc, hash, verbose = FALSE) {
  res <- bfcquery(bfc, query = hash, field = "rname")
  as.logical(nrow(res))
}

#' Unexported functions
#' Report the \emph{\code{biomaRt}} cache location
#' @description \command{.biomartCacheLocation()} returns the path held in the
#' BIOMART_CACHE environment variable, falling back to the \code{rappdirs}
#' user cache directory for the \emph{\code{biomaRt}} application when that
#' variable is unset. Copied from the \emph{\code{biomaRt}} package to avoid a
#' ':::' operator.
#' @return (\code{character}) of length one. The cache directory path.
#' @seealso \code{\link[rappdirs]{user_cache_dir}}
#' @keywords internal
.biomartCacheLocation <- function() {
  Sys.getenv(x = "BIOMART_CACHE",
             unset = rappdirs::user_cache_dir(appname="biomaRt"))
}

## from 'ensembl_ssl_settings.R'.

#' Unexported functions
#' Probe the Ensembl hosts over HTTPS
#' @description \command{.test_ensembl()} issues two short \command{GET()}
#' requests against the main and US East Ensembl sites. It exists so that
#' \command{.checkEnsemblSSL()} can provoke, and then inspect, whatever SSL
#' error the local system produces. Copied from the \emph{\code{biomaRt}}
#' package to avoid a ':::' operator.
#' @return The response of the second request, invisibly; the function is
#' called for the error it raises, not for its value.
#' @seealso \code{\link[httr]{GET}}
#' @keywords internal
.test_ensembl <- function() {
  httr::GET("https://www.ensembl.org/index.html", timeout(5), set_cookies(redirect_mirror = "no"))
  httr::GET("https://useast.ensembl.org/index.html", timeout(5), set_cookies(redirect_mirror = "no"))
}

#' Unexported functions
#' Determine the CURL SSL options an Ensembl connection needs
#' @description \command{.checkEnsemblSSL()} probes the Ensembl hosts and, when
#' the connection fails, inspects the error to decide which \command{httr}
#' configuration works around it: a lowered cipher security level for the
#' 'sslv3 alert handshake failure' seen on Ubuntu 20.04 and relatives, or
#' disabled peer verification for missing issuer certificates. A timeout, or an
#' error it does not recognise, ends the probing. The function is a modified
#' version of \code{.checkEnsemblSSL} from the \emph{\code{biomaRt}} package.
#' @return A \code{list} of \code{request} objects holding the CURL options
#' that were needed; empty when the plain connection already works.
#' @seealso \code{\link[httr]{config}}, \code{\link{.test_ensembl}}
#' @keywords internal
.checkEnsemblSSL <- function() {

  ## we'll modify the global httr setting in this function
  ## get the current settings so we can reset it later
  prev_config <- getOption("httr_config")
  on.exit(options("httr_config" = prev_config))

  ensembl_config <- list()
  i <- 1
  test <- try(.test_ensembl(), silent = TRUE)
  while(is(test, "try-error")) {

    if(grepl(test[1], ## This address problems with Ubuntu 20.04 et al and the Ensembl https certificates
             pattern = "sslv3 alert handshake failure")) {
      ensembl_config[[i]] <- httr::config(ssl_cipher_list = "DEFAULT@SECLEVEL=1")
    } else if (grepl(x = test[1], ## two reported error messages solved with the same setting
                     pattern = "(unable to get local issuer certificate)|(server certificate verification failed)")) {
      ensembl_config[[i]] <- new_config <- httr::config(ssl_verifypeer = FALSE)
    } else if (grepl(x = test[1], ## We end up here if the test timed out
                     pattern = "Timeout was reached")) {
      ## Time out is unfortunate, but lets not inform the user since it might not be a problem.
      ## Quit the testing and proceed
      break;
    } else {
      message("Possible SSL connectivity problems detected.\n",
              "Please report this issue at https://github.com/grimbough/biomaRt/issues\n",
              test[1])
      ## We can't fix this, so just quit the tests
      break;
    }
    httr::set_config(ensembl_config[[i]], override = FALSE)
    i <- i + 1
    test <- try(.test_ensembl(), silent = TRUE)
  }
  return(ensembl_config)
}
