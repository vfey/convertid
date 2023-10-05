#' Unexported functions
#' Set the location for the biomaRt cache
#' @description \command{.setBiomaRtCacheLocation()} attempts to set the cache location
#' used by the functions in the 'biomaRt' package and defined in the BIOMART_CACHE
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
#' Create and initialise a BioC cache at a given location.
#' @description \command{.create.cache()} uses \command{BiocFileCache::BiocFileCache()} to set up
#' a BioC file cache at a location defined by the user. See the \command{BiocFileCache} package
#' for more information.
#' @param cache.path (\code{character}). The path to use for the cached files.
#' @return TRUE if the location was successfully set up, FALSE if not.
#' @seealso \code{\link[BiocFileCache]{BiocFileCache}}
#' @examples
#' \dontrun{.create.cache(rappdirs::user_cache_dir("biomaRt"))}
#' @keywords internal
.create.cache <-
  function(cache.path = rappdirs::user_cache_dir("biomaRt")) {
    assertthat::assert_that(is.character(cache.path), length(cache.path) == 1L, !is.na(cache.path))
    if (!file.exists(cache.path)) {
      dir.create(cache.path, recursive = TRUE, showWarnings = FALSE)
    }
    .cache.writable(cache.path)
  }
