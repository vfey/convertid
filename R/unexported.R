#' Unexported functions
#' Set the location for the biomaRt cache
#' @description \command{.setBiomaRtCacheLocation()} attempts to set the cache location
#' used by the functions in the 'biomaRt' package and defined in the BIOMART_CACHE
#' environment variable.
#' If that variable is set and the defined location exists and is writable nothing is done.
#' If the system default cache location exists and is writable a sub-folder \code{app} is used (and created if necessary).
#' If not the function attempts to create \code{cache.dir} and a sub-folder \code{app}, whether or not those exist, already.
#' If all of the above fail the function attempts to create \code{file.path(tempdir(), app)}. If tat fails, too,
#' an exception is thrown.
#' @param cache.dir (\code{character}). Optional user-defined path used as cache parent directory. Defaults to \code{rappdirs::user_cache_dir()}.
#' @param app (\code{character}). Optional application-specific cache sub-directory, i.e., the actually used cache location. Defaults to "biomaRt".
#' @return A path name with the BIOMART_CACHE location.
#' @seealso \code{\link[rappdirs]{user_cache_dir}}, \code{\link[BiocFileCache]{BiocFileCache}}
#' @examples
#' \dontrun{.setCacheLocation()}
#' @keywords internal
.setCacheLocation <-
  function(cache.dir = rappdirs::user_cache_dir(), app = "biomaRt") {
    system.cache <- rappdirs::user_cache_dir()
    current.cache <- Sys.getenv(x = "BIOMART_CACHE")
    if (.cache.writable(current.cache)) {
      cache.path <- current.cache
    } else if (.cache.writable(system.cache)) {
      cache.path <- rappdirs::user_cache_dir(app)
    } else {
      if (app == basename(cache.dir)) {
        cache.path <- cache.dir
      } else {
        cache.path <- file.path(cache.dir, app)
      }
    }
    path.ok <- .create.cache(cache.path)
    if (!path.ok)
      tmp.ok <- .create.cache(file.path(tempdir(), app))
    if (!tmp.ok) {
      stop("Could not create cache directory. Make sure the location in 'cache.dir' is writable.")
    }
    Sys.setenv(BIOMART_CACHE = cache.path)
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
  function(cache.path = BiocFileCache::getBFCOption("CACHE")) {
    cache.ok <- try(BiocFileCache::BiocFileCache(cache.path, ask = FALSE))
    return(!is(cache.ok, "try-error"))
  }
