#' Unexported functions
#' Set the location for the biomaRt cache
#' @description \command{.setBiomaRtCacheLocation()} determines and sets the cache location
#' used by the functions in the 'biomaRt' package and defined in the BIOMART_CACHE
#' environment variable.
#' If the system default cache location is under the user's home directory a sub-folder biomaRt is used.
#' If not and 'x' is an empty string (the default) and the variable is unset it will be set to
#' "~/.caches/biomaRt", which is the default location in Linux. The tilde will be
#' expanded to match the user's home directory.
#' If the variable has been set externally to a location under the user's home that path is returned.
#' @param x (\code{character}). Optional user-defined path used as cache location. Defaults to "".
#' @return A character vector with the BIOMART_CACHE location.
#' @seealso \code{\link[rappdirs]{user_cache_dir}}, \code{\link[base]{Sys.getenv}}
#' @examples
#' \dontrun{.setBiomaRtCacheLocation()}
#' @keywords internal
.setBiomaRtCacheLocation <-
  function(x = "") {
    system.cache <- path.expand(rappdirs::user_cache_dir())
    user <- Sys.getenv("USER")
    if (!length(grep(user, path.expand(Sys.getenv(x = "BIOMART_CACHE"))))) {
      Sys.unsetenv("BIOMART_CACHE")
    }
    current.cache <- path.expand(Sys.getenv(x = "BIOMART_CACHE"))
    if (length(grep(user, system.cache)) && current.cache == "") {
      Sys.getenv(x = "BIOMART_CACHE", unset = rappdirs::user_cache_dir(appname = "biomaRt"))
    } else if (current.cache == "") {
      if (x != "") Sys.setenv(BIOMART_CACHE = x)
      Sys.getenv(x = "BIOMART_CACHE", unset = path.expand("~/.caches/biomaRt"))
    } else {
      current.cache
    }
  }
