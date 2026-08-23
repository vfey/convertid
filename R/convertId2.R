#' @title Convert Gene IDs Between Each Other and Fetch Annotations from Biomart
#' @name convertid
#' @description Gene Symbols, Ensembl Gene IDs or Entrez Gene IDs are converted
#'     using the Bimap interface in 'AnnotationDbi' in convertId2() for the most
#'     common use cases in data analysis. The main function in the package is
#'     convert.bm() which queries Biomart using the full capacity of the API
#'     provided through the 'biomaRt' package. Presets and defaults are provided
#'     for convenience but all "marts", "filters" and "attributes" can be set by
#'     the user. Function convert.alias() converts Gene Symbols to Aliases and
#'     vice versa and function likely_symbol() attempts to determine the most
#'     likely current Gene Symbol.
#' @author Vidal Fey <vidal.fey@gmail.com>
#' Maintainer: Vidal Fey <vidal.fey@gmail.com>
#' @details \tabular{ll}{
#' Package: \tab convertid\cr
#' Type: \tab Package\cr
#' Initial version: \tab 0.1-0\cr
#' Created: \tab 2021-08-18\cr
#' License: \tab GPL-3\cr
#' LazyLoad: \tab yes\cr
#' }
#'
#' @keywords package
#' @keywords internal
"_PACKAGE"
#' @import AnnotationDbi
#' @import plyr
#' @import stringr
#' @import biomaRt
#' @import xml2
#' @importFrom methods is
#' @importFrom assertthat assert_that
#' @importFrom stats na.omit setNames
#' @importFrom utils read.delim
#' @importFrom rappdirs user_cache_dir
#' @importFrom BiocFileCache BiocFileCache bfcadd bfcquery
#' @importFrom httr config set_config set_cookies GET timeout
NULL

#' Unexported functions
#' Strict single-value lookup in an AnnotationDbi map
#' @description \command{.safe_get1()} returns the value stored under
#' \code{key}, but only when that value is unique. A missing key, or a key
#' mapping to more than one value, yields \code{NA_character_}. It is used for
#' the intermediate Entrez hop in \command{convertId2()}, where an ambiguous
#' result has to be discarded to preserve the one-to-one correspondence
#' between input and output, independently of \code{multi2NA}.
#' @param key (\code{character}) of length one. The identifier to look up.
#' @param env An AnnotationDbi \code{Bimap} or environment to look in.
#' @return A \code{character} of length one, or \code{NA_character_} when the
#' key is absent or ambiguous.
#' @keywords internal
.safe_get1 <- function(key, env) {
  if (AnnotationDbi::exists(key, envir = env)) {
    val <- get(key, envir = env)
    if (length(val) > 1L) NA_character_ else val
  } else {
    NA_character_
  }
}

#' Unexported functions
#' Resolve a lookup result to a single character value
#' @description \command{.resolve_multi()} collapses the result of a terminal
#' identifier lookup to one element. Nothing found becomes
#' \code{NA_character_}; a one-to-many mapping becomes either
#' \code{NA_character_} or a \code{" /// "}-separated string, depending on
#' \code{multi2NA}.
#' @param val The looked-up value; a vector or list of any length.
#' @param multi2NA (\code{logical}). Should a one-to-many mapping yield
#' \code{NA_character_} instead of a collapsed string? Defaults to
#' \code{FALSE}.
#' @return A \code{character} of length one, possibly \code{NA_character_}.
#' @keywords internal
.resolve_multi <- function(val, multi2NA = FALSE) {
  if (length(val) == 0L) return(NA_character_)
  if (length(val) > 1L) {
    return(if (multi2NA) NA_character_ else paste(val, collapse = " /// "))
  }
  if (is.na(val[[1L]])) NA_character_ else as.character(val[[1L]])
}

#' Unexported functions
#' Terminal single-key lookup in an AnnotationDbi map
#' @description \command{.get_term()} looks \code{key} up and hands the result
#' to \command{.resolve_multi()}. A missing key yields \code{NA_character_}.
#' This is the terminal step of a conversion, as opposed to the strict
#' intermediate hop performed by \command{.safe_get1()}.
#' @param key (\code{character}) of length one. The identifier to look up.
#' @param env An AnnotationDbi \code{Bimap} or environment to look in.
#' @param multi2NA (\code{logical}). Passed to \command{.resolve_multi()}.
#' Defaults to \code{FALSE}.
#' @return A \code{character} of length one, possibly \code{NA_character_}.
#' @seealso \code{\link{.resolve_multi}}, \code{\link{.safe_get1}}
#' @keywords internal
.get_term <- function(key, env, multi2NA = FALSE) {
  if (AnnotationDbi::exists(key, envir = env)) {
    .resolve_multi(get(key, envir = env), multi2NA = multi2NA)
  } else {
    NA_character_
  }
}

#' Convert Gene IDs Between Ensembl, Symbol and Entrez
#' @description \command{convertId2()} is a fast 1-to-1 gene identifier converter
#'     using AnnotationDbi organism packages. It converts between Ensembl gene IDs,
#'     gene symbols and Entrez gene IDs. The output vector is always the same
#'     length and in the same order as the input vector, with unresolved entries
#'     returned as \code{NA}.
#' @param id (\code{character}). Vector of gene identifiers to convert. Can be
#'     Ensembl gene IDs (e.g. \code{"ENSG00000075624"}), gene symbols
#'     (e.g. \code{"ACTB"}), or Entrez gene IDs (e.g. \code{"60"}).
#'     Input type is detected automatically from the first non-NA element.
#' @param species (\code{character}). One of \code{"Human"} or \code{"Mouse"}.
#'     Defaults to \code{"Human"}.
#' @param output (\code{character}). One of \code{"auto"}, \code{"symbol"},
#'     \code{"ensembl"}, or \code{"entrez"}. Controls the return type:
#'     \describe{
#'       \item{\code{"auto"}}{Automatic: Ensembl input returns symbols, symbol
#'         input returns Ensembl IDs, Entrez input returns symbols.}
#'       \item{\code{"symbol"}}{Always return gene symbols regardless of input type.}
#'       \item{\code{"ensembl"}}{Always return Ensembl gene IDs regardless of input
#'         type. Useful when converting Entrez IDs to Ensembl IDs for downstream
#'         processing.}
#'       \item{\code{"entrez"}}{Always return Entrez gene IDs regardless of input
#'         type. Useful for offline, network-independent retrieval of Entrez IDs
#'         (e.g. for KEGG pathway enrichment) as a fast alternative to a BioMart
#'         query.}
#'     }
#'     Defaults to \code{"auto"}. Note that \code{"auto"} never yields Entrez IDs;
#'     request them explicitly with \code{output = "entrez"}.
#' @param multi2NA (\code{logical}). Controls how one-to-many mappings are
#'     resolved in the returned vector. When a query ID maps to more than one
#'     target ID (e.g. one gene symbol resolving to several Ensembl or Entrez
#'     IDs), \code{multi2NA = FALSE} (the default) collapses all hits into a
#'     single string separated by \code{" /// "}, whereas \code{multi2NA = TRUE}
#'     returns \code{NA} for that entry. Set \code{multi2NA = TRUE} when the
#'     result must contain at most one ID per input for downstream use (at the
#'     cost of possible \code{NA}s). This setting applies to all output types.
#'     Note that the \emph{intermediate} Entrez step used for
#'     Ensembl-to-symbol and symbol-to-Ensembl conversion is always resolved
#'     strictly (an ambiguous intermediate is dropped to \code{NA}) regardless
#'     of \code{multi2NA}. Defaults to \code{FALSE}.
#' @details
#' Conversion routes, depending on input type and \code{output}:
#' \itemize{
#'   \item Ensembl \eqn{\rightarrow} Entrez \eqn{\rightarrow} symbol
#'   \item Ensembl \eqn{\rightarrow} Entrez (direct)
#'   \item Symbol \eqn{\rightarrow} Entrez \eqn{\rightarrow} Ensembl
#'   \item Symbol \eqn{\rightarrow} Entrez (direct)
#'   \item Entrez \eqn{\rightarrow} symbol (direct)
#'   \item Entrez \eqn{\rightarrow} Ensembl (direct)
#' }
#' When input and output type coincide (e.g. Ensembl input with
#' \code{output = "ensembl"}) the input is returned unchanged.
#'
#' Entries are returned as \code{NA} when the input ID is not found in the
#' annotation database. One-to-many mappings in the final lookup are governed
#' by \code{multi2NA}: collapsed into a \code{" /// "}-separated string by
#' default, or returned as \code{NA} when \code{multi2NA = TRUE}. The
#' intermediate Entrez step used for Ensembl-to-symbol and symbol-to-Ensembl
#' conversion is always resolved strictly: an input mapping to more than one
#' Entrez gene ID at that step is dropped to \code{NA} to preserve the 1-to-1
#' correspondence between input and output vectors.
#' Input type is detected automatically from the first non-NA element of
#' \code{id}: IDs matching the species Ensembl prefix (\code{ENSG} for Human,
#' \code{ENSMU} for Mouse) are treated as Ensembl gene IDs; purely numeric
#' strings are treated as Entrez gene IDs; all others are treated as gene
#' symbols. All elements of \code{id} are assumed to be of the same type.
#' The function is limited to Human and Mouse annotations and is provided mainly
#' as fast conversion mechanism for the most common use cases in data analysis.
#' @return A named character vector of the same length and order as \code{id},
#'     named by the input IDs. Entries that could not be converted are
#'     \code{NA}.
#' @seealso \code{\link{convert.bm}} for BioMart-based conversion which
#'     returns richer annotations but does not guarantee output length or order.
#' @examples
#' \dontrun{
#' # Ensembl -> symbol (auto)
#' convertId2("ENSG00000075624")
#' convertId2(c("ENSG00000075624", "ENSG00000111640"))
#'
#' # Symbol -> Ensembl (auto)
#' convertId2("ACTB")
#' convertId2(c("ACTB", "GAPDH"))
#'
#' # Entrez -> symbol (auto)
#' convertId2("60")
#' convertId2(c("60", "2597"))
#'
#' # Entrez -> Ensembl (explicit output)
#' convertId2("60", output = "ensembl")
#' convertId2(c("60", "2597"), output = "ensembl")
#'
#' # Ensembl -> Entrez (explicit output)
#' convertId2("ENSG00000075624", output = "entrez")
#' convertId2(c("ENSG00000075624", "ENSG00000111640"), output = "entrez")
#'
#' # Symbol -> Entrez (explicit output)
#' convertId2(c("ACTB", "GAPDH"), output = "entrez")
#'
#' # Return NA instead of a "///"-collapsed string for one-to-many mappings
#' convertId2(c("ACTB", "GAPDH"), output = "entrez", multi2NA = TRUE)
#' }
#' @export
convertId2 <- function(id, species = c("Human", "Mouse"),
                       output = c("auto", "symbol", "ensembl", "entrez"),
                       multi2NA = FALSE)
{
  species <- match.arg(species)
  output  <- match.arg(output)

  # ---------------------------------------------------------------------------
  # Load organism package environments
  # ---------------------------------------------------------------------------
  if (species == "Human") {
    if (!requireNamespace("org.Hs.eg.db", quietly = TRUE))
      stop(
        "The Bioconductor package 'org.Hs.eg.db' is required for this function.\n",
        "Install it via BiocManager::install('org.Hs.eg.db').",
        call. = FALSE
      )
    orgdb    <- getNamespace("org.Hs.eg.db")
    ensg2eg  <- orgdb[["org.Hs.egENSEMBL2EG"]]
    sym      <- orgdb[["org.Hs.egSYMBOL"]]
    sym2eg   <- orgdb[["org.Hs.egSYMBOL2EG"]]
    ensg     <- orgdb[["org.Hs.egENSEMBL"]]
    ensg_pfx <- "ENSG"
  }
  if (species == "Mouse") {
    if (!requireNamespace("org.Mm.eg.db", quietly = TRUE))
      stop(
        "The Bioconductor package 'org.Mm.eg.db' is required for this function.\n",
        "Install it via BiocManager::install('org.Mm.eg.db').",
        call. = FALSE
      )
    orgdb    <- getNamespace("org.Mm.eg.db")
    ensg2eg  <- orgdb[["org.Mm.egENSEMBL2EG"]]
    sym      <- orgdb[["org.Mm.egSYMBOL"]]
    sym2eg   <- orgdb[["org.Mm.egSYMBOL2EG"]]
    ensg     <- orgdb[["org.Mm.egENSEMBL"]]
    ensg_pfx <- "ENSMU"
  }

  # ---------------------------------------------------------------------------
  # Detect input type from first non-NA element
  # ---------------------------------------------------------------------------
  first_id   <- id[!is.na(id)][1L]
  is_ensembl <- length(grep(ensg_pfx, first_id)) > 0L
  is_entrez  <- !is_ensembl && grepl("^[0-9]+$", first_id)
  # otherwise: symbol

  # Resolve effective output type under "auto"
  eff_output <- if (output == "auto") {
    if (is_ensembl || is_entrez) "symbol" else "ensembl"
  } else {
    output
  }

  # ---------------------------------------------------------------------------
  # Single-element fast path
  # ---------------------------------------------------------------------------
  if (length(id) == 1L) {
    result <- if (is_ensembl) {
      if (eff_output == "symbol") {
        eg <- .safe_get1(id, ensg2eg)                 # intermediate hop (strict)
        if (is.na(eg)) NA_character_ else .get_term(eg, sym, multi2NA)
      } else if (eff_output == "entrez") {
        .get_term(id, ensg2eg, multi2NA)                        # ENSG -> Entrez (direct)
      } else {
        id                                           # ensembl -> ensembl: as-is
      }
    } else if (is_entrez) {
      if (eff_output == "symbol") {
        .get_term(id, sym, multi2NA)                            # Entrez -> symbol (direct)
      } else if (eff_output == "entrez") {
        id                                           # entrez -> entrez: as-is
      } else {
        .get_term(id, ensg, multi2NA)                           # Entrez -> ENSG (direct)
      }
    } else {
      # symbol input
      if (eff_output == "ensembl") {
        eg <- .safe_get1(id, sym2eg)                  # intermediate hop (strict)
        if (is.na(eg)) NA_character_ else .get_term(eg, ensg, multi2NA)
      } else if (eff_output == "entrez") {
        .get_term(id, sym2eg, multi2NA)                         # Symbol -> Entrez (direct)
      } else {
        id                                           # symbol -> symbol: as-is
      }
    }
    names(result) <- id
    return(result)
  }

  # ---------------------------------------------------------------------------
  # Batch path
  # ---------------------------------------------------------------------------
  result <- if (is_ensembl) {
    if (eff_output == "symbol") {
      # ENSG -> Entrez -> symbol; intermediate multi/NA Entrez dropped to "---"
      entrez <- mget(id, envir = ensg2eg, ifnotfound = NA)
      entrez <- sapply(entrez, function(x) {
        if (length(x) > 1L || is.na(x)) "---" else x
      })
      hugo <- mget(entrez, envir = sym, ifnotfound = NA)
      sapply(hugo, .resolve_multi, multi2NA = multi2NA)
    } else if (eff_output == "entrez") {
      # ENSG -> Entrez (direct)
      sapply(mget(id, envir = ensg2eg, ifnotfound = NA), .resolve_multi, multi2NA = multi2NA)
    } else {
      setNames(id, id)   # ensembl -> ensembl
    }

  } else if (is_entrez) {
    if (eff_output == "symbol") {
      # Entrez -> symbol (direct)
      sapply(mget(id, envir = sym, ifnotfound = NA), .resolve_multi, multi2NA = multi2NA)
    } else if (eff_output == "entrez") {
      setNames(id, id)   # entrez -> entrez
    } else {
      # Entrez -> ENSG (direct)
      sapply(mget(id, envir = ensg, ifnotfound = NA), .resolve_multi, multi2NA = multi2NA)
    }

  } else {
    # Symbol input
    if (eff_output == "ensembl") {
      # symbol -> Entrez -> ENSG; intermediate multi/NA Entrez dropped to "---"
      entrez <- mget(id, envir = sym2eg, ifnotfound = NA)
      entrez <- sapply(entrez, function(x) {
        if (length(x) > 1L || is.na(x)) "---" else x
      })
      result_ensg <- mget(entrez, envir = ensg, ifnotfound = NA)
      sapply(result_ensg, .resolve_multi, multi2NA = multi2NA)
    } else if (eff_output == "entrez") {
      # symbol -> Entrez (direct)
      sapply(mget(id, envir = sym2eg, ifnotfound = NA), .resolve_multi, multi2NA = multi2NA)
    } else {
      setNames(id, id)   # symbol -> symbol
    }
  }

  names(result) <- id
  result
}

#' Convert Symbols to Aliases and Vice Versa.
#' @description \command{convert.alias()} attempts to find all possible symbol-alias combinations for a given gene symbol, i.e.,
#'     it assumes the input ID to be either an Alias or a Symbol and performs multiple queries to find all possible
#'     counterparts. The input IDs are converted to title and upper case before querying and all possibilities are tested.
#'     There are species presets for Human and Mouse annotations.
#' @param id (\code{character}). Vector of gene symbols.
#' @param species (\code{character}). One of "Human" and "Mouse". Defaults to "Human".
#' @param db (\code{AnnotationDb object}). Annotation package object.
#' @return A \code{data.frame} with two columns:
#' \tabular{ll}{
#' \tab 'SYMBOL': The official gene symbol.\cr
#' \tab 'ALIAS': All possible aliases.\cr
#' }
#' @seealso \code{\link[AnnotationDbi]{select}}
#' @examples
#' \donttest{
#' if (requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
#'   convert.alias("TRPV4")
#' }
#' }
#' @export
convert.alias <-
  function(id, species = c("Human", "Mouse"), db = NULL)
  {
    if (missing(id))
      stop("Need input ID vector!")
    if (is.null(db)) {
      species <- match.arg(species)
      orgpkg <- switch(species,
                       Human = "org.Hs.eg.db",
                       Mouse = "org.Mm.eg.db"
      )
      if (!requireNamespace(orgpkg, quietly = TRUE)) {
        stop(
          sprintf("The Bioconductor package '%s' is required for this function.\n", orgpkg),
          sprintf("Install it via BiocManager::install('%s').", orgpkg),
          call. = FALSE
        )
      }
      # Get the namespace without using ::
      db <- getNamespace(orgpkg)[[orgpkg]]
    }
    syms <- plyr::ldply(id, function(i) {
      i1 <- stringr::str_to_title(i)
      i2 <- stringr::str_to_upper(i)
      kdf <- plyr::ldply(c("ALIAS", "SYMBOL"), function(k) {
        validkeys <- keys(db, k)
        plyr::ldply(c(i1, i2), function(x) {
          if (any(validkeys %in% x))
            suppressMessages(AnnotationDbi::select(db, keys=x, columns=c("SYMBOL","ALIAS"), keytype=k))
          else
            data.frame(ALIAS=NA_character_, SYMBOL=NA_character_)
        })
      })
      if (all(is.na(kdf[[1]]))) {
        data.frame(ALIAS=i2, SYMBOL=NA_character_)
      } else {
        kdf <- na.omit(kdf)
        plyr::ddply(kdf, "SYMBOL", function(s) {
          suppressMessages(AnnotationDbi::select(db, keys=s$SYMBOL, columns=c("SYMBOL","ALIAS"), keytype="SYMBOL"))
        })
      }
    })
    syms <- syms[!duplicated(syms$ALIAS), ]
    return(syms)
  }

#' Retrieve Additional Annotations from Biomart
#' @description \command{convert.bm()} is a wrapper for \command{get.bm()} which in turn makes use of \command{getBM()} from the \emph{biomaRt} package.
#' It takes a matrix or data frame with the IDs to be converted in one column or as row names as input and returns a data frame with additional
#' annotations after cleaning the fetched annotations and merging them with the input data frame.
#' @param dat \code{matrix} or \code{data.frame}. Matrix or data frame with the ids to be converted in a column or as row names.
#' @param id \code{character}. Name of the column with the ids to be converted, special name "rownames" will use the row names.
#' @param biom.data.set \code{character} of length one. Biomart data set to use.
#' @param biom.mart \code{character} vector. Biomart to use (uses the first element of the vector), defaults to "ensembl".
#' @param host \code{character} of length one. Host URL.
#' @param biom.filter \code{character} of length one. Name of biomart filter, i.e., type of query ids, defaults to "ensembl_gene_id".
#' @param biom.attributes \code{character} vector. Biomart attributes, i.e., type of desired result(s);
#'   if \code{biom.filter} is missing from this it will be added internally as it is needed for merging query result and input data.
#' @param biom.cache \code{character}. Path name giving the location of the cache \command{getBM()} uses if \code{use.cache=TRUE}. Defaults to the value in the \emph{BIOMART_CACHE} environment variable.
#' @param use.cache (\code{logical}). Should \command{getBM()} use the cache? Defaults to \code{TRUE} as in the \command{getBM()} function and is passed on to that.
#' @param biomart.fallback \code{character} vector. Fallback host URLs to try if the primary
#'   \code{host} fails. Set to \code{NULL} to disable fallback. Defaults to Ensembl mirror sites.
#' @param chunk.size \code{integer} of length one. Maximum number of IDs per BioMart query.
#'   Large ID lists are split into chunks of this size to avoid server timeouts.
#'   Set to \code{Inf} to disable chunking. Defaults to \code{500}.
#' @param sym.col \code{character}. Name of the column in the query result with gene symbols.
#' @param rm.dups \code{logical}. Should duplicated input IDs (\option{biom.filter}) be removed from the result?
#' @param verbose (\code{logical}). Should verbose output be written to the console? Defaults to \code{FALSE}.
#' @details Wrapped around `get.bm`.
#' @return  A data frame with the retrieved information.
#' @author Vidal Fey
#' @seealso \command{\link[biomaRt]{getBM}}
#' @keywords utilities
#' @examples
#' \dontrun{
#' dat <- data.frame(ID=c("ENSG00000111199", "ENSG00000134121", "ENSG00000176102", "ENSG00000171611"))
#' bm <- convert.bm(dat)
#' bm
#' }
#' @export
convert.bm <-
  function(dat, id="ID",
           biom.data.set = c("human", "mouse"),
           biom.mart=c("ensembl", "mouse", "snp", "funcgen", "plants"),
           host="https://www.ensembl.org",
           biom.filter="ensembl_gene_id",
           biom.attributes=c("ensembl_gene_id","hgnc_symbol","description"),
           biom.cache = rappdirs::user_cache_dir("biomaRt"),
           use.cache = TRUE,
           biomart.fallback = c("https://useast.ensembl.org",
                                "https://uswest.ensembl.org",
                                "https://asia.ensembl.org"),
           chunk.size = 500L,
           sym.col="hgnc_symbol",
           rm.dups=FALSE,
           verbose = FALSE)
  {
    if (id=="row.names") {
      values <- rownames(dat)
    } else {
      values <- dat[[id]]
    }
    if (!biom.filter %in% biom.attributes) {
      biom.attributes <- c(biom.filter, biom.attributes)
    }
    biom.ids <- get.bm(values           = values,
                       biom.data.set    = biom.data.set,
                       biom.mart        = biom.mart,
                       host             = host,
                       biom.filter      = biom.filter,
                       biom.attributes  = biom.attributes,
                       biom.cache       = biom.cache,
                       use.cache        = use.cache,
                       biomart.fallback = biomart.fallback,
                       chunk.size       = chunk.size,
                       verbose          = verbose)
    gene.lab <- merge(biom.ids, dat, by.x=biom.filter, by.y=id, all.y=TRUE, all.x=FALSE, sort=TRUE)
    if (rm.dups) {
      if (verbose) message("  Removing ", length(which(duplicated(gene.lab[[biom.filter]]))), " duplicated row(s)...")
      gene.lab <- gene.lab[!duplicated(gene.lab[[biom.filter]]), ]
    }
    missing <- gene.lab[[sym.col]] == "" | is.na(gene.lab[[sym.col]])
    if (any(missing)) {
      if (verbose) message("  Replacing ", sum(missing),
                           " missing Gene Symbols by ", sQuote(biom.filter), "...")
      gene.lab[[sym.col]][missing] <- gene.lab[[biom.filter]][missing]
    }
    return(gene.lab)
  }

#' Convenience Function to Convert Ensembl Gene IDs to Gene Symbols
#' @description \command{todisp2()} uses Biomart by employing \command{get.bm()} to retrieve Gene Symbols for a set of Ensembl
#'     Gene IDs. It is mainly meant as a fast way to convert IDs in standard gene expression analysis output to Symbols,
#'     e.g., for visualisation, which is why the input ID type is hard-coded to ENSG IDs. If Biomart is not available
#'     the function can fall back to use \command{convertId2()} or a user-provided data frame with corresponding ENSG IDs and
#'     Symbols.
#' @param ensg (\code{character}). Vector of Ensemble Gene IDs. Other ID types are not yet supported.
#' @param lab (\code{data.frame}). A data frame with Ensembl Gene IDs as row names and Gene Symbols in the only column.
#' @param biomart (\code{logical}). Should Biomart be used? Defaults to \code{TRUE}.
#' @param biom.data.set \code{character} of length one. Biomart data set to use. Defaults to 'hsapiens_gene_ensembl'
#' @param biom.mart \code{character} vector. Biomart to use (uses the first element of the vector), defaults to "ensembl".
#' @param host \code{character} of length one. Host URL.
#' @param biom.filter \code{character} of length one. Name of biomart filter, i.e., type of query ids, defaults to "ensembl_gene_id".
#' @param biom.attributes \code{character} vector. Biomart attributes, i.e., type of desired result(s); make sure query id type is included!
#' @param biom.cache \code{character}. Path name giving the location of the cache \command{getBM()} uses if \code{use.cache=TRUE}. Defaults to the value in the \emph{BIOMART_CACHE} environment variable.
#' @param use.cache (\code{logical}). Should \command{getBM()} use the cache? Defaults to \code{TRUE} as in the \command{getBM()} function and is passed on to that.
#' @param keep.original (\code{logical}). Should the order and length of the input vector be preserved, i.e., should also IDs missing after conversion be kept? Defaults to \code{TRUE}.
#' @param verbose (\code{logical}). Should verbose output be written to the console? Defaults to \code{FALSE}.
#' @return A character vector of Gene Symbols.
#' @seealso \command{\link[convertid]{get.bm}}
#' @examples
#' \dontrun{
#' val <- c("ENSG00000111199", "ENSG00000134121", "ENSG00000176102", "ENSG00000171611")
#' sym <- todisp2(val)
#' sym
#' }
#' @keywords utilities
#' @export
todisp2 <- function(ensg,
                    lab=NULL,
                    biomart=TRUE,
                    biom.data.set = "hsapiens_gene_ensembl",
                    biom.mart = "ensembl",
                    host = "https://www.ensembl.org",
                    biom.filter = "ensembl_gene_id",
                    biom.attributes = c("ensembl_gene_id", "hgnc_symbol"),
                    biom.cache = rappdirs::user_cache_dir("biomaRt"),
                    use.cache = TRUE,
                    keep.original = TRUE,
                    verbose = FALSE)
{
  if (biomart) {
    if (!length(grep("^ENS[A-Z]{0,}[0-9]{11}", ensg[1]))) {
      if (verbose) message("    Input is not Ensembl Gene IDs. Doing nothing.")
      return(ensg)
    }
    sym <- get.bm(ensg, biom.data.set=biom.data.set, biom.mart=biom.mart, host=host,
                  biom.filter=biom.filter, biom.attributes=biom.attributes,
                  biom.cache = biom.cache, use.cache = use.cache,
                  verbose = verbose)
  } else if(!is.null(lab)) {
    if (verbose) message("  Using input data frame for ID conversion...")
    sym <- data.frame(ensembl_gene_id=rownames(lab), hgnc_symbol=lab[, 1], stringsAsFactors=FALSE)
  } else {
    if (verbose) message("  Using 'AnnotationDbi framework for ID conversion...")
    sym <- convertId2(ensg)
    if (length(sym) == 1 && is.na(sym)) {
      return(ensg)
    } else {
      if (length(grep("^ENS[A-Z]{0,}[0-9]{11}", na.omit(sym)[1]))) {
        if (verbose) message("    Input was Gene Symbol")
        return(names(sym))
      } else if(length(grep("^ENS[A-Z]{0,}[0-9]{11}", names(sym)[1]))) {
        if (verbose) message("    Input was Ensemble Gene ID")
        sym <- data.frame(ensembl_gene_id=names(sym), hgnc_symbol=as.character(sym), stringsAsFactors=FALSE)
      } else {
        return(ensg)
      }
    }
  }
  if (verbose) message("    Merging input IDs and converted IDs...")
  if (keep.original) {
    gene.lab <- merge(data.frame(ensembl_gene_id=ensg, stringsAsFactors=FALSE), sym, by="ensembl_gene_id", all.x = TRUE, sort=FALSE)
    gene.lab <- gene.lab[match(ensg, gene.lab$ensembl_gene_id), ]
  } else {
    gene.lab <- merge(data.frame(ensembl_gene_id=ensg, stringsAsFactors=FALSE), sym, by="ensembl_gene_id", sort=FALSE)
  }
  if (verbose) message("done")
  if (any(gene.lab$hgnc_symbol=="" | is.na(gene.lab$hgnc_symbol))) {
    if (verbose) message("    Replacing ", length(which(gene.lab$hgnc_symbol=="" | is.na(gene.lab$hgnc_symbol))), " missing Gene Symbol(s) by Ensembl IDs...")
    replace <- gene.lab$hgnc_symbol=="" | is.na(gene.lab$hgnc_symbol)
    gene.lab$hgnc_symbol[replace] <- gene.lab$ensembl_gene_id[replace]
  }
  return(gene.lab$hgnc_symbol)
}
