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
#' @param output (\code{character}). One of \code{"auto"}, \code{"symbol"}, or
#'     \code{"ensembl"}. Controls the return type:
#'     \describe{
#'       \item{\code{"auto"}}{Automatic: Ensembl input returns symbols, symbol
#'         input returns Ensembl IDs, Entrez input returns symbols.}
#'       \item{\code{"symbol"}}{Always return gene symbols regardless of input type.}
#'       \item{\code{"ensembl"}}{Always return Ensembl gene IDs regardless of input
#'         type. Useful when converting Entrez IDs to Ensembl IDs for downstream
#'         processing.}
#'     }
#'     Defaults to \code{"auto"}.
#' @details
#' Conversion is performed via Entrez gene IDs as an intermediate step:
#' \itemize{
#'   \item Ensembl \eqn{\rightarrow} Entrez \eqn{\rightarrow} symbol
#'   \item Symbol \eqn{\rightarrow} Entrez \eqn{\rightarrow} Ensembl
#'   \item Entrez \eqn{\rightarrow} symbol (direct)
#'   \item Entrez \eqn{\rightarrow} Ensembl (direct)
#' }
#' Entries are returned as \code{NA} in two cases:
#' \itemize{
#'   \item The input ID is not found in the annotation database.
#'   \item The input ID maps to more than one Entrez gene ID. Such ambiguous
#'     mappings are discarded rather than returning multiple values, in order
#'     to strictly preserve the 1-to-1 correspondence between input and output
#'     vectors.
#' }
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
#' }
#' @export
convertId2 <- function(id, species = c("Human", "Mouse"),
                       output = c("auto", "symbol", "ensembl"))
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
  # Helper: safe single-value lookup; returns NA on miss or if length > 1
  # (length > 1 = ambiguous multi-Entrez mapping, discarded to preserve 1-to-1)
  # ---------------------------------------------------------------------------
  safe_get1 <- function(key, env) {
    if (AnnotationDbi::exists(key, envir = env)) {
      val <- get(key, envir = env)
      if (length(val) > 1L) NA_character_ else val
    } else {
      NA_character_
    }
  }

  # Helper: multi-value lookup collapsed to /// string; NA on miss
  safe_get_multi <- function(key, env) {
    if (AnnotationDbi::exists(key, envir = env)) {
      paste(get(key, envir = env), collapse = " /// ")
    } else {
      NA_character_
    }
  }

  # ---------------------------------------------------------------------------
  # Single-element fast path
  # ---------------------------------------------------------------------------
  if (length(id) == 1L) {
    result <- if (is_ensembl) {
      if (eff_output == "symbol") {
        eg <- safe_get1(id, ensg2eg)
        if (is.na(eg)) NA_character_ else safe_get1(eg, sym)
      } else {
        id   # ensembl -> ensembl: return as-is
      }
    } else if (is_entrez) {
      if (eff_output == "symbol") {
        safe_get1(id, sym)
      } else {
        safe_get_multi(id, ensg)
      }
    } else {
      # symbol input
      if (eff_output == "ensembl") {
        eg <- safe_get1(id, sym2eg)
        if (is.na(eg)) NA_character_ else safe_get_multi(eg, ensg)
      } else {
        id   # symbol -> symbol: return as-is
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
      # ENSG -> Entrez -> symbol (original logic preserved exactly)
      entrez <- mget(id, envir = ensg2eg, ifnotfound = NA)
      entrez <- sapply(entrez, function(x) {
        if (length(x) > 1L || is.na(x)) "---" else x
      })
      hugo <- mget(entrez, envir = sym, ifnotfound = NA)
      sapply(hugo, function(x) {
        if (length(x) > 1L) paste(x, collapse = " /// ") else x
      })
    } else {
      setNames(id, id)   # ensembl -> ensembl
    }

  } else if (is_entrez) {
    if (eff_output == "symbol") {
      # Entrez -> symbol (direct)
      result_sym <- mget(id, envir = sym, ifnotfound = NA)
      sapply(result_sym, function(x) {
        if (length(x) > 1L)      paste(x, collapse = " /// ")
        else if (is.na(x[[1L]])) NA_character_
        else                     x[[1L]]
      })
    } else {
      # Entrez -> ENSG (direct)
      result_ensg <- mget(id, envir = ensg, ifnotfound = NA)
      sapply(result_ensg, function(x) {
        if (length(x) > 1L)      paste(x, collapse = " /// ")
        else if (is.na(x[[1L]])) NA_character_
        else                     x[[1L]]
      })
    }

  } else {
    # Symbol input
    if (eff_output == "ensembl") {
      # symbol -> Entrez -> ENSG (original logic preserved exactly)
      entrez <- mget(id, envir = sym2eg, ifnotfound = NA)
      entrez <- sapply(entrez, function(x) {
        if (length(x) > 1L || is.na(x)) "---" else x
      })
      result_ensg <- mget(entrez, envir = ensg, ifnotfound = NA)
      sapply(result_ensg, function(x) {
        if (length(x) > 1L) paste(x, collapse = " /// ") else x
      })
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
#' convert.alias("TRPV4")
#' @export
convert.alias <-
  function(id, species = c("Human", "Mouse"), db = NULL)
  {
    if (missing(id))
      stop("Need input ID vector!")
    if (is.null(db)) {
      species <- match.arg(species)
      if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
        stop(
          "The Bioconductor package 'org.Hs.eg.db' is required for this function.\n",
          "Install it via BiocManager::install('org.Hs.eg.db').",
          call. = FALSE
        )
      }
      # Get the namespace without using ::
      orgdb <- switch(species,
                      Human=getNamespace("org.Hs.eg.db"),
                      Mouse=getNamespace("org.Mm.eg.db")
      )
      db <- switch(species,
                   Human=orgdb[["org.Hs.eg.db"]],
                   Mouse=orgdb[["org.Mm.eg.db"]]
      )
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
    biom.ids <- get.bm(values, biom.data.set, biom.mart, host, biom.filter, biom.attributes, biom.cache, use.cache, biomart.fallback, chunk.size, verbose = verbose)
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
