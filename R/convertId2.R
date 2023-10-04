#' @title Convert Gene IDs Between Each Other and Fetch Annotations from Biomart
#' @docType package
#' @name convertid
#' @description Gene Symbols or Ensembl Gene IDs are converted using the Bimap interface in 'AnnotationDbi' in convertId2() but
#'     that function is only provided as fallback mechanism for the most common use cases in data analysis. The main function
#'     in the package is convert.bm() which queries Biomart using the full capacity of the API provided through the
#'     'biomaRt' package. Presets and defaults are provided for convenience but all "marts", "filters" and "attributes"
#'     can be set by the user. Function convert.alias() converts Gene Symbols to Aliases and vice versa and function
#'     likely_symbol() attempts to determine the most likely current Gene Symbol.
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
#' @import org.Hs.eg.db
#' @import org.Mm.eg.db
#' @import AnnotationDbi
#' @import plyr
#' @import stringr
#' @import biomaRt
#' @import xml2
#' @importFrom stats na.omit
#' @importFrom utils read.delim
#' @importFrom rappdirs user_cache_dir
NULL
#' Convert Gene Symbols to Ensembl Gene IDs or vice versa
#' @description \command{convertId2()} uses the Bimap interface in AnnotationDbi to extract information from
#'     annotation packages. The function is limited to Human and Mouse annotations and is provided only as
#'     fallback mechanism for the most common use cases in data analysis. Please use the Biomart interface
#'     function \code{convert.bm()} for more flexibility.
#' @param id (\code{character}). Vector of gene symbols.
#' @param species (\code{character}). One of "Human" and "Mouse". Defaults to "Human".
#' @return A named character vector where the input IDs are the names and the query results the values.
#' @seealso \code{\link[AnnotationDbi]{Bimap-envirAPI}}
#' @examples
#' convertId2("ENSG00000111199")
#' convertId2("TRPV4")
#' @export
convertId2 <-
  function (id, species = c("Human", "Mouse"))
  {
    species <- match.arg(species)
    if (species == "Human") {
      ensg2eg.env <- org.Hs.eg.db::org.Hs.egENSEMBL2EG
      sym.env <- org.Hs.eg.db::org.Hs.egSYMBOL
      sym2eg.env <- org.Hs.eg.db::org.Hs.egSYMBOL2EG
      ensg.env <- org.Hs.eg.db::org.Hs.egENSEMBL
      ensg <- "ENSG"
    }
    if (species == "Mouse") {
      ensg2eg.env <- org.Mm.eg.db::org.Mm.egENSEMBL2EG
      sym.env <- org.Mm.eg.db::org.Mm.egSYMBOL
      sym2eg.env <- org.Mm.eg.db::org.Mm.egSYMBOL2EG
      ensg.env <- org.Mm.eg.db::org.Mm.egENSEMBL
      ensg <- "ENSMU"
    }
    if (length(id) == 1) {
      if (length(grep(ensg, id)) > 0) {
        if (AnnotationDbi::exists(id, envir = ensg2eg.env)) {
          entrez <- get(id, envir = ensg2eg.env)
          if (length(entrez) > 1) {
            sym <- NA_character_
          }
          else {
            if (AnnotationDbi::exists(entrez, envir = sym.env)) {
              sym <- paste(get(entrez, envir = sym.env),
                           collapse = " /// ")
            } else {
              sym <- NA_character_
            }
          }
        }
        else {
          sym <- NA_character_
        }
      }
      else {
        if (AnnotationDbi::exists(id, envir = sym2eg.env)) {
          entrez <- get(id, envir = sym2eg.env)
          if (length(entrez) > 1) {
            sym <- NA_character_
          }
          else {
            if (AnnotationDbi::exists(entrez, envir = ensg.env)) {
              sym <- paste(get(entrez, envir = ensg.env),
                           collapse = " /// ")
            } else {
              sym <- NA_character_
            }
          }
        }
        else {
          sym <- NA_character_
        }
      }
      names(sym) <- id
      return(sym)
    }
    else {
      if (length(grep(ensg, id[1])) > 0) {
        entrez <- mget(id, envir = ensg2eg.env, ifnotfound = NA)
        entrez <- sapply(entrez, function(x) {
          if (length(x) > 1 || is.na(x)) {
            "---"
          }
          else {
            x
          }
        })
        hugo <- mget(entrez, envir = sym.env, ifnotfound = NA)
        hugo <- sapply(hugo, function(x) {
          if (length(x) > 1) {
            paste(x, collapse = " /// ")
          }
          else {
            x
          }
        })
        names(hugo) <- id
        return(hugo)
      }
      else {
        entrez <- mget(id, envir = sym2eg.env, ifnotfound = NA)
        entrez <- sapply(entrez, function(x) {
          if (length(x) > 1 || is.na(x)) {
            "---"
          }
          else {
            x
          }
        })
        ensg <- mget(entrez, envir = ensg.env, ifnotfound = NA)
        ensg <- sapply(ensg, function(x) {
          if (length(x) > 1) {
            paste(x, collapse = " /// ")
          }
          else {
            x
          }
        })
        names(ensg) <- id
        return(ensg)
      }
    }
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
      db <- switch(species,
                   Human=org.Hs.eg.db::org.Hs.eg.db,
                   Mouse=org.Mm.eg.db::org.Mm.eg.db
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
#' @param biom.attributes \code{character} vector. Biomart attributes, i.e., type of desired result(s); make sure query id type is included!
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
  function(dat, id="ID", biom.data.set = c("human", "mouse"),
           biom.mart=c("ensembl", "mouse", "snp", "funcgen", "plants"),
           host="https://www.ensembl.org", biom.filter="ensembl_gene_id",
           biom.attributes=c("ensembl_gene_id","hgnc_symbol","description"),
           sym.col="hgnc_symbol", rm.dups=FALSE, verbose = FALSE)
  {
    if (id=="row.names") {
      values <- rownames(dat)
    } else {
      values <- dat[[id]]
    }
    biom.ids <- get.bm(values, biom.data.set, biom.mart, host, biom.filter, biom.attributes, verbose = verbose)
    gene.lab <- merge(biom.ids, dat, by.x=biom.filter, by.y=id, all.y=TRUE, all.x=FALSE, sort=TRUE)
    if (rm.dups) {
      if (verbose) message("  Removing ", length(which(duplicated(gene.lab[[biom.filter]]))), " duplicated row(s)...")
      gene.lab <- gene.lab[!duplicated(gene.lab[[biom.filter]]), ]
    }
    if (any(gene.lab[[sym.col]]=="") || any(is.na(gene.lab[[sym.col]]))) {
      if (verbose) message("  Replacing ", length(which(gene.lab[[sym.col]]=="" | is.na(gene.lab[[sym.col]]))), " missing Gene Symbols by", sQuote(biom.filter), "...")
      gene.lab[[sym.col]][gene.lab[[sym.col]]=="" | is.na(gene.lab[[sym.col]])] <- gene.lab[[biom.filter]][gene.lab[[sym.col]]=="" | is.na(gene.lab[[sym.col]])]
    }
    return(gene.lab)
  }

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
           biom.cache = Sys.getenv(x = "BIOMART_CACHE"),
           use.cache = TRUE,
           verbose = FALSE)
  {
    biom <- match.arg(biom.mart)
    if (biom=="plants" && host == "https://www.ensembl.org") {
      if (verbose) message(sQuote("Plants"), "mart requested. Setting host to ", sQuote("https://plants.ensembl.org"), "...")
      host <- "https://plants.ensembl.org"
    }
    marts <- biomaRt::listMarts(host=host)[["biomart"]]
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

    if (use.cache) {
      .setBiomaRtCacheLocation(biom.cache)
    }
    if (verbose) message("  Information requested: ", sQuote(setdiff(biom.attributes, biom.filter)), "...")
    biomaRt::getBM(attributes=biom.attributes, filters=biom.filter, values=values, mart=mart, useCache = use.cache)
  }

#'
#' Convenience Function to Convert Ensembl Gene IDs to Gene Symbols
#' @description \command{todisp2()} uses Biomart by employing \command{get.bm()} to retrieve Gene Symbols for a set of Ensembl
#'     Gene IDs. It is mainly meant as a fast way to convert IDs in standard gene expression analysis output to Symbols,
#'     e.g., for visualisation, which is why the input ID type is hard coded to ENSG IDs. If Biomart is not available
#'     the function can fall back to use \command{convertId2()} or a user-provided data frame with corresponding ENSG IDs and
#'     Symbols.
#' @param ensg (\code{character}). Vector of Ensemble Gene IDs. Other ID types are not yet supported.
#' @param lab (\code{data.frame}). A data frame with Ensembl Gene IDs as row names and Gene Symbols in the only column.
#' @param biomart (\code{logical}). Should Biomart be used? Defaults to \code{TRUE}.
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
todisp2 <- function(ensg, lab=NULL, biomart=TRUE, verbose = FALSE)
{
  if (biomart) {
    if (!length(grep("^ENS[A-Z]{0,}[0-9]{11}", ensg[1]))) {
      if (verbose) message("    Input is not Ensembl Gene IDs. Doing nothing.")
      return(ensg)
    }
    sym <- get.bm(ensg, biom.data.set="hsapiens_gene_ensembl", biom.mart="ensembl", host="https://www.ensembl.org",
                  biom.filter="ensembl_gene_id", biom.attributes=c("ensembl_gene_id","hgnc_symbol"),
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
  gene.lab <- merge(data.frame(ensembl_gene_id=ensg, stringsAsFactors=FALSE), sym, by="ensembl_gene_id", sort=FALSE)
  if (verbose) message("done")
  if (any(gene.lab$hgnc_symbol=="" | is.na(gene.lab$hgnc_symbol))) {
    if (verbose) message("    Replacing ", length(which(gene.lab$hgnc_symbol=="" | is.na(gene.lab$hgnc_symbol))), " missing Gene Symbol(s) by Ensembl IDs...")
    replace <- gene.lab$hgnc_symbol=="" | is.na(gene.lab$hgnc_symbol)
    gene.lab$hgnc_symbol[replace] <- gene.lab$ensembl_gene_id[replace]
  }
  return(gene.lab$hgnc_symbol)
}
