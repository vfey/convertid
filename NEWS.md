# convertid 0.4.0

## New features

* the `output` argument of `convertId2()` gained the choice `"entrez"`, now
  fully supporting conversion to Entrez Gene IDs in addition to Gene Symbols and
  Ensembl Gene IDs. `output = "auto"` keeps the previous behaviour: Ensembl
  input returns symbols, symbol input returns Ensembl IDs, Entrez input returns
  symbols.; `"symbol"`, `"ensembl"` and `"entrez"` force the return type
  regardless of the input type. Entrez output is useful for offline,
  network-independent retrieval, e.g. for KEGG pathway enrichment.  

* `todisp2()` gained `biomart.fallback` and `chunk.size` as well, and now
  degrades gracefully: if every BioMart host fails it falls back to the data
  frame in `lab`, or to `convertId2()`, instead of aborting.  

* `get.bm()` now detects versioned Ensembl gene IDs (e.g. `ENSG00000111199.5`)
  in the input and stops with an informative message, since BioMart's
  `ensembl_gene_id` filter expects unversioned IDs.  

## Bug fixes

* `convert.alias()` checked for `org.Hs.eg.db` even when `species = "Mouse"`,
  producing a misleading error message.

* `likely_symbol()` returned different results above and below
  `index_threshold`: the indexed path recorded only the first alias of a
  matched gene, so a query naming any later alias was returned unresolved. The
  query itself is now trimmed as well, so a padded symbol still matches.

* `convertId2()` aborted on input consisting entirely of `NA`, and on empty
  input. It now returns an all-`NA` result of the same length.

* `unify_gene_ids()` returned `NA` for every `hgnc_symbol` when all BioMart
  hosts failed but the AnnotationDbi lookup succeeded. The AnnotationDbi result
  was recorded in a separate column and never promoted into `hgnc_symbol`, so
  the degradation to "AnnotationDbi results only" that the documentation and
  the verbose output describe did not actually take place.

* `unify_gene_ids()` could abort with "NAs are not allowed in subscripted
  assignments" when the gene symbol column contained `NA`.

* The deduplication filters could fabricate a row: an `NA` in a row-selection
  condition does not drop the row but yields one filled with `NA`, which could
  then be taken for the single surviving candidate of its group. One filter
  could also abort outright on an `NA` gene name.

* The Ensembl-placeholder fix in the deduplication step never ran in ENSG-only
  mode, where there is no gene symbol column.

* `convert.bm()` accepts both `"row.names"` and `"rownames"` as the special
  value of `id` selecting the row names. Only the former worked, while the
  documentation named the latter.

* `.with_biomart_fallback()` treated a `NULL` return value from a successful
  query as a failed host, and reported an unhelpful error when no host was
  configured at all.

* `get.bm()` aborted on R >= 4.2 when `biom.data.set` was given as a
  multi-element vector of user-supplied data set names.

* `likely_symbol()` aborted when a query containing `|` could not be resolved.
  Such a query is now returned unchanged, as any other unresolved symbol is.

## Documentation

* Documented the `previous_symbol` column of `likely_symbol()` and corrected the
  column list for `output = "all"`.

* Corrected the return value description of `unify_gene_ids()`: the
  intermediate lookup columns are dropped unless `keep_intermediates = TRUE`.

## Internal

* Extracted the deduplication filters, the AnnotationDbi lookup helpers and the
  HGNC field tokeniser into documented, individually testable unexported
  functions.

* Test helpers no longer live in `R/`, and BioMart availability is probed with
  an actual connection rather than a port-level check.

* `AnnotationDbi` moved from `Depends` to `Imports`. The package imports
  everything it needs from it, so nothing changes for the functions in this
  package, but `library(convertid)` no longer attaches `AnnotationDbi` to the
  search path. Code that relied on that side effect should call
  `library(AnnotationDbi)` itself.

* The test suite uses the third edition of `testthat`. `URL` and `BugReports`
  fields were added to `DESCRIPTION`.
