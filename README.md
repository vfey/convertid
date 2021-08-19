# R package `convertid`
Convert Gene IDs Between Each Other and Fetch Annotations from Biomart

Gene Symbols or Ensembl Gene IDs are converted using the Bimap interface in AnnotationDbi in 'convertId2()' but that function is only provided as fallback mechanism for the most common use cases in data analysis. The main function in the package is 'convert.bm()' which queries Biomart using the full capacity of the API. Presets and defaults are provided for convenience but all 'marts', 'filters' and 'attributes' can be set by the user. Function 'convert.alias()' converts Gene Symbols to Aliases and vice versa.
