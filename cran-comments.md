---
title: CRAN package convertid
---

## Resubmission 2021-08-20
This is a resubmission. The version was increased to 0.1.1 after addressing the comments by CRAN staff member Uwe Ligges:

* single-quoted package and software, e.g. 'AnnotationDbi'
* unquoted function names, e.g., convertId2()

The word "Ensembl" is not misspelled but refers to the Ensembl project. It is used, e.g., in the DESCRIPTION of the `biomaRt` package.

In addition to the above corrections, I added a co-author and made minor modifications to the package description.

## Notes
The first win.builder check gave the following timings when running examples:  

```
Examples with CPU (user + system) or elapsed time > 10s
            user system elapsed
get.bm     10.09   0.11   32.23
convert.bm  9.61   0.05   12.73
todisp2     9.33   0.01   13.49
```

For this reason, examples in these functions were wrapped in `\dontrun{}`.

## Test environments
* local OS X install: x86_64-apple-darwin17.0, R 4.0.2
* win-builder (devel and release)
* CentOS Linux release 7.9.2009 (Core) [:core-4.1-amd64:core-4.1-noarch], R 4.0.4

## R CMD check results
There were no ERRORs or WARNINGs.

There was 1 NOTE:

```
R CMD check --as-cran convertid_0.1.0.tar.gz
.
* checking CRAN incoming feasibility ... NOTE     
Maintainer: ‘Vidal Fey <vidal.fey@gmail.com>’
.
```
