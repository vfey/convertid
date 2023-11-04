---
title: CRAN package convertid
---

## Resubmission 2023-11-04
This is a resubmission of a maintenance release of the package. The version number was increased to 0.1.7 after addressing a bug related to *dbplyr's* database input/output framework. The fix implements a workaround allowing the user to optionally disable the use of the Bioc file cache when getting/setting CURL SSL options which in turn requires *dbplyr*.

### Test environments (2023-11-04 - )
* local OS X install: x86_64-apple-darwin22.6.0, R 4.3.1
* win-builder (devel, release and oldrelease)
* Red Hat Enterprise Linux release 9.2 (Plow), R 4.3.1

## Resubmission 2023-10-06
This is a resubmission of a maintenance release of the package. The version number was increased to 0.1.6 after addressing a bug that was was discovered shortly after submitting the previous version. The bug was related to an inconsistent cache handling in the biomaRt package. In consequence, the fixes attempt to make testing and setting the cache location more robust to failure.

### R CMD check results
There were no ERRORs, WARNINGs or NOTEs.  

## Resubmission 2023-10-04
This is a resubmission of a maintenance release of the package. The version number was increased to 0.1.5 after adding 'https://' to all host URLs to comply with Ensembl requirements. In addition, an internal function, `.setBiomaRtCacheLocation()`, was added to address a problem with the biomaRt cache on non-standard Linux installations.  

### R CMD check results
There were no ERRORs, WARNINGs or NOTEs.  

### Test environments (2023-10-04 - 2023-10-06)
* local OS X install: x86_64-apple-darwin22.6.0, R 4.3.1
* win-builder (devel, release and oldrelease)
* CentOS Linux release 7.9.2009 (Core) [:core-4.1-amd64:core-4.1-noarch], R 4.3.1

## Resubmission 2021-09-13
This is a resubmission of a maintenance release of the package. The version number was increased to 0.1.3 after correcting a bug in function `todisp2()`.  
In addition, argument `verbose` was added to functions generating progress info messages to control their printing.  

### R CMD check
There were no ERRORs, WARNINGs or NOTEs.  

## Resubmission 2021-09-02
This is a resubmission of a new release of the package. The version number was increased to 0.1.2 after making corrections to the default behavior in function `todisp2()`. The functionality has not changed but the function is more resistant to wrong user input. In addition, a new function was added, `likely_symbol`.  

### Notes
Since both the local OS X check as well as the win.builder check returned NOTEs on example timings those examples were wrapped in `\dontrun{}`.  

There was 1 NOTE during the local OS X check:  

```
R CMD check --as-cran convertid_0.1.2.tar.gz
.
* checking examples ... NOTE
Examples with CPU (user + system) or elapsed time > 5s
               user system elapsed
likely_symbol 5.934  0.117   7.848
.
```

There were 2 NOTEs during win.builder check:  

```
** running examples for arch 'i386' ... [63s] NOTE
Examples with CPU (user + system) or elapsed time > 10s
               user system elapsed
likely_symbol 47.52   7.55   55.39
** running examples for arch 'x64' ... [53s] NOTE
Examples with CPU (user + system) or elapsed time > 10s
              user system elapsed
likely_symbol 35.3   7.53    43.3
```

## Resubmission 2021-08-20
This is a resubmission. The version was increased to 0.1.1 after addressing the comments by CRAN staff member Uwe Ligges:

* single-quoted package and software, e.g. 'AnnotationDbi'
* unquoted function names, e.g., convertId2()

The word "Ensembl" is not misspelled but refers to the Ensembl project. It is used, e.g., in the DESCRIPTION of the `biomaRt` package.

In addition to the above corrections, I added a co-author and made minor modifications to the package description.

### Notes
The first win.builder check gave the following timings when running examples:  

```
Examples with CPU (user + system) or elapsed time > 10s
            user system elapsed
get.bm     10.09   0.11   32.23
convert.bm  9.61   0.05   12.73
todisp2     9.33   0.01   13.49
```

For this reason, examples in these functions were wrapped in `\dontrun{}`.

## Test environments (2021-08-18 - 2021-09-13)
* local OS X install: x86_64-apple-darwin17.0, R 4.0.2
* win-builder (devel and release)
* CentOS Linux release 7.9.2009 (Core) [:core-4.1-amd64:core-4.1-noarch], R 4.0.4
