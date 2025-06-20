poppr 2.9.7
===========

MISC
----

* compatibility fix for ggplot2 4.0.0 (reported: @teunbrand, #263, 
  fixed: @zkamvar, #264)
* CRAN maintenance: documentaiton fixes

poppr 2.9.6
===========

MISC
----

* Poppr itself will no longer accidentally modify the random seed when it is
  attached (found by @timtaylor,
  https://github.com/grunwaldlab/poppr/issues/259)

poppr 2.9.5
===========

CRAN MAINTENANCE
----------------

* The `usenames` option from xcolor has been removed from the algorithms and
  equations vignette.
* Old documentation has been fixed to conform with CRAN standards.

poppr 2.9.4
===========

CRAN MAINTENANCE
----------------

* a function declaration was added for `SEXP omp_test()`
* failing tests were fixed.
* code where object classes were compared with `==` was fixed to use `inherits()`

DEPENDENCIES
------------

* {RClone} has been removed as a suggested package as it was archived on CRAN
  some time ago.

poppr 2.9.3
===========

CRAN MAINTENANCE
----------------

* C headers were re-arranged to place R-specific headers _after_ OMP headers to
  avoid problems in clang 13 (@zkamvar, https://github.com/grunwaldlab/poppr/pull/246).

poppr 2.9.2
===========

CRAN MAINTENANCE
----------------

* A test for `bruvo.between()` was fixed in a superficial way to avoid an error
  on CRAN R-devel (@zkamvar, https://github.com/grunwaldlab/poppr/pull/242).

poppr 2.9.1
===========

DEPENDENCY UPDATE
----------------

* The {phangorn} package is no longer imported. We only used one line of code
  from the package (upgma), so we copied over the implementation with the 
  author's permission (@KlausVigo, https://github.com/grunwaldlab/poppr/pull/237).

poppr 2.9.0
===========

NEW FEATURES
------------

* `bruvo.between()` will calculate bruvo's distances between a query dataset
  and a reference dataset (@davefol, #223)

DEPRECATION
-----------

* The argument `blacklist` has been deprecated in favor of `exclude` for the
  following functions: `bruvo.msn()`, `poppr.msn()`, `clonecorrect()`, `poppr()`, 
  `mlg.table()`, `mlg.crosspop()`, and `popsub()`. It will be removed in the
  poppr version 2.10 (@zkamvar, #218)

GITHUB
------

* The default branch for the repository is now "main" (@zkamvar, #218)

BUG FIX
-------

* `genind2genalex()` no longer converts diploid sequence data to zeros on export
  This fixes #231 (@zkamvar, #233).
* `bitwise.ia()` will no longer have integer overflows early on Windows
  (@zkamvar, #235)

poppr 2.8.7
===========

CRAN MAINTENANCE
----------------

* Links to other packages are now formatted correctly.

poppr 2.8.6
===========

BUG FIX
-------

* `read.genalex()` now can import columns with entirely "T" alleles as "T" 
  instead of "TRUE". (See https://github.com/grunwaldlab/poppr/issues/214
  for details).

NEW IMPORTS
-----------

* Progress bars are now implemented via the `{progressr}` package,
  which gives the user control over what should be displayed (or not). For
  example, to get auditory updates instead of a progress bar, you can use the
  "beepr" package and set `progressr::handlers("beepr")`. This will play an alert
  for each step (~50) and a final sound. To suppress all progress bars entirely,
  you can use `progressr::handlers("void")`. These progress bars have replaced
  the `utils::txtProgressBar()` and `dplyr::progress_estimated()` bars.

poppr 2.8.5
===========

CRAN MAINTENANCE
----------------

* The output from the `poppr()` function will no longer contain factor columns
  for the population names or the file names. This is due to R 4.0.0 changing
  behavior with the `stringsAsFactors` default in `data.frame()`. (See 
  https://github.com/grunwaldlab/poppr/issue/212 for details).

poppr 2.8.4
===========

CRAN MAINTENANCE
----------------

* Update documentation for S4 method

poppr 2.8.3
===========

BUG FIX
-------

* `read.genalex()` now correctly parses strata when the user imports data that
  contains duplicated data AND has some individuals named as integers less than
  the number of samples in the data (prepended by zeroes) 
  (See https://github.com/grunwaldlab/poppr/pull/202).

NEW FEATURES
------------

* MSN functions: nodes with single populations displayed as circles instead of 
  pies. (@fdchevalier, https://github.com/grunwaldlab/poppr/pull/203)

MISC
----

* `mlg.vector()` is now safer as it now uses a for loop instead of a
  function with the out-of-scope operator (`<<-`) (see #205)
* `shufflepop()` is now safer as it now uses a for loop instead of a
  function with the out-of-scope operator (`<<-`) (see #205)
* The MLG class gains a new `distenv` slot, which will store the environment
  where the distance function or matrix exists. This is accompanied by an 
  accessor of the same name (see #206).
* `"mlg.filter<-"()` replacement methods will no longer search the global
  environment when evaluating the distance function or matrix (see #206).
* Tests for `mlg.filter()` no longer assign objects to the global environment
* DOIs for the publications have been added to the DESCRIPTION

poppr 2.8.2
===========

MISC
----

* Tests relying on randomization were updated before R 3.6.0, which fixes a
  biased randomization. This has no visible impact on users.
  See https://github.com/grunwaldlab/poppr/pull/198 for details.

poppr 2.8.1
===========

BUG FIX
-------

* An error that appeared in some AMOVA calls with genind objects with character-
  based alleles was fixed 
  (see https://github.com/grunwaldlab/poppr/issues/190 for details)

DOCUMENTATION
-------------

* `aboot()` documentation was updated to add the citation and make clear its
  purpose and limitations.

MISC
----

* DOIs have been hyperlinked to doi.org instead of dx.doi.org (#188, @katrinleinweber)

poppr 2.8.0
===========

BUG FIX
-------

* `win.ia()` now has more consistent behavior with chromosome structure and will
  no longer result in an integer overflow. 
  (see https://github.com/grunwaldlab/poppr/issues/179). Thanks to @MarisaMiller
  for the detailed bug report.
* `plot_filter_stats()` will plot stats if supplied a list of thresholds.

ALGORITHMIC CHANGE
------------------

* `win.ia()` may result in slightly different results because of two changes:
  1. The windows will now always start at position one on any given chromosome.
  This will result in some windows at the beginning of chromosomes having a 
  value of `NA` if the first variant starts beyond the first window. 
  2. Windows are now calculated for each chromosome independently. The previous
  version first concatenated chromosomes with at least a window-sized gap 
  between the chromosomes, but failed to ensure that the window always started
  at the beginning of the chromosome. This version fixes that issue. 
  (see https://github.com/grunwaldlab/poppr/issues/179).
  
DEPRECATION
-----------

* The `chromosome_buffer` argument for `win.ia()` has been permanently set to
  `TRUE` and deprecated as it is no longer used. 

NEW FEATURES
------------

* `poppr.amova()` will now handle genlight/snpclone objects.
   See https://github.com/grunwaldlab/poppr/pull/185 for details.

* `bitwise.dist()` now has two new options: `euclidean` and `scale_missing`.
  When both of these are set to `TRUE`, the distance measured will be Euclidean
  scaled for the amount of missing data in each comparison. This matches the 
  output of base R's `dist()` function at a fraction of time and memory. 
  See https://github.com/grunwaldlab/poppr/pull/176 for details.

* `make_haplotypes()` is now a generic defined for both genind and genlight.
  
* `genind2genalex()` will no longer write to "genalex.csv" by default. Instead,
  it will warn the user and write to a temporary file. 
  See https://github.com/grunwaldlab/poppr/issues/175 for details.

* `genind2genalex()` now has an `overwrite` parameter set to `FALSE` to prevent
  accidental overwriting of files. 

* `win.ia()` has a new argument `name_window`, which will give each element in 
  the result the designation of the terminal position of that window. Thanks to
  @MarisaMiller for the suggestion!

* `pair.ia()` can now calculate p-values via permutations.
   (See https://github.com/grunwaldlab/poppr/issues/180 for details)

DOCUMENTATION
-------------

* `cutoff_predictor()` was added to the MLG vignette

MISC
----

* The RClone package has been added to suggests
  (see https://github.com/grunwaldlab/poppr/issues/187)

poppr 2.7.1
============

MISC
----

* Missing documentation for `poppr.amova()` has been added.
* Polysat is now listed in imports.

poppr 2.7.0
============

NEW FUNCTIONS
-------------

* `make_haplotypes()` will split your data into pseudo-haplotypes for use in
  AMOVA-like analyses. This was a previously internal function, but has been 
  promoted to a user-facing function in this version. 

* `as.genambig()` will convert genind/genclone objects to Polysat's "genambig"
  class. Note that polysat must be installed for this to work. 

ALGORITHMIC CHANGE
------------

* AMOVA will now default to using euclidean distance. This affects all 
  calculations made with `within = FALSE` or `filter = TRUE` without a 
  user-supplied distance. This will not have affect those with haploid or 
  diploid data using `within = TRUE`. The dissimilarity distance is equivalent 
  to a squared euclidean distance for haploid genotypes, but not for any higher 
  ploidy. Those using `filter = TRUE` without specifying a distance should use
  a euclidean threshold. This should not be an issue for those who simply want
  to group isolates with missing data, however as a zero distance is the same
  for euclidean and dissimilarity. Thanks goes to Patrick Meirmans for alerting
  me to this error.


NEW FEATURES
------------

* AMOVA will now calculate within-individual variance for polyploid data.

MISC
----

* printing of AMOVA will now better handle any changes in methods from pegas or
  ade4. 

poppr 2.6.1
===========

BUG FIX
-------

* An out-of-bounds memory access error in `bitwise.dist()` was fixed.
  See https://github.com/grunwaldlab/poppr/issues/169 for details.

poppr 2.6.0
===========

NEW FUNCTIONS
-------------

* The new function `boot.ia()` is conceptually similar to `resample.ia()`,
  except it resamples with replacement. 
  
NEW FEATURES
------------

* The function `resample.ia()` now can resample individuals weighted by their
  Psex value. 
* The minimum spanning networks will now scale nodes by area instead of radius.
  This gives a more accurate picture of the differences between MLGs. See
  https://github.com/grunwaldlab/poppr/issues/154 for details.
* A legend for samples/node is now added to all minimum spanning networks. See
  https://github.com/grunwaldlab/poppr/issues/158 for details.
* The imsn() option for node size scale has been changed to a slider.
  
BUG FIX
-------

* An issue where data with sample names containing apostrophes could not be
  imported was fixed (Identified in 
  https://github.com/grunwaldlab/poppr/issues/156).
* a bug in `imsn()` where custom MLGs would result in an error was fixed. See
  https://github.com/grunwaldlab/poppr/issues/155 for details.
* a bug in `plot_poppr_msn()` where setting `scale.leg = FALSE` would result in a
  very small MSN plot was fixed.
* `mlg()` now works properly for snpclone and genlight objects. See 
  https://github.com/grunwaldlab/poppr/issues/155 for details.

DEPENDENCIES
------------

* The minimum version of igraph has been set to 1.0.0.

MISC
----

* The MSN is now plotted last in `plot_poppr_msn()` so additional legends can
  be added if necessary.

poppr 2.5.0
===========

ALGORITHMIC CHANGE
------------------

* Identified in https://github.com/grunwaldlab/poppr/issues/139, Bruvo's 
  distance will now consider all possible combinations of ordered alleles in the
  calculation under the genome addition and loss models for missing data. This
  will affect those who have polyploid data that contain more than one missing 
  allele at any genotype
  
  To facilitate comparison, the global option old.bruvo.model, has been created.
  By default it is set to FALSE, indicating that poppr should use the ordered 
  allele combinations. If the user wants to use the method considering unorderd 
  allele combinations, they can set options(old.bruvo.model = TRUE)
  
  It must be repeated that this does not affect haploid or diploid comparisons, 
  those that use the infinite alleles model, or those who do not have more than 
  one missing allele at any genotype.

DEPRECATION
-----------

* The warning for a short repeat length vector for Bruvo's distance is 
  deprecated and will become an error in the future
* jack.ia is deprecated in favor of resample.ia for clarity.
  
BUG FIX
-------

 * A bug in read.genalex() where removed samples would have incorrect strata
   labels was fixed. Thanks to Hernán Dario Capador-Barreto for identifying it.
   See https://github.com/grunwaldlab/poppr/issues/147.
  
MISC
----

* The internal plotting function for mlg.table now uses tidy evaluation for
  dplyr versions > 0.5.0
* The package reshape2 was removed from imports and replaced with base functions
  (see https://github.com/grunwaldlab/poppr/issues/144 for details)

NEW IMPORTS
-----------

* Due to the migration to dplyr version 0.7.0, poppr now imports the "!!"
  operator from the rlang package

poppr 2.4.1
===========

BUG FIX
-------

* A corner case where repeat length vectors out of order would be erroneously
  subset with `test_replen()` and `fix_replen()` has been fixed. See
  https://github.com/grunwaldlab/poppr/issues/136 for details.
* All functions that perform filtering will now run serially due to a bug on
  Fedora machines with at least two threads. Details can be found at
  https://github.com/grunwaldlab/poppr/issues/138.


poppr 2.4.0
===========

NEW FUNCTIONS
-------------

* `jack.ia()` will randomly jackknife your sample to a specified n (default
  is the number of MLG), and calculate the index of association over 
  multiple iterations, giving a distribution of possible values at a given
  sample size.

NEW FEATURES
------------

* The function `mlg.table()` gains new parameters, "color" and "background". The
  "color" parameter will create a single barplot with colors representing 
  populations while the "background" parameter will create a background plot 
  showing the abundance of MLGs across populations within the facets.
* The function `win.ia()` will now take into consideration chromosomal
  coordinates when constructing windows. It has additionally acquired a new
  parameter `chromosome_buffer`, which allows the user to specify whether or not
  the window should be limited to within chromosomes.
  
NOTABLE CHANGES
---------------

* calculation of MLGs for snpclone and genlight objects will be performed via 
  distance-based methods by default. This is in contrast to the previous behavior
  where individuals were assumed to have unique genotypes. 
  see https://github.com/grunwaldlab/poppr/issues/125 for details.
* An error will be thrown when attempting to use `mlg.crosspop()` with an object
  that has < 2 populations.
* `genotype_curve()` will now remove monomorphic loci before calculation by
  default as these loci misleadingly influence the shape of the curve. This will
  change the shape of the curve if you have monomorphic loci. This change IS
  optional via the drop and dropna parameters, but it is not recommended to 
  change these parameters. 
* The calculation for `psex()` has changed to be more accurate when using
  method = "multiple". It also gains the ability to use several values of G,
  one for each population. Documentation for `psex()` has also been improved.
  For details of the change, see https://github.com/grunwaldlab/poppr/issues/101

BUG FIX
-------

* The error given when a genlight object is passed to `poppr()` now correctly
  identifies the substitute function as `diversity_stats()` and not diversity
  table (see https://github.com/grunwaldlab/poppr/issues/123). 
* GenAlEx files imported with duplicate loci will generate a warning telling the
  user that the duplicated loci have been renamed (usually to <locus name>_1)
  (see https://github.com/grunwaldlab/poppr/issues/122).
* Haplodiploids imported in genalex files will be properly treated 
  (see https://github.com/grunwaldlab/poppr/issues/124).
* `read.genalex()` will now implicitly check for the correct number of 
  individuals in the data (see https://github.com/grunwaldlab/poppr/issues/128).
* The function `poppr()` no longer throws an error if the sample > 0 and the 
  data has no population (see https://github.com/grunwaldlab/poppr/issues/130).
* A bug where round-robin allele frequencies calculated with by_pop = TRUE were
  inaccurate for all but the first population was fixed. For details, see
  https://github.com/grunwaldlab/poppr/issues/132.
* A potential integer overflow was fixed in `SEXP association_index_haploid`.
  This was a ghost from https://github.com/grunwaldlab/poppr/issues/100.
* PROTECT statements were placed around allocation statements. For details, see
  https://github.com/grunwaldlab/poppr/issues/133.
  
MISC
----

* The documentation for `bitwise.dist()` clarifies the role of the
  `differences_only` flag (see https://github.com/grunwaldlab/poppr/issues/119).
* Interruptions in C code is now handled gracefully via `R_CheckUserInterrupt()`.
  The benefit is that long-running calculations are interrupted near instantly, 
  but at the cost of a few more milliseconds of computation time.
  (see https://github.com/grunwaldlab/poppr/issues/86)
* Bruvo's distance now has complete tests for recursion as of commit 4e4fa40d16

poppr 2.3.0
===========

NEW FUNCTIONS
------------

* The function `bootgen2genind()` will help users take advantage of 
  bootstrapping distance functions from other packages that require genind 
  objects. For details, see https://github.com/grunwaldlab/poppr/issues/112 and
  https://github.com/grunwaldlab/poppr/issues/111

NEW FEATURES
------------

* There is now a `plot` parameter for the genotype curve to enable or suppress
  plotting.
* Progress bars are now automatically suppressed when running non-interactively.
  to turn them on when running non-interactively, use 
  `options(poppr.debug = TRUE)`.
* The progress bar for `ia()` and `poppr()` will now show estimated time. This
  is from dplyr's `progress_estimated()`.

MISC
-----------

* The `hist` argument in the `ia()` is deprecated in favor of `plot`.
* The x axis for the `genotype_curve()` plot is now numeric, allowing you to
  fit a smoothing function over the points without having to use the hack
  `geom_smooth(aes(group = 1))`. This is thanks to Kara Woo for pointing this
  out on twitter (https://twitter.com/kara_woo/status/783336540407685120).
* The "show" method for genclone objects now delimits populations and strata
  by a comma, avoiding confusion with multi-word population names. Thanks to
  @knausb for the fix in https://github.com/grunwaldlab/poppr/pull/116.  
* Documentation for `poppr.amova` now contains a note about significance testing
  with the ade4 function `randtest.amova`.

BUG FIX
-----------

* The subsetting methods will now properly handle mlgs when using sample names
  to subset genclone and snpclone objects. See
  https://github.com/grunwaldlab/poppr/issues/114 for details.
* A plotting bug for `mlg.table()` was fixed so that the plots now show the
  maximum value.
* Bugs with subsetting bootgen and bruvomat objects with no loci specified were 
  fixed. See https://github.com/grunwaldlab/poppr/issues/118 for details.

poppr 2.2.1
===========

BUG FIX
-----------

* A problem exporting haploid sequence data to genalex format was fixed 
  (see https://github.com/grunwaldlab/poppr/issues/108)
* Data without population structure no longer throws an error in `imsn()`

NEW FEATURES
-----------

* You can now specify layout parameters in `imsn()`
* The *.msn functions will now take into account whether or not your data is set
  to collapsed MLGs. (see https://github.com/grunwaldlab/poppr/issues/107)

MISC
----

* Errors from mlg.filter now make more sense. 
  (see https://github.com/grunwaldlab/poppr/issues/109)
* The vignette poppr_manual has been converted to HTML format.
  (see https://github.com/grunwaldlab/poppr/issues/113)

poppr 2.2.0
===========

NEW FUNCTIONS
-----------

* `incomp()` will check your data to see if there are any incomparable samples.

NEW FEATURES
-----------

* Threshold argument added to `filter_stats()` 
  (see https://github.com/grunwaldlab/poppr/issues/94) 
* The pipe operator (`%>%`) is now exported from magrittr to make chaining 
  commands easier.
* `mlg.filter()` can now return multiple statistics.
* The user can now control the size of the labels in the index of association 
  plots with the labsize and linesize arguments in the plot method.
* `private_alleles()` gains a drop argument.
* `recode_polyploids()` can now take haplodiploid data.
  
  
BUG FIX
-----------

* An issue that caused errors in the `imsn()` code output was fixed
  (see https://github.com/grunwaldlab/poppr/issues/93)
* Caught bug where the `mlg.filter()` assignment method was using `nei.dist()`
  instead of `diss.dist()` when no distance was specified.
* A bug that resulted in significantly negative values from `bitwise.ia()` with
  large sample sizes was fixed. Spotted by @knausb
  (see https://github.com/grunwaldlab/poppr/issues/100)
* Fixed issue where `plot_filter_stats()` wasn't displaying the full range of 
  MLGs
* Color vectors are now correctly parsed when passed to msn functions 
  (see https://github.com/grunwaldlab/poppr/issues/55) for details. Long color
  vectors are now accepted, albiet with warning.
* An issue where `mll.reset()` did not reset non-MLG class objects in the mlg
  slot was fixed. 
  
MISC
-----------

* Documentation for `mlg.filter()` was clarified and updated with more examples.
* The vignette "Migration from poppr version 1" has been removed.
* Previously deprecated `*hierarchy()` functions have been removed.

poppr 2.1.1
===========

NEW FEATURES
-----------

* `imsn()` now has collapsible side panels
* `nmll()` and `mll()` will now handle genind and genlight objects
* `rraf()` now gives options for minor allele correction encompassed in the 
  internal function `rare_allele_correction()`. This extends also to `pgen()` 
  and `psex()`, which must correct minor allele frequencies by default. See 
    https://github.com/grunwaldlab/poppr/issues/81 for details.

MISC
-----------

* `mlg.filter()` now defaults to using `diss.dist()`
* default threshold for `filter_stats()` is now 1e+6
* `mlg.filter()` now returns a list instead of a pairlist
* the "hist" argument for `poppr()` is now deprecated in favor of "plot"
* documentation improvements
* the show method for snpclone objects now looks distinct from genlight
* genlight objects no longer get passed to `missingno()` in `filter_stats()`
* `filter_stats()` now returns invisibly when plot = TRUE; see 
  https://github.com/grunwaldlab/poppr/issues/87 for details.

BUG FIX
------------

* When supplied an object with no strata, `clonecorrect()` will default to
  `strata = NA`.
* `read.genalex()` will no longer fail if missing data is not coded as zero;
  see https://github.com/grunwaldlab/poppr/issues/84 for details
* `missingno()` no longer removes genotypes AT specified threshold;
  see https://github.com/grunwaldlab/poppr/issues/90 for details

poppr 2.1.0
===========

IMPROVEMENTS
-----------

* `win.ia` and `samp.ia` gain a significant speedup thanks to Jonah Brooks
  implementing the code in C.
* The internal code for the `genotype_curve` has been implemented in C for a 10x
  increase in speed.

NEW FEATURES
-----------

* `poppr.msn`, `bruvo.msn`, and `plot_poppr_msn` gain the ability to take
  character vectors for color palettes. See issue #55
  (https://github.com/grunwaldlab/poppr/issues/55) for details.
* `plot_poppr_msn` returns the modified graph.
* All functions related to Bruvo's distance can now take a named vector of
  repeat lengths in any order. See issue #61
  (https://github.com/grunwaldlab/poppr/issues/61) for details.
* `aboot` gains the argument strata so that you can automatically convert genind
  to genpop.
* `genotype_curve` can now take in loci objects from pegas.
* You can now specify the maximum number of loci to analyze in `genotype_curve`.
* `filter_stats` can now optionally plot a histogram in the background.
* `bruvo.dist` can now optionally return distance matrices by locus. This is
  addresed in issue #60 (https://github.com/grunwaldlab/poppr/issues/60)
* `aboot` can now handle matrices as previously specified in the documentation.
* `aboot` can now take custom functions to calculate distance for genlight 
  objects.
* `poppr.amova` can now perform amova using the pegas implementation.

NEW FUNCTIONS
-----------

* `rrmlg` will calculate round-robin multilocus genotypes for each locus.
* `rraf` will calculate round-robin allele frequencies for each locus.
* `pgen` will calculate the probabilities of observed genotypes.
* `psex` will calculate the probability that an observed genotype will be
  observed more than once by chance.
  
NEW TRADITIONALISTS
-----------

* because we're through being cool.

MISC
-----------

* Documentation for genclone and snpclone classes are more coherent.
* Accessors added for internal MLG objects (for developers).
* The genotype accumulation curve displays the iteration and locus number
  instead of a progress bar.
* The genotype accumulation plot is now scaled from 0 to the number of observed
  mlgs.
* Documentation for `poppr.amova` no longer references the "hierarchy" slot.

BUG FIX
-----------

* Single locus data sets can now be read in with `read.genalex`. This was
  brought up in issue [#58](https://github.com/grunwaldlab/poppr/issues/58)
  (Thanks to Nick Wong for spotting it).
* A bug in `informloci` where the MAF argument wasn't being applied to P/A data
  has been fixed.
* Printing of genclone objects with mixed ploidies previously reported
  erroneously due to a sorting error. This has been fixed.
* Edge case where a missing cell in the genalex matrix was interpreted as a
  literal "NA" was fixed.
* msn functions now return nodes that are correctly named. See Issue #66
  (https://github.com/grunwaldlab/poppr/issues/66) for details.

poppr 2.0.2
===========

BUG FIX
-----------

* Definition of Hexp was fixed. It originally was mis-calculated, inflating the
  metric. It is now correctly calculated and documented. More information at
  issue [#47](https://github.com/grunwaldlab/poppr/issues/47)

poppr 2.0.1
===========

BUG FIX
-----------

* Memory leak with Bruvo's distance was fixed by @JonahBrooks in 90facb4 (issue
  #40)
* Cutoff field now works for distances other than dissimilarity in `imsn` (issue
  #41)
* Switching between data sets no longer shows an error in `imsn`
* `read.genealex` can now correctly import missing data for diploids (issue #42)

NEW FEATURES
-----------

* Startup message now tells you if poppr was compiled with OMP support.

poppr 2.0.0
===========

COMPATIBILITY UPDATE
-----------

* poppr has moved to version 2.0 due to adegenet's recent update. The hierarchy
  slot introduced in version 1.1 is now being moved to adegenet and renamed
  `strata`. For maximum backwards compatibility, all of the hierarchy methods
  still exist, but they are deprecated and will print a warning with the proper
  function to use. If you were accessing the hierarchy slot without using the
  `*hierarchy()` methods, your code will fail as the hierarchy slot now should
  only contain a formula object.

NEW DEPENDENCIES
-----------

* poppr now imports elements from dplyr and shiny
* As required as of 2015-06-29, poppr now explicitly imports: stats, graphics,
  grDevices, and utils

NEW SUGGESTS
-----------

* poppr now suggests the cowplot, poweRlaw, and polysat packages.

NEW DATA
-----------

* A data set called `Pram` containing SSR genotypes from the Sudden Oak Death
  pathogen _Phytophthora ramorum_ (Kamvar et.al., 2015)

NEW MINT FLAVOR
-----------

* refreshing!

NEW FEATURES
-----------

* The default plot for the index of association will now be a single histogram.
  The user has the option to visualize the standardized index of association
  (`index = "rbarD"`, default) or the classic index of association (`index =
  "Ia"`). If the user uses the function `ia` with the argument `valuereturn =
  TRUE`, then the resulting object can be plotted with the plot function.
* The function `poppr` will now plot all populations in a single faceted plot
  instead of one plot per population.
* `aboot` and `bruvo.boot` will now be able to utilize any function to generate
  trees (suggested in issue #18).
* The mlg slot in the genclone object can now optionally hold an MLG class
  object. This object will contain different definitions of multilocus
  genotypes, allowing the user to switch between observed, custom, and mlgs
  defined given a genetic distance threshold.
* minimum spanning network functions gain the ability to include reticulations
  using the option `include.ties`
* minimum spanning network functions gain the ability to collapse multilocus
  genotypes by genetic distance with the option `threshold`.
* `poppr.amova` gains the ability to filter multilocus genotypes before
  calculation.
* `informloci` gains the argument "MAF", which allows the specification of a
  minor allele frequency cutoff in addition to the cutoff argument. Examples
  have been updated.
* `aboot` can now take genlight objects.
* `poppr.msn` and `plot_poppr_msn` can now take genlight objects.
* `plot_poppr_msn` now gives users the option to exclude the legends.
* Default plotting for `mlg.table` will no longer produce one plot per
  population. It will now produce a single ggplot object for all populations.
  Note that the bars are no longer colored by count.
* `poppr` will no longer calculate "Hexp". Instead, Simpson's index will be
  calumniated, but the old index can be retrieved by using (N/(N - 1))*lambda.
* `poppr` can now take any statistic that can be calculated from a table of
  multilocus genotype counts.
* `genind2genalex` gains the ability to selectively write different strata.

NEW FUNCTIONS
-----------

* `mlg.filter` will contract multilocus genotypes given a genetic distance and
  threshold using one of three algorithms. It can report statistics such as the
  multilocus genotypes returned, the number of samples within each multilocus
  genotype, the thresholds at which multilocus genotypes were collapsed, and the
  genetic distance matrix that represents the new multilocus genotypes.
* `filter_stats` will show you graphical output of all the algorithms in
  `mlg.filter`.
* `cutoff_predictor` will predict the cutoff threshold from `mlg.filter`.
* `bitwise.dist` can efficiently calculate absolute genetic distance for
  genlight objects.
* `mll` "multilocus lineages" is a new replacement for `mlg.vector` which gains
  the functionality of selecting the multilocus genotype definition from the mlg
  slot.
* `nmll` counts the number of multilocus lineages
* `mll.custom` allows the user to define custom multilocus genotypes.
* `mll.levels` allows the user to edit the names of custom multilocus genotypes.
* `poppr_has_parallel` will return `TRUE` if poppr was built with OpenMP
  parallel library.
* `win.ia` calculates windows of \bar{r}_d along genlight chromosomes.
* `samp.ia` calculates \bar{r}_d for genlight object by randomly sampling a
  user-defined number of SNPs.
* `test_replen` will test repeat lengths of microsatellite markers for
  consistency.
* `fix_replen` will fix inconsistent repeat lengths for microsatellite markers.
* `diversity_stats` returns a matrix containing diversity statistics. Defaults
  to 4 found in `poppr`, but can be extended to any statistic that can be
  calculated on a vector of MLG counts.
* `diversity_boot` will bootstrap a MLG matrix over the statistics specified for
  `get_stats`. Can also perform rarefaction bootstrap.
* `diversity_ci` will calculate and plot confidence intervals for bootstrap
  resampling of an MLG matrix. This includes rarefaction to the smallest sample
  size.
* `imsn` provides an interactive shiny interface for construction of minimum
  spanning networks.
* `pair.ia` will calculate the index of association for pairs of loci and plot
  heatmaps.

NEW CLASSES
-----------

* snpclone is an extension of the genlight object that acts very much like
  genclone. It contains an mlg slot.
* MLG is an internal class that lives inside the mlg slot of snpclone and
  genclone objects. It allows the user to easily switch between multilocus
  genotype definitions.

BUG FIXES
-----------

* I was going too fast to count.

poppr 1.1.5
===========

BUG FIX
-----------

* Fixed internal bug for `fix_negative_branch` when only one branch had a
  negative edge.
* Fixed bug in `diss.dist` where a single locus would return an error.
* Fixed bug in `poppr.amova` where a single locus would return an error due to
  `repool_haplotypes`.
* Fixed bug from the future! mlg.table will now return a matrix all the time
  (Fix #25).

poppr 1.1.4
===========

BUG FIX
-----------

* Fixed an internal bug that fails only on Windows OS.

poppr 1.1.3
===========

NEW FEATURES
-----------

* new arguments to `plot_poppr_msn` to allow for easier manipulation of node
  sizes and of labeling
* read.genalex can now take read text connections as input. Addresses issue #8
* users can now specify cutoff for missing values in `aboot`

BUG FIX
-----------

* Fixed issue where monomorphic loci would cause an error in `recode_polyploids`
* Fixed logical error that would cause The infinite alleles model of Bruvo's
  distance to inflate the distance. (Found by Michael Metzger. Addresses issue
  #5).
* AMOVA can now take subset genclone objects (Addresses issue #7).
* in `mlg.table`, the mlgsub argument will now subset by name instead of index
  (fixed in #7).
* Fixed issue for neighbor-joining trees where the internal function to fix
  negative branch lengths was accidentally shuffling the corrected branches.
  Addresses issue #11.
* `diss.dist` can now be used with `aboot`

MISC
-----------

* `info_table` will print a discrete scale as opposed to colorbar when type =
  "ploidy"
* attempted to make model choices for Bruvo's distance more clear in the
  documentation

poppr 1.1.2
===========

BUG FIX
-----------

* Fixed memory allocation bug (Further addresses issue #2).
* Memory allocated in C function bruvo_dist is now properly freed.

poppr 1.1.1
===========

BUG FIX
-----------

* Fixed bug where the loss and add options for Bruvo's distance were switched.
* Fixed illegal memory access error by UBSAN. Made memory management of internal
  C functions more sane. (Addresses issue #2).
* Fixed directional quotes and em-dashes produced by Mavericks (Addresses issue
  #3).

poppr 1.1.0
===========

NEW FEATURES
-----------

* Polyploids with ambiguous genotypes are now supported in poppr. See
  documentation for `recode_polyploids` for details.
* Calculations of Bruvo's distance now features correction for partial missing
  data utilizing genome addition and genome loss models as presented in Bruvo et
  al. 2004.
* `diss.dist` now has options to return raw distances and a matrix instead of a
  dist object.
* `read.genalex` now has the option to import as a genclone object. This is the
  default action.
* `poppr.all` will be able to analyze lists of genind or genclone objects.
* `ia` now has the argument valuereturn which will return the sampled data.
* `[bruvo,poppr].msn` functions now give the user the choice to show the graph.
* `bruvo.boot` has a cleaner plot style.

NEW DATA CLASSES
-----------

* The `genclone` object is a new extension of the `genind` object from adegenet.
  This object contains slots containing population hierarchies and multilocus
  genotype definitions and will work with all analyses in adegenet and poppr.

NEW FUNCTIONS
-----------

* [get,set,name,split,add]hierarchy - functions that will manipulate the
  hierarchy slot in a `genclone` object utilizing hierarchical formulae as
  arguments for simplification.
* `setpop` will set the population of a `genclone` object utilizing model
  formulae regarding the hierarchy slot.
* `as.genclone` will automatically convert genind objects to genclone objects.
* `is.genclone` checks the validity of genclone objects.
* `poppr.amova` will run amova on any hierarchical level. This also includes the
  feature to run amova on clone censored data sets. It utilizes the ade4 version
  of amova.
* `info_table` will calculate missing data per population per locus or ploidy
  per individual per locus and gives the user the option to visualize this as a
  heatmap.
* `locus_table` will calculate diversity and evenness statistics over all loci
  in a genind or genclone object.
* `*.dist` functions will calculate Nei's distance, Rogers' Distance, Edwards'
  Distance, Reynolds' Distance, and Provestis' Distance.
* `aboot` will allow the user to create bootstrapped dendrograms for ANY
  distance that can be calculated on genind or genpop objects.
* `plot_poppr_msn` will plot minimum spanning networks produced with poppr.
* `private_alleles` will give information about the presence of private alleles
  within a genind or genclone object.
* `recode_polyploids` will take in a polyploid genind/genclone object (with
  missing alleles coded as extra zero-value allele) and recode them to have
  frequencies relative to the observed number of alleles.
* `genotype_curve` will create a genotype accumulation curve for increasing
  number of loci.
* `mlg.id` will return a list indicating the samples belonging to a specific
  multilocus genotype.

NEW DATA SETS
-----------

* Pinf - a data set of 86 isolates from different populations of the late blight
  pathogen, Phytophthora infestans. Provided by Erica Goss
* monpop - a large data set of 694 Monilinia fructicola isolates from a single
  orchard over three years. Provided by Sydney E. Everhart

NEW CAR
-----------

* Not really.

NAMESPACE CHANGES
-----------

* poppr no longer depends on pegas.
* ade4 and reshape2 are now explicitly required.

IMPROVEMENTS
-----------

* default shuffling algorithm has been implemented in C to increase speed.
* output of the mlg functions are now represented as integers to decrease their
  size in memory.
* `mlg.matrix` is now calculated faster utilizing R's internal tabulating
  capabilities.
* The function `poppr` will no longer return rounded results, but rather is
  printed with three significant digits.

MISC
-----------

* Added unit tests.
* The poppr user manual has been shortened to only include instructions on data
  manipulation.
* A new vignette, "Algorithms and Equations" gives algorithmic details for
  calculations performed in poppr.

poppr 1.0.7
===========

UPDATE
-----------

* Updated README to include link to poppr google group.

BUG FIX UPDATE
-----------

* Made last bug fix more stable (corrected on ape side).

poppr 1.0.6
===========

BUG FIX
-----------

* Fixed bug for users who have downloaded ape version 3.1 or higher where
  bruvo.boot would throw an error.

MISC
-----------

* Updated citation information.

poppr 1.0.5
===========

NOTABLE CHANGE
-----------

* The default shuffling algorithm for calculating the index of association has
  changed from multilocus-style sampling to permutation of alleles. All of the 4
  methods are available, but new assignments are as follows: Method 1: permute
  alleles, Method 2: parametric bootstrap, Method 3: non-parametric bootstrap,
  Method 4: Multilocus-style sampling. Previously, Multilocus was 1 and the rest
  followed in the same order. There should be no compatibility issues with this
  change. Functions affected: `ia`, `poppr` `shufflepop`

BUG FIX
-----------

* Bootstrapping algorithm for `bruvo.boot` function was not shuffling the repeat
  lengths for each locus resulting in potentially erroneous bootstrap support
  values. This has been fixed by implementing an internal S4 class that will
  allow direct bootstrapping of the data and repeat lengths together.
* An occasional error, "INTEGER() can only be applied to a 'integer', not a
  'NULL'" in `bruvo.boot` or `bruvo.dist` fixed.

IMPROVEMENTS
-----------

* Changes to `bruvo.boot` allow for ever so slightly faster bootstrapping.

MISC
-----------

* Permutations for I_A and \bar{r}_d are now visualized as a progress bar as
  opposed to dots.

poppr 1.0.4
===========

BUG FIX
-----------

* A previous error where bootstrap values greater than 100 were reported from
  `bruvo.boot` on UPGMA trees has been fixed.
* Fixed correction of negative branch lengths using Kuhner and Felsenstein
  (1994) normalization for NJ trees.

MISC
-----------

* github repository for poppr has changed from github.com/poppr/poppr to
  github.com/grunwaldlab/poppr

poppr 1.0.3
===========

IMPROVEMENTS
-----------

* Optimized internal sampling function to run up to 2x faster.
* Utilized rmultinom function to increase speed of bootstrap sampling methods
  for shufflepop and ia.

NEW FEATURES
-----------

* Function `informloci` will remove phylogenetically uninformative loci.

NAMESPACE
-----------

* Now importing specific functions from igraph and ape due to dependency issues.
* Removed igraph, ape, ggplot2, and phangorn form "Dependencies", but keeping
  them in "Imports".

BUG FIXES
-----------

* `read.genalex` will no longer insert an "X" in front of loci with numeric
  names.

poppr 1.0.2
===========

BUG FIXES
-----------

* Fixed bug in diss.dist function that would return an inflated distance for
  haploids.

DOCUMENTATION
-----------

* Added explanation for the index of association in poppr_manual.
* Expanded installation section to include installation instructions from
  github.

MISC
-----------

* internal permutation algorithm no longer lists permutations in reverse order

poppr 1.0.1
===========

IMPROVEMENTS
-----------

* Algorithm for the index of association was updated to increase speed.

BUG FIXES
-----------

* Removed unnecessary rounding factor for missing data in `bruvo.dist`.
* Corrected handling of duplicate entries for `read.genind`.
* Input values that are not multiples of the specified repeat length for Bruvo's
  distance are now rounded (as opposed to being forced as integers).

MISC
-----------

* Vignette updated for aesthetics and to reflect algorithmic changes.

poppr 1.0.0
===========

MISC
-----------

* Poppr has been confirmed to work on Linux, Mac, and Windows systems with R
  3.0.0.
* Vignette `poppr_manual` now has cross-references to different sections.
* Vignette `poppr_manual` is quicker loading.

BUG FIXES
-----------

* removed alpha channel from plot for resampled values of I_A and \bar{r}_d due
  to warnings.

poppr 0.4.1
===========

NEW FEATURES
-----------

* `getfile` has a new argument, "combine", which will automatically add the path
  to the list of files, so they can be read without switching working directory.
* information printed to screen from `missingno` and `mlg.crosspop` will now be
  wrapped to 80 characters.

BUG FIXES
-----------

* `poppr` will now be able to correctly recognize GenAlEx files with both
  geographic and regional data.
* calculation of the index of association on P/A data with missing values will
  no longer return an error.

poppr 0.4
===========

BUG FIXES
-----------

* mistake in Bruvo's distance where it did not correctly check for ploidy level
  was fixed.
* `read.genalex` will be able to correctly distinguish between SNP and AFLP
  data.
* `read.genalex` can now correctly recognize regional formatting without an
  extra column.

NEW FEATURES
-----------

* `read.genalex` will now be able to take in a file that is formatted with both
  regional and geographic data.
* `genind2genalex` can now export xy coordinates into the GenAlEx format.
* `poppr_manual` vignette now contains images of example GenAlEx files.

NEW FILES
-----------

* `rootrot2.csv` is an example of a GenAlEx file formatted with regional data.

OTHER UPDATES
-----------

* function for guessing repeat lengths for Bruvo's distance moved into internal
  file.
* redundancy in `read.genalex` was removed.
* changed instructions in README

poppr 0.3.1
===========

BUG FIXES
-----------

* `read.genalex` will now give a warning whenever the input file is not comma
  delimited.

poppr 0.3
===========

NEW FUNCTIONS
-----------

* `poppr.msn` will draw a minimum spanning network for any distance matrix
  derived from your data set.

NEW FEATURES
-----------

* vignette now has sections describing `poppr.msn`, `diss.dist`, `greycurve`,
  and a section discussing how to export graphics.

BUG FIXES
-----------

* The graphs output by `poppr` and `ia` will now display \bar{r}_d instead of
  \bar{r}_D.
* `bruvo.boot` now has a dedicated `quiet` argument.

poppr 0.2.2
===========

NEW FEATURES
-----------

* index of association distributions will now feature a rug plot at the bottom
  as a better way to visualize the distribution of the index of association from
  the shuffled data sets.

poppr 0.2.1
===========

NEW FUNCTIONS
-----------

* `diss.dist` will produce a distance matrix based on discreet distances.
* `greycurve` will produce a grey scale adjusted to user-supplied parameters.
  This will be useful for future minimum spanning network functions.

NEW FEATURES
-----------

* `bruvo.msn` can now adjust the edge grey level to be weighted toward either
  closely or distantly weighted individuals.
* `bruvo.msn` will now return a list giving the user the graph with all of the
  color, label, and weight properties so that they can plot it themselves. The
  legend arguments are also returned.

BUG FIXES
-----------

* fixed shufflepop so that it will now shuffle PA markers with a specific method
* fixed warning message mistakes in clonecorrect function.

poppr 0.2
=========
NEW FEATURES
-----------

* Added NEWS file and will now be incrementing version number (3/15/2013)


poppr 0.1
=========
* First development version of poppr (2012 - 3/2013)
