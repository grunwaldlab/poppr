## ----echo=FALSE, warning=FALSE-------------------------------------------
knitr::opts_knit$set(out.format = "latex")
thm <- knitr::knit_theme$get("acid")
knitr::knit_theme$set(thm)
knitr::opts_chunk$set(concordance=TRUE)
knitr::opts_chunk$set(size = 'small', message = FALSE, warning = FALSE)
knitr::opts_chunk$set(out.width = '0.5\\linewidth', fig.align = "center", fig.show = 'asis')

## ----poppr_funk, eval = TRUE, echo = FALSE-------------------------------
print_command <- function(funk){
  fargs <- formals(funk)

  lapply(names(fargs), function(arg_name, fargs){
    arg <- fargs[[arg_name]]
    if (missing(arg)){
      fargs[[arg_name]] <<- as.symbol(arg_name)
      names(fargs)[names(fargs) == arg_name] <<- ""
    }
  }, fargs)
  fargs$call <- as.symbol(funk)
  fargs <- fargs[c(length(fargs), 1:(length(fargs) - 1))]
  return(as.call(fargs))
}

## ----migration_vignette, eval = FALSE------------------------------------
## vignette("how_to_migrate", package = "poppr")

## ----popprcite, eval = TRUE, size="normalsize"---------------------------
citation(package = "poppr")

## ----install, eval=FALSE-------------------------------------------------
## install.packages("poppr", dependencies=TRUE)

## ----install_devtools, eval = FALSE--------------------------------------
## install.packages("devtools")

## ----install_github, eval = FALSE----------------------------------------
## devtools::install_github("grunwaldlab/poppr")

## ----install_devel, eval = FALSE-----------------------------------------
## devtools::install_github("grunwaldlab/poppr@devel")

## ----install_depend, eval=FALSE, tidy = FALSE----------------------------
## pkgs <- c("adegenet", "pegas", "vegan", "ggplot2", "phangorn", "ape",
##           "igraph", "reshape2", "dplyr", "shiny")
## install.packages(pkgs)

## ----install_source, eval=FALSE------------------------------------------
## install.packages("/path/to/poppr.tar.gz", type="source", repos=NULL)

## ----echo=FALSE----------------------------------------------------------
x <- list(files="/path/to/R/poppr/files/rootrot.csv", path="/path/to/R/poppr/files")

## ----message = TRUE------------------------------------------------------
library("poppr")

## ----getfilefunk, eval=FALSE---------------------------------------------
## x <- getfile()

## ----getfilex------------------------------------------------------------
x

## ----firstpoppr, eval=FALSE----------------------------------------------
## myData  <- read.genalex(x$files)
## myData

## ----echo = FALSE--------------------------------------------------------
options(width=90)
myData <- read.genalex(system.file("files/rootrot.csv", package="poppr"))

## ----aflp----------------------------------------------------------------
popdata <- poppr(myData)

## ----firstpoppr2---------------------------------------------------------
popdata

## ----microsave, eval=FALSE-----------------------------------------------
## library("poppr")
## data(microbov)
## strata(microbov) <- data.frame(other(microbov)) # set the strata
## microbov
## genind2genalex(microbov, file = "~/Desktop/microbov.csv")

## ----microsave_message, echo=FALSE---------------------------------------
cat("Extracting the table ... Writing the table to ~/Desktop/microbov.csv ... Done.")

## ----read.genalex_command, echo = FALSE, comment = NA, size = 'normalsize'----
funk <- "read.genalex"
print_command(funk)

## ----system_file_genalex, eval=FALSE-------------------------------------
## system.file("files/rootrot.csv", package="poppr")

## ----system_file_echo, echo=FALSE----------------------------------------
paste("/path/to/R/library/poppr/files/rootrot.csv")

## ----read.genalex_ex-----------------------------------------------------
rootrot <- read.genalex(system.file("files/rootrot.csv", package="poppr"))

## ----read.genalex_ex2----------------------------------------------------
rootrot

## ----genind2genalex_command, echo = FALSE, comment = NA, size = 'normalsize'----
funk <- "genind2genalex"
print_command(funk)

## ----genind2genalex, eval=FALSE------------------------------------------
## genind2genalex(rootrot, "~/Desktop/rootrot.csv")

## ----genind2genalex_cat, echo=FALSE--------------------------------------
cat("Extracting the table ... Writing the table to ~/Desktop/rootrot.csv ... Done.\n")

## ----nancyxy-------------------------------------------------------------
data(nancycats)
other(nancycats)$xy

## ----genind2genalex_nancy, eval=FALSE------------------------------------
## genind2genalex(nancycats, "~/Desktop/nancycats_pop_xy.csv", geo = TRUE)

## ----genind2genalex_cat2, echo=FALSE-------------------------------------
cat("Extracting the table ... Writing the table to ~/Desktop/nancycats_pop_xy.csv ... Done.\n")

## ----nancy_grow_xy-------------------------------------------------------
nan2 <- nancycats
other(nan2)$xy <- other(nan2)$xy[pop(nan2), ]
tail(other(nan2)$xy)

## ----genind2genalex_nancy_grow, eval=FALSE-------------------------------
## genind2genalex(nan2, "~/Desktop/nancycats_inds_xy.csv", geo = TRUE)

## ----genind2genalex_cat3, echo=FALSE-------------------------------------
cat("Extracting the table ... Writing the table to ~/Desktop/nancycats_inds_xy.csv ... Done.\n")

## ----ex_data_picture, echo = FALSE, fig.width=14, fig.height=4, out.width="7in", out.height = "2in"----
library("adegenet")
df <- data.frame(list(locus1=c("101/101", "102/103", "102/102"),
                      locus2=c("201/201", "202/203", "203/204"),
                      locus3=c("301/302", "301/303", "304/305")))
dat <- tab(df2genind(df, sep="/"))
tdat <- dat
tdat[] <- 1
x <- barplot(tdat, axes = FALSE, axisnames = FALSE)
barplot(rep(3, 12), col = rep(rainbow(3, alpha = 0.5), 3:5), axes = FALSE, add = TRUE)
axis(2, at = 1:3 - 0.5, labels = 1:3, tick = FALSE)
axis(3, at = x, labels = colnames(tdat), tick = FALSE)
axis(1, at = c(2, 6.125, 11.5), labels = names(df), tick = FALSE)

## ----ex_genind, echo=FALSE-----------------------------------------------
dat

## ----other_replen--------------------------------------------------------
data(nancycats) # Load the data
other(nancycats) # geographical coordinates
repeats <- rep(2, nLoc(nancycats)) #nLoc = number of loci
repeats
other(nancycats)$repeat_lengths <- repeats
other(nancycats) # two items named xy and repeat_lengths

## ----as.genclone_command, echo = FALSE, comment = NA, size = 'normalsize'----
funk <- "as.genclone"
print_command(funk)

## ----show_genind---------------------------------------------------------
library("poppr")
data(Aeut)
Aeut

## ----add_strata_Aeut-----------------------------------------------------
strata(Aeut) <- other(Aeut)$population_hierarchy[-1]
Aeut

## ----as_genclone---------------------------------------------------------
agc <- as.genclone(Aeut)
agc

## ----genclone_compare----------------------------------------------------
c(is.genind(Aeut), is.genclone(Aeut), is.genind(agc), is.genclone(agc))

# Adegenet functions work the same, too
c(nInd(Aeut), nInd(agc))

## ----genclone2genind-----------------------------------------------------
genclone2genind(agc)

## ----pinf_table, fig.width = 5, fig.height = 15--------------------------
data(Pinf)
Pinf
ptab <- info_table(Pinf, type = "ploidy", plot = TRUE)

## ----pinf_show-----------------------------------------------------------
tail(tab(Pinf[loc = locNames(Pinf)[9:10]]))

## ----pinf_show_df--------------------------------------------------------
Pinfdf <- genind2df(Pinf, sep = "/")
tail(Pinfdf[10:11])

## ----pinf_recode---------------------------------------------------------
Pinf_rc <- recode_polyploids(Pinf, newploidy = TRUE)
Pinf_rc # Notice that the new ploidy is accounted for.
tail(tab(Pinf_rc[loc = locNames(Pinf_rc)[9:10]]))

## ----pinf_recode_df------------------------------------------------------
Pinfrcdf <- genind2df(Pinf_rc, sep = "/")
tail(Pinfrcdf[10:11])

## ----pinf_rerecode-------------------------------------------------------
tail(tab(recode_polyploids(Pinf_rc[loc = locNames(Pinf_rc)[9:10]], addzero = TRUE)))

## ----missingo_command, echo = FALSE, comment = NA, size = 'normalsize'----
funk <- "missingno"
print_command(funk)

## ----initializing_poppr, out.width=".8\\linewidth", fig.height = 6, fig.width = 10----
library("poppr")
data(nancycats)
info_table(nancycats, plot = TRUE)

## ----nancy_indiv---------------------------------------------------------
tab(nancycats)[1:5, 8:13]

## ----missingno_exclude---------------------------------------------------
nanloci <-  missingno(nancycats, "loci")
nangeno <-  missingno(nancycats, "geno")
tab(nanloci)[1:5, 8:13]

## ----missingno_loci------------------------------------------------------
nInd(nanloci)     # Individuals
locNames(nanloci) # Names of the loci

## ----missingno_geno------------------------------------------------------
tab(nangeno)[1:5, 8:13]
nInd(nangeno)     # Individuals
locNames(nangeno) # Names of the loci

## ----popsub_command, echo = FALSE, comment = NA, size = 'normalsize'-----
funk <- "popsub"
print_command(funk)

## ----load_H3N2-----------------------------------------------------------
data("H3N2", package = "adegenet")
strata(H3N2) <- data.frame(other(H3N2)$x)
H3N2

## ----popsub_sublist------------------------------------------------------
setPop(H3N2) <- ~country
popNames(H3N2) # Only two countries from North America.
v_na <- popsub(H3N2, sublist = c("USA", "Canada"))
popNames(v_na)

## ----popsub_sizes--------------------------------------------------------
c(NorthAmerica = nInd(v_na), Total = nInd(H3N2))

## ----popsub_blacklist----------------------------------------------------
v_na_minus <- popsub(H3N2, blacklist = c("USA", "Canada"))
popNames(v_na_minus)

## ----length_test---------------------------------------------------------
(nInd(v_na_minus) + nInd(v_na)) == nInd(H3N2)

## ----popsub_combine------------------------------------------------------
vsort <- sort(popNames(H3N2))[1:10]
vsort
valph <- popsub(H3N2, sublist = vsort, blacklist = c("USA", "Canada"))
popNames(valph)

## ----clonecorrect_command, echo = FALSE, comment = NA, size = 'normalsize'----
funk <- "clonecorrect"
print_command(funk)

## ----clonecorrect--------------------------------------------------------
data(Aeut)
strata(Aeut) <- data.frame(other(Aeut)$population_hierarchy[-1])
Aeut

## ----clonecorrect_genclone-----------------------------------------------
aphan <- as.genclone(Aeut)
nameStrata(Aeut) <- ~field/sample

## ----clonecorrect2-------------------------------------------------------
clonecorrect(aphan,  strata = ~Pop/Subpop)
# Your turn: Use the same stratification and use combine = TRUE and then
# keep = 1:2. Is there any difference?

## ----clonecorrect3-------------------------------------------------------
clonecorrect(aphan, strata = ~Pop)

## ----clonecorrectx-------------------------------------------------------
clonecorrect(aphan, strata = NA)

## ----shufflepop_command, echo = FALSE, comment = NA, size = 'normalsize'----
funk <- "shufflepop"
print_command(funk)

## ----shuffle_bd----------------------------------------------------------
data(nancycats)
nan1 <- popsub(nancycats, 1)
reps <- rep(2, 9) # Assuming dinucleotide repeats.
observed <- mean(bruvo.dist(nan1, replen = reps))
observed

## ----shuffle_bd_replicate_dummy, eval=FALSE------------------------------
## set.seed(9999)
## bd.test <- replicate(999, mean(bruvo.dist(shufflepop(nan1, method = 2), replen = reps)))

## ----shuffle_bd_replicate, echo = FALSE----------------------------------
load("nancybruvo.rda")

## ----bd_histogram, resolution = 300--------------------------------------
hist(bd.test, xlab = "Bruvo's Distance", main = "Average Bruvo's distance over 999 randomizations")
abline(v = observed, col = "red")
legend('topleft', legend="observed", col="red", lty = 1)

## ----informloci_command, echo = FALSE, comment = NA, size = 'normalsize'----
funk <- "informloci"
print_command(funk)

## ----inform.H3N2.1, eval = FALSE-----------------------------------------
## H.five <- informloci(H3N2, cutoff = 0.05)

## ----inform.H3N2.2, echo = FALSE, message = TRUE-------------------------
#res <- c(157,177,233,243,262,267,280,303,313,327,357,382,384,399,412,418,424,425,429,433,451,470,529,546,555,557,564,576,592,595,597,602,612,627,642,647,648,654,658,663,667,681,717,806,824,837,882)
# cat("cutoff value: 5 percent ( 95 individuals ).\n","47 uninfomative loci found:", res, fill = 80)
msg <- readLines("msg.txt")
for (i in msg) message(i)

## ----inform.nancy, message = TRUE----------------------------------------
data(nancycats)
naninform <- informloci(nancycats, cutoff = 0.05)

## ----eval = FALSE--------------------------------------------------------
## vignette("mlg", package = "poppr")

## ----view_mlg------------------------------------------------------------
H3N2

## ----mlg_genind----------------------------------------------------------
H3N2_mlg <- mlg(H3N2)
H3N2_mlg

## ----mlg.crosspop_command, echo = FALSE, comment = NA, size = 'normalsize'----
funk <- "mlg.crosspop"
print_command(funk)

## ----crosspop, eval=FALSE------------------------------------------------
## setPop(H3N2) <- ~country
## v.dup <- mlg.crosspop(H3N2, quiet=TRUE)

## ----crosspopout, echo=FALSE---------------------------------------------
setPop(H3N2) <- ~country
v.dup <- structure(list(MLG.3 = structure(c(4L, 8L), .Names = c("USA",
"Denmark")), MLG.9 = structure(c(1L, 13L, 1L, 1L), .Names = c("Japan",
"USA", "Finland", "Denmark")), MLG.31 = structure(c(2L, 7L), .Names = c("Japan",
"Canada")), MLG.75 = structure(c(2L, 8L, 2L, 1L, 6L, 2L, 1L,
1L), .Names = c("Japan", "USA", "Finland", "Norway", "Denmark",
"Austria", "Russia", "Ireland")), MLG.80 = structure(c(1L, 1L
), .Names = c("USA", "Denmark")), MLG.86 = structure(3:4, .Names = c("Denmark",
"Austria")), MLG.95 = structure(c(1L, 1L), .Names = c("USA",
"Bangladesh")), MLG.97 = structure(c(1L, 5L, 1L, 1L), .Names = c("USA",
"Austria", "Bangladesh", "Romania")), MLG.104 = structure(1:2, .Names = c("USA",
"France")), MLG.110 = structure(c(2L, 3L, 11L), .Names = c("Japan",
"USA", "China"))), .Names = c("MLG.3", "MLG.9", "MLG.31", "MLG.75",
"MLG.80", "MLG.86", "MLG.95", "MLG.97", "MLG.104", "MLG.110"))
printthings <- function(ind, x){
  cat(paste0(names(x)[ind], ":"),
      paste0("(", sum(x[[ind]])," inds)"),
      names(x[[ind]]),
      "\n")
}
invisible(lapply(1:10, printthings, v.dup))

## ----crosspop2-----------------------------------------------------------
head(v.dup)
v.num <- sapply(v.dup, length) # count the number of populations each MLG crosses.
head(v.num)

## ----mlg.table_command, echo = FALSE, comment = NA, size = 'normalsize'----
funk <- "mlg.table"
print_command(funk)

## ----mlgbar, eval = FALSE------------------------------------------------
## v.tab <- mlg.table(H3N2, plot = TRUE)
## v.tab[1:10, 1:10] # Showing the first 10 columns and rows of the table.

## ----mlgbarshow, echo = FALSE--------------------------------------------
v.tab <- mlg.table(H3N2, plot = FALSE)
v.tab[1:10, 1:10] # Showing the first 10 columns and rows of the table.

## ----mlgbarplot, results = 'hide', out.width = "0.8\\linewidth"----------
mlg.table(H3N2, sublist = "Norway", plot = TRUE)

## ----mlgrare1_dummy, eval = FALSE----------------------------------------
## setPop(H3N2) <- ~year
## summary(H3N2) # Check the data to make sure it's correct.

## ----mlgrare1, echo = FALSE----------------------------------------------
setPop(H3N2) <- ~year
res <- c("", " # Total number of genotypes:  1903 ", "", " # Population sample sizes:  ",
         "2002 2003 2004 2005 2006 ", " 158  415  399  469  462 ", "",
         " # Number of alleles per locus:  ", "  6  17  39  42  45  51  60  72  73  90 108 123 129 134 145 148 149 157 168 171 177 225 ",
         "  3   3   4   2   4   2   3   2   4   3   4   2   4   3   2   2   3   3   2   2   3   3 ",
         "233 243 247 262 267 280 303 313 317 327 334 345 351 357 376 382 384 391 396 399 412 418 ",
         "  3   2   2   2   2   2   2   2   2   2   2   4   4   3   3   3   4   2   2   2   4   3 ",
         "424 425 429 430 433 434 435 451 463 464 468 470 476 483 490 517 529 546 555 557 561 562 ",
         "  2   3   4   2   3   2   3   2   2   2   4   2   2   2   2   2   2   2   4   4   4   3 ",
         "564 566 576 577 578 582 592 594 595 597 600 602 604 612 627 642 647 648 654 658 663 664 ",
         "  3   2   3   4   3   2   3   3   3   3   2   3   2   4   2   3   2   2   3   3   3   3 ",
         "666 667 673 674 676 679 681 685 717 763 806 807 824 837 882 897 906 910 915 929 933 936 ",
         "  2   2   2   2   3   2   3   2   3   2   3   2   3   3   2   2   2   3   2   2   2   3 ",
         "939 940 957 961 962 963 966 967 969 973 975 977 978 979 980 ",
         "  3   3   2   2   3   3   3   3   4   2   3   3   4   3   2 ",
         "", " # Number of alleles per population:  ", "2002 2003 2004 2005 2006 ",
         " 203  255  232  262  240 ", "", " # Percentage of missing data:  ",
         "[1] 2.363426", "", " # Observed heterozygosity:  ", "[1] 0",
         "", " # Expected heterozygosity:  ", "[1] 0")
cat(res, sep = "\n")

## ----mlgrare2, eval=FALSE, tidy=FALSE------------------------------------
## library("vegan")
## H.year <- mlg.table(H3N2, plot = FALSE)
## rarecurve(H.year, ylab="Number of expected MLGs", sample=min(rowSums(H.year)),
##           border = NA, fill = NA, font = 2, cex = 1, col = "blue")

## ----mlgrareplot, echo=FALSE, results='hide'-----------------------------
# library("vegan")
H.year <- mlg.table(H3N2, plot = FALSE)
vegan::rarecurve(H.year, ylab="Number of expected MLGs", sample=min(rowSums(H.year)), border = NA, fill = NA, font = 2, cex = 1, col = "blue")

## ----subcross------------------------------------------------------------
setPop(H3N2) <- ~country
UGNN.list <- c("United Kingdom", "Germany", "Netherlands", "Norway")
UGNN <- mlg.crosspop(H3N2, sublist=UGNN.list, indexreturn=TRUE)

## ----subtable------------------------------------------------------------
UGNN # Note that we have three numbers here. This will index the columns for us.
UGNN.list # And let's not forget that we have the population names.
v.tab[UGNN.list, UGNN]

## ----mlg.vector_first----------------------------------------------------
v.vec <- mlg.vector(H3N2)
str(v.vec) # Analyze the structure.

## ----mlg.vector_second---------------------------------------------------
length(unique(v.vec)) # count the number of MLGs
H3N2 # equal to the first number in this output.

## ----mlg.vector_match----------------------------------------------------
UGNN # Show what we are looking for
UGNN_match <- v.vec %in% UGNN
table(UGNN_match) # How many individuals matched to those three MLGs?

## ----mlg.vector_inds-----------------------------------------------------
indNames(H3N2)[UGNN_match]

## ----mlg.id--------------------------------------------------------------
H3N2.id <- mlg.id(H3N2)
H3N2.id[as.character(UGNN)]

## ----mlgsub_flag, eval=FALSE---------------------------------------------
## mlg.table(H3N2, mlgsub = UGNN)

## ----mlgsub_flagshow, results='hide', echo=FALSE, out.width = "0.8\\linewidth"----
mlg.table(H3N2, mlgsub = UGNN)

## ----Aeut_MLG_Aeut, fig.width = 7, fig.height = 3, out.width = "0.8\\linewidth"----
library("poppr")
library("ggplot2")
data(Aeut)
Aeut.tab <- mlg.table(Aeut)
p <- last_plot()

## ----Aeut_MLG_title, fig.width = 7, fig.height = 3, out.width = "0.8\\linewidth"----
myTitle <- expression(paste(italic("Aphanomyces euteiches"), " multilocus genotype distribution"))
(pt <- p +
   ggtitle(myTitle) +
   xlab("Multilocus genotype")) # We can label the x axis, too

## ----Aeut_MLG_theme, fig.width = 7, fig.height = 3, out.width = "0.8\\linewidth"----
(ptt <- pt + theme_bw())

## ----Aeut_MLG_theme_axis, fig.width = 7, fig.height = 3, out.width = "0.8\\linewidth"----
(ptta <- ptt + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()))

## ----Aeut_silly, fig.width = 7, fig.height = 3, out.width = "0.8\\linewidth"----
(ptttaf <- ptta + aes(fill = count))

## ----Aeut_MLG_data-------------------------------------------------------
head(p$data)

## ----ggsave1, eval=FALSE-------------------------------------------------
## data(nancycats) # Load the data set.
## poppr(nancycats, sample=999) # Produce a single plot.
## ggsave("nancycats.pdf")

## ----png_save, eval=FALSE------------------------------------------------
## data(nancycats)
## ####
## png("nancy_pair%02d.png", width = 14, height = 14, units = "in", res = 300)
## poppairs <- lapply(seppop(nancycats), pair.ia, limits = c(-0.25, 1))
## dev.off()
## ####

## ----pdf_save, eval=FALSE------------------------------------------------
## pdf("nancy_pair.png", width = 14, height = 14, compress = FALSE)
## poppairs <- lapply(seppop(nancycats), pair.ia, limits = c(-0.25, 1))
## dev.off()

