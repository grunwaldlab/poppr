context("Round Robin Tests")

x <- structure(c("B", "C", "B", "C", "B", "C", "A", "C", "A", "A", 
"C", "B", "C", "A", "A"), .Dim = c(5L, 3L))
suppressWarnings(x <- df2genind(x, ploidy = 1))
rrx_m <- rrmlg(x)
mlg_truth <- structure(c(3L, 2L, 3L, 1L, 1L, 2L, 4L, 2L, 3L, 1L, 2L, 3L, 2L, 
                         3L, 1L), 
                       .Dim = c(5L, 3L), 
                       .Dimnames = list(NULL, c("L1", "L2", "L3"))
                       )
rrx_f <- rraf(x)
freq_truth <- structure(list(L1 = structure(c(0.333333333333333, 0.666666666666667), 
                                            .Names = c("L1.B", "L1.C")), 
                             L2 = structure(c(0.25, 0.75), 
                                            .Names = c("L2.C", "L2.A")), 
                             L3 = structure(c(0.333333333333333, 0.333333333333333, 
                                              0.333333333333333), 
                                            .Names = c("L3.C", "L3.B", "L3.A"))), 
                        .Names = c("L1", "L2", "L3")
                        )

test_that("rrmlg produces the correct mlgs", {
  expect_equivalent(rrx_m, mlg_truth)
})

test_that("rrmlg will not work on genlight objects", {
  skip_on_cran()
  expect_error(rrmlg(glSim(10, 10, 10, parallel = FALSE)))
})

test_that("rrmlg can take loci objects", {
  skip_on_cran()
  expect_equivalent(rrmlg(pegas::as.loci(x)), rrx_m)
})

test_that("rraf produces correct allele frequencies", {
  skip_on_cran()
  expect_equivalent(rrx_f, freq_truth)
  expect_equivalent(rraf(x, res = "vector"), unlist(freq_truth))
})

test_that("rare_allele_correction gets angry if called from global environment", {
  skip_on_cran()
  expect_warning(do.call(poppr:::rare_allele_correction, 
                         args = list(rrx_f, rrx_m), envir = .GlobalEnv))
})

test_that("correction is properly applied in rraf", {
  skip_on_cran()
  data(monpop)
  mll(monpop) <- "original"
  monrrmlg    <- 1/colSums(!apply(rrmlg(monpop), 2, duplicated))
  
  monc       <- rraf(monpop)                    # default
  monnc      <- rraf(monpop, correction = FALSE)# No correction
  monc_sto   <- rraf(monpop, sum_to_one = TRUE) # summed to one
  monc_mlg   <- rraf(monpop, d = "mlg")         # by multilocus genotype
  monc_rrmlg <- rraf(monpop, d = "rrmlg")       # by round-robin multilocus genotype
  monc_m     <- rraf(monpop, m = 0.5)           # assume diploid (not true here, though)
  monc_e     <- rraf(monpop, e = pi)            # A ridiculous correction
  monc_e2    <- rraf(monpop, e = pi, d = "rrmlg", m = 0.5) # e overrides

  SER <- function(x) x$SER["SER.147"]
  SED <- function(x) x$SED["SED.187"]
  
  # Correction can be turned off and on
  monc_vec  <- vapply(monc, sum, numeric(1))
  monnc_vec <- vapply(monnc, sum, numeric(1))
  expect_gt(sum(monc_vec), sum(monnc_vec))
  expect_equivalent(sum(monnc_vec), nLoc(monpop))
  
  # sum_to_one argument augments does what it says
  monc_sum2one <- vapply(monc_sto, sum, numeric(1))
  expect_equal(sum(monc_sum2one), nLoc(monpop))
  expect_true(all(monc_vec == monc_sum2one), 
    label = paste0(paste(capture.output(
      data.frame(vec = monc_vec, sum2one = monc_sum2one)
    ), collapse = "\n"), "\n")
  )
  
  # The default is 1/n
  expect_equivalent(SER(monc), 1/nInd(monpop))
  expect_true(identical(SER(monc)[[1]], SED(monc)[[1]]))
  
  # The multiplier works
  expect_equivalent(SER(monc)/2, SER(monc_m))
  
  # Works by MLG
  expect_equivalent(SER(monc_mlg), 1/nmll(monpop))
  
  # Works by rrMLG
  expect_equivalent(SER(monc_rrmlg), monrrmlg["SER"])
  # We would not expect cross-loci corrections to be equal
  expect_false(identical(SER(monc_rrmlg)[[1]], SED(monc_rrmlg)[[1]]))
  
  # Works by custom value
  expect_equivalent(SER(monc_e), pi)
  
  # Custom value overrides anything else
  expect_equivalent(SED(monc_e), SED(monc_e2))
  
})

test_that("rraf produces a data frame", {
  skip_on_cran()
  rrx_df <- rraf(x, res = "data.frame")
  expect_is(rrx_df, "data.frame")
  expect_equivalent(names(rrx_df), c("frequency", "locus", "allele"))
  # The length of the data should be equql
  expect_equal(nrow(rrx_df), sum(lengths(rrx_f)))
  # The content should be equal
  expect_equal(rrx_df$frequency, unlist(rrx_f, use.names = FALSE))
})

test_that("rraf calculates per population", {
  skip_on_cran()
  data(Pram)
  rrx_matrix <- rraf(Pram, by_pop = TRUE)
  nout <- numeric(ncol(tab(Pram)))
  exp_matrix <- t(vapply(seppop(Pram), rraf, nout, res = "vector"))
  expect_is(rrx_matrix, "matrix")
  expect_equivalent(dim(rrx_matrix), c(nPop(Pram), ncol(tab(Pram))))
  expect_equivalent(rownames(rrx_matrix), popNames(Pram))
  # The correction is 1/N per pop. PistolRSF_OR has a pop size of 4
  expect_equivalent(rrx_matrix[nrow(rrx_matrix), ncol(rrx_matrix)], 1/4)
  expect_equal(rrx_matrix, exp_matrix)
})

test_that("rraf calculates per population when supplied with a population factor", {
  skip_on_cran()
  data(Pram)
  rrx_matrix    <- rraf(Pram, by_pop = TRUE, correction = FALSE)
  rr_pop_factor <- rraf(Pram, pop = as.character(pop(Pram)), correction = FALSE)
  rr_formula    <- rraf(Pram, pop = ~SOURCE/STATE, correction = FALSE)
  expect_equivalent(rrx_matrix, rr_pop_factor)
  expect_equivalent(rrx_matrix, rr_formula)
})

test_that("psex and pgen internals produce expected results", {
  skip_on_cran()

  # setup -------------------------------------------------------------------

  # From Parks and Werth, 1993
  x <- "
   Hk Lap Mdh2 Pgm1 Pgm2 X6Pgd2
  54 12 12 12 23 22 11
  36 22 22 11 22 33 11
  10 23 22 11 33 13 13"
    
  xtab <- read.table(text = x, header = TRUE, row.names = 1)
  xgid <- df2genind(xtab[rep(rownames(xtab), as.integer(rownames(xtab))), ], ncode = 1)
  afreq <- c(Hk.1 = 0.167, Hk.2 = 0.795, Hk.3 = 0.038, 
             Lap.1 = 0.190, Lap.2 = 0.798, Lap.3 = 0.012,
             Mdh2.0 = 0.011, Mdh2.1 = 0.967, Mdh2.2 = 0.022,
             Pgm1.2 = 0.279, Pgm1.3 = 0.529, Pgm1.4 = 0.162, Pgm1.5 = 0.029,
             Pgm2.1 = 0.128, Pgm2.2 = 0.385, Pgm2.3 = 0.487,
             X6Pgd2.1 = 0.526, X6Pgd2.2 = 0.051, X6Pgd2.3 = 0.423)
  freqs   <- afreq[colnames(tab(xgid))]
  mfreq   <- matrix(freqs, nrow = 1, dimnames = list(NULL, names(freqs)))
  pNotGen <- psex(xgid, by_pop = FALSE, freq = freqs, G = 45, method = "single")
  pGen   <- exp(rowSums(pgen(xgid, by_pop = FALSE, freq = freqs)))
  pGenm  <- exp(rowSums(pgen(xgid, by_pop = FALSE, freq = mfreq)))
  res    <- matrix(c(unique(pGen), unique(pNotGen)), ncol = 2)
  expected_result <- structure(c(4.14726753733799e-05, 0.00192234266749449, 0.000558567714432373, 
                                 0.00186456862059092, 0.0829457730256321, 0.024829127718338), 
                               .Dim = c(3L,2L))
  
  # tests -------------------------------------------------------------------
  # Testing that pgen can take both a matrix and vector of af
  expect_equal(pGen, pGenm)
  expect_error(pgen(xgid, by_pop = FALSE, freq = mfreq[, -1, drop = FALSE]), "frequency matrix")
  expect_error(pgen(xgid, by_pop = FALSE, freq = mfreq[, -1, drop = TRUE]), "frequencies")
  # Testing single encounter
  expect_equivalent(res, expected_result)
})


test_that("psex produces a vector", {
  skip_on_cran()
  data(Pram)
  psexpram <- psex(Pram)
  expect_is(psexpram, "numeric")
  expect_equal(length(psexpram), nInd(Pram))
})

test_that("psex can use the multiple method (not testing for accuracy)", {
  skip_on_cran()
  data(Pram)
  psexpram <- psex(Pram, method = "multiple")
  expect_is(psexpram, "numeric")
  expect_equal(length(psexpram), nInd(Pram))
})
test_that("psex can take population factors", {
  skip_on_cran()
  data(Pram)
  psexpram <- psex(Pram)
  psexpop  <- psex(Pram, pop = as.character(pop(Pram)))
  psexform <- psex(Pram, pop = ~SOURCE/STATE)
  expect_equivalent(psexpram, psexpop)
  expect_equivalent(psexpram, psexform)
})

test_that("pgen produces a matrix", {
  skip_on_cran()
  data(Pram)
  pgenlog <- pgen(Pram, by_pop = FALSE)
  expect_is(pgenlog, "matrix")
  expect_equivalent(exp(pgenlog), pgen(Pram, by_pop = FALSE, log = FALSE))
})

test_that("psex can take population factors", {
  skip_on_cran()
  data(Pram)
  pgenpram <- pgen(Pram)
  pgenpop  <- pgen(Pram, pop = as.character(pop(Pram)))
  pgenform <- pgen(Pram, pop = ~SOURCE/STATE)
  expect_equivalent(pgenpram, pgenpop)
  expect_equivalent(pgenpram, pgenform)
})

test_that("pgen and psex work for haploids", {
  skip_on_cran()
  data(monpop)
  expect_is(psex(monpop, by_pop = FALSE), "numeric")
})

test_that("pgen can't work with polyploids", {
  skip_on_cran()
  data(Pinf)
  expect_error(pgen(Pinf))
})

test_that("psex is accurate", {
  skip_on_cran()
  skip_if_not_installed("RClone")
  # Setup
  # 
  data("zostera", package = "RClone")
  rzos <- RClone::convert_GC(zostera[, -c(1:3)], num = 3)
  pzos <- as.genclone(df2genind(zostera[, -c(1:3)], ncode = 3, pop = zostera[, 1],
                      strata = zostera[, 1, drop = FALSE]))
  extract_psex <- function(x){
    as.numeric(x$psex)
  }
  #
  # Testing basic equality, no population
  rzpsex <- RClone::psex(rzos, RR = TRUE) %>% 
    extract_psex %>% 
    setNames(indNames(pzos))
  pzpsex <- poppr::psex(pzos, by_pop = FALSE, method = "multiple")
  nas    <- is.na(rzpsex)
  
  expect_equal(rzpsex[!nas], pzpsex[!nas])
  #
  # Testing equality with population
  rzpsex <- RClone::psex(rzos, RR = TRUE, vecpop = pop(pzos)) %>% 
    lapply(extract_psex) %>%
    unlist() %>%
    setNames(indNames(pzos))
  gtab    <- table(pop(pzos))
  pzpsex  <- poppr::psex(pzos, by_pop = TRUE, method = "multiple")
  pzpsex2 <- poppr::psex(pzos, by_pop = TRUE, G = gtab, method = "multiple")
  expect_error(poppr::psex(pzos, by_pop = TRUE, G = as.vector(gtab), method = "multiple"))
  nas <- is.na(rzpsex)
  expect_equal(rzpsex[!nas], pzpsex[!nas])
  expect_equal(rzpsex[!nas], pzpsex2[!nas])
})
