context("Poppr table tests")

data(Aeut, package = "poppr")
data(partial_clone, package = "poppr")
strata(Aeut) <- other(Aeut)$population_hierarchy[-1]
A.tab <- poppr(Aeut, quiet = TRUE)
afile <- system.file("files/rootrot.csv", package = "poppr")
sims <- system.file("files/simulated.dat", package = "poppr")

Aeut_comparison <- structure(list(Pop = structure(1:3, 
    .Label = c("Athena", "Mt. Vernon", "Total"), 
    class = "factor"), 
    N = c(97, 90, 187), 
    MLG = c(70, 50, 119), 
    eMLG = c(65.9808393377379, 50, 68.4525682280634), 
    SE = c(1.24567688244012, 0, 2.98857840353058), 
    H = c(4.06272002528149, 3.66843094399907, 4.55798828426928), 
    G = c(42.1928251121076, 28.7234042553191, 68.9723865877712), 
    lambda = c(0.976299287915825, 0.965185185185185, 0.985501444136235), 
    E.5 = c(0.721008688944842, 0.725926650260449, 0.720112175857993), 
    Hexp = c(0.169711892488954, 0.158016764758338, 0.365053352719387), #c(0.408951651286696, 0.33656220322887, 0.726202391505946),
    Ia = c(2.90602921191748, 13.3024309367662, 14.3707995986407), 
    rbarD = c(0.0723700801886747, 0.281642324983496, 0.270617053778004), 
    File = structure(c(1L, 1L, 1L), class = "factor", .Label = "rootrot.csv")), 
    .Names = c("Pop", "N", "MLG", "eMLG", "SE", "H", "G", "lambda", "E.5", "Hexp", "Ia", "rbarD", "File"), 
    row.names = c(NA, -3L), 
    class = c("popprtable", "data.frame"))

pc_comparison <- structure(list(Pop = structure(1:5, .Label = c("1", "2", "3", 
"4", "Total"), class = "factor"),
    N = c(13, 13, 12, 12, 50), 
    MLG = c(10, 12, 11, 9, 26),
    eMLG = c(9.46153846153846, 11.1538461538462, 11, 9, 9.93701353621932), 
    SE = c(0.498518515262143, 0.36080121229411, 0, 0, 1.12881059579593),
    H = c(2.24503527412618, 2.45831132968308, 2.36938211969468, 2.09472904752765, 3.07152395656842), 
    G = c(8.89473684210526, 11.2666666666667, 10.2857142857143, 7.2, 17.8571428571429 ), 
    lambda = c(0.887573964497041, 0.911242603550296, 0.902777777777778, 
             0.861111111111111, 0.944), 
    E.5 = c(0.935312405238733, 0.960842907662783, 0.958200460752105, 0.870390481833875, 0.819311895784525), 
    Hexp = c(0.617846153846154, 0.614461538461539, 0.485144927536232, 0.505072463768116, 
             0.559212121212121), #c(0.849679487179487, 0.83629191321499, 0.714756944444444, 0.750347222222222, 0.78384),
    Ia = c(2.1580763424628, 1.87492360969648, 1.15572679509632, 1.157153633392, 1.93513179817012),
    rbarD = c(0.243225877705591, 0.212786561587854, 0.132460530412697, 0.13328661732193, 0.217470007471919),
    File = structure(c(1L, 1L, 1L, 1L, 1L), class = "factor", .Label = "partial_clone.dat")), 
    .Names = c("Pop", "N", "MLG", "eMLG", "SE", "H", "G", "lambda", "E.5", "Hexp", "Ia", "rbarD", "File"),
    row.names = c(NA, -5L),
    class = c("popprtable", "data.frame"))

p.tab <- poppr(partial_clone, quiet = TRUE)

test_that("poppr returns expected PA values", { 
  expect_that(A.tab$Pop, is_equivalent_to(Aeut_comparison$Pop))
  expect_that(A.tab$N, equals(Aeut_comparison$N))
  expect_that(A.tab$MLG, equals(Aeut_comparison$MLG))
  expect_that(A.tab$eMLG, equals(Aeut_comparison$eMLG))
  expect_that(A.tab$SE, equals(Aeut_comparison$SE))
  expect_that(A.tab$H, equals(Aeut_comparison$H))
  expect_that(A.tab$G, equals(Aeut_comparison$G))
  expect_that(A.tab$lambda, equals(Aeut_comparison$lambda))
  expect_that(A.tab$E.5, equals(Aeut_comparison$E.5))
  expect_that(A.tab$Ia, equals(Aeut_comparison$Ia))
  expect_that(A.tab$rbarD, equals(Aeut_comparison$rbarD))
})

test_that("poppr returns expected codominant values", {
  expect_that(p.tab$Pop, is_equivalent_to(pc_comparison$Pop))
  expect_that(p.tab$N, equals(pc_comparison$N))
  expect_that(p.tab$MLG, equals(pc_comparison$MLG))
  expect_that(p.tab$eMLG, equals(pc_comparison$eMLG))
  expect_that(p.tab$SE, equals(pc_comparison$SE))
  expect_that(p.tab$H, equals(pc_comparison$H))
  expect_that(p.tab$G, equals(pc_comparison$G))
  expect_that(p.tab$lambda, equals(pc_comparison$lambda))
  expect_that(p.tab$E.5, equals(pc_comparison$E.5))
  expect_that(p.tab$Ia, equals(pc_comparison$Ia))
  expect_that(p.tab$rbarD, equals(pc_comparison$rbarD))
})

test_that("poppr perform clone correction", {
  skip_on_cran()
  res_na  <- poppr(Aeut, clonecorrect = TRUE, strata = NA, quiet = TRUE)
  res_p   <- poppr(Aeut, clonecorrect = TRUE, strata = ~Pop, quiet = TRUE)
  res_ps  <- poppr(Aeut, clonecorrect = TRUE, strata = ~Pop/Subpop, quiet = TRUE)
  res_ps2 <- poppr(Aeut, clonecorrect = TRUE, strata = ~Pop/Subpop, keep = 1:2, 
                   quiet = TRUE)
  res_naf  <- poppr(afile, clonecorrect = TRUE, strata = ~Pop, quiet = TRUE)
  res_mv  <- poppr(Aeut, clonecorrect = TRUE, strata = ~Pop/Subpop, 
                   quiet = TRUE, sublist = "Mt. Vernon")
  expect_that(nrow(res_na), equals(nrow(Aeut_comparison)))
  expect_that(nrow(res_p), equals(nrow(Aeut_comparison)))
  expect_that(nrow(res_ps), equals(nrow(Aeut_comparison)))
  expect_gt(nrow(res_ps2), nrow(Aeut_comparison))

  expect_true(all(res_ps$N < Aeut_comparison$N))
  expect_true(all(res_p$N < res_ps$N))
  expect_true(all(res_na$N <= res_p$N))
  expect_identical(res_ps2[-c(1, length(res_na))], res_naf[-c(1, length(res_naf))])

  expect_true(all.equal(res_mv[, -c(4:5)], res_ps[2, -c(4:5)], check.attributes = FALSE))
})

test_that("poppr will correct missing values as mean for P/A markers", {
  skip_on_cran()
  a <- tab(Aeut)
  a[117, 18] <- NA # this locus has a mean of 0.5, happily.
  A <- df2genind(a, ind.names = indNames(Aeut), 
                 pop = pop(Aeut), ploidy = 1, 
                 type = "PA", sep = "/")
  expect_warning(pm <- poppr(A, missing = "mean", quiet = TRUE, total = FALSE), "adegenet_2.0-0")
  expect_lt(pm$rbarD[2], A.tab$rbarD[2]) # smaller ia value because of homogenization
  expect_gt(pm$MLG[2], A.tab$MLG[2])     # more MLGs because of extra missing value
})

test_that("poppr skips over sample sizes less than three", {
  skip_on_cran()
  expect_warning(plt <- poppr(partial_clone[1:8], sample = 10, quiet = TRUE), "could not be plotted")
  expect_is(plt, "popprtable")
  expect_equivalent(signif(plt$Ia, 3), c(rep(NA, 4), 0.167))
  expect_equivalent(signif(plt$rbarD, 3), c(rep(NA, 4), 0.0195))
  expect_output(print(plt, digits = 2), "0.69")
  expect_true(is.na(ia(partial_clone[1:2])[2]))
  # Presence/absence case
  a <- Aeut
  pop(a) <- c("A", "A", as.character(pop(Aeut))[-(1:2)])
  expect_true(is.na(poppr(a, sublist = "A")$rbarD[1]))
  expect_true(is.na(ia(a[1:2])[2]))
})

test_that("poppr can produce output from input file", {
  skip_on_cran()
  expect_message(out <- poppr(afile, legend = TRUE, quiet = TRUE), "Simpson")
  expect_message(outs <- poppr(sims, legend = TRUE, quiet = TRUE), "Simpson")
  expect_is(out, "popprtable")
  expect_is(outs, "popprtable")
})

test_that("poppr.all works on list of files", {
  skip_on_cran()
  expect_output(out <- poppr.all(list(afile, pc = partial_clone)), " | File: rootrot.csv")
})

test_that("poppr without sampling works on data without populations", {
  skip_on_cran()
  pop(partial_clone) <- NULL
  out <- poppr(partial_clone, sample = 9, quiet = TRUE)
  expect_is(out, "data.frame")
  expect_equal(nrow(out), 1L)
})