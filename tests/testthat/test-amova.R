{
data("Aeut", package = "poppr")
strata(Aeut) <- other(Aeut)$population_hierarchy[-1]

agc    <- as.genclone(Aeut)
Athena <- popsub(agc, "Athena")
res    <- poppr.amova(Aeut, ~Pop/Subpop, quiet = TRUE)
rescc  <- poppr.amova(Aeut, ~Pop/Subpop, quiet = TRUE, clonecorrect = TRUE)
resper <- c(70.0067859292295, 8.40748251295027, 21.5857315578203, 100)
ressig <- c(11.0634458464745, 1.3286673034988, 
            3.41127747798475, 15.8033906279581)
resccper <- c(66.7776803325885, 6.10535452127678, 27.1169651461347, 100)
resccsig <- c(10.4131525344064, 0.952054452775842, 
              4.22855500416477, 15.593761991347)

meirmans_polyploid <- data.frame(stringsAsFactors=FALSE,
                                 Locus827 = c("053053053097097060", "053060060025025025",
                                              "053053053053053060", "053053053053053025",
                                              "053053053053060025", "053053053025025064", "053053053053097097",
                                              "053053053053053097", "053053053060060025", "053053053053025025",
                                              "053025025025025083", "053053025025083083", "053097097025083059",
                                              "053053025025025083", "053097025025083083", "053097025083083083",
                                              "053053097025025083", "053053025083083083", "053053053053097083",
                                              "053053053053097097", "053053053053053041", "053053053053053053",
                                              "053053053053053053", "053053053083059041",
                                              "053053053041041041", "053053053053053041", "053053053083041041",
                                              "053053053053053064", "053053053053097064", "053053053064083041",
                                              "064064064076076008", "053064064064008008", "025064064064064076",
                                              "064064064064064076", "053053025064064064", "053064076008008008",
                                              "064064064064064008", "053064064064064064", "053064064076008008",
                                              "064064064076076076", "064064064064064008", "053064008008008008",
                                              "053053064064064064", "064064064064064008",
                                              "064064064064064008", "053064064064008008", "053053053064064029",
                                              "064064064064008008", "053053064064064064", "064064064008008008",
                                              "076008008008008008", "064064064076076019", "076076008008008075",
                                              "064076008008019075", "076076076076008008", "076076008008008008",
                                              "053064064076008075", "053064076008008075", "064064064008008008",
                                              "064064076008008029"),
                                 Locus991 = c("098010010010067066", "010010010067066066",
                                              "098010010067066066", "098098010010010067",
                                              "098010010010066066", "098010010010067066", "010010010010067026",
                                              "010010010067067066", "098010010010010066", "010010010010067066",
                                              "010010041041041008", "010010010010067067", "010010067067067008",
                                              "010010010010041008", "010010010010067041", "010010010010067041",
                                              "010010010067067041", "010067067066066041", "010067067067008008",
                                              "010067066041041008", "010010066041041008", "098010010010010041",
                                              "010066066041041041", "010010010041041041",
                                              "010010066066041041", "010010010041041008", "010010066066066008",
                                              "010066066041041008", "010010010041041008", "010010010066041008",
                                              "067008008008008088", "067008088036036007", "067067067007007016",
                                              "067008007007007007", "008007090090090090", "008008007007007016",
                                              "067067067088090090", "067067008008007090", "067088007007024024",
                                              "067008008008007090", "066066024024024043", "066066066066066024",
                                              "066066066024024024", "066066066024024043",
                                              "067066066066066066", "066066066066024043", "066066066066066024",
                                              "067066066066024043", "066066066036036043", "066066066016024043",
                                              "066016016016016016", "067067066008016016", "067067066066016016",
                                              "066066016016016024", "067066066016024024", "067067067066016016",
                                              "067066066016016016", "067066016016016016", "067066066066016024",
                                              "067067067066066066"),
                                 Individual = c("Ind001_0001", "Ind001_0015", "Ind001_0018", "Ind001_0026",
                                                "Ind001_0061", "Ind001_0066", "Ind001_0071", "Ind001_0078",
                                                "Ind001_0080", "Ind001_0098", "Ind002_0005", "Ind002_0027",
                                                "Ind002_0035", "Ind002_0047", "Ind002_0048", "Ind002_0056",
                                                "Ind002_0077", "Ind002_0083", "Ind002_0090", "Ind002_0100",
                                                "Ind003_0001", "Ind003_0009", "Ind003_0010", "Ind003_0023",
                                                "Ind003_0034", "Ind003_0076", "Ind003_0084", "Ind003_0092", "Ind003_0094",
                                                "Ind003_0096", "Ind011_0018", "Ind011_0020", "Ind011_0038",
                                                "Ind011_0044", "Ind011_0047", "Ind011_0070", "Ind011_0077",
                                                "Ind011_0081", "Ind011_0095", "Ind011_0096", "Ind012_0010",
                                                "Ind012_0022", "Ind012_0023", "Ind012_0044", "Ind012_0047",
                                                "Ind012_0051", "Ind012_0061", "Ind012_0096", "Ind012_0098", "Ind012_0100",
                                                "Ind013_0005", "Ind013_0014", "Ind013_0016", "Ind013_0020",
                                                "Ind013_0029", "Ind013_0031", "Ind013_0065", "Ind013_0070",
                                                "Ind013_0085", "Ind013_0097"),
                                 pop = as.factor(c("A_1_Ind001_0001", "A_1_Ind001_0015",
                                                   "A_1_Ind001_0018", "A_1_Ind001_0026",
                                                   "A_1_Ind001_0061", "A_1_Ind001_0066", "A_1_Ind001_0071",
                                                   "A_1_Ind001_0078", "A_1_Ind001_0080", "A_1_Ind001_0098",
                                                   "A_2_Ind002_0005", "A_2_Ind002_0027", "A_2_Ind002_0035",
                                                   "A_2_Ind002_0047", "A_2_Ind002_0048",
                                                   "A_2_Ind002_0056", "A_2_Ind002_0077", "A_2_Ind002_0083",
                                                   "A_2_Ind002_0090", "A_2_Ind002_0100", "A_3_Ind003_0001",
                                                   "A_3_Ind003_0009", "A_3_Ind003_0010", "A_3_Ind003_0023",
                                                   "A_3_Ind003_0034", "A_3_Ind003_0076", "A_3_Ind003_0084",
                                                   "A_3_Ind003_0092", "A_3_Ind003_0094",
                                                   "A_3_Ind003_0096", "B_4_Ind011_0018", "B_4_Ind011_0020",
                                                   "B_4_Ind011_0038", "B_4_Ind011_0044", "B_4_Ind011_0047",
                                                   "B_4_Ind011_0070", "B_4_Ind011_0077", "B_4_Ind011_0081",
                                                   "B_4_Ind011_0095", "B_4_Ind011_0096", "B_5_Ind012_0010",
                                                   "B_5_Ind012_0022", "B_5_Ind012_0023",
                                                   "B_5_Ind012_0044", "B_5_Ind012_0047", "B_5_Ind012_0051",
                                                   "B_5_Ind012_0061", "B_5_Ind012_0096", "B_5_Ind012_0098",
                                                   "B_5_Ind012_0100", "B_6_Ind013_0005", "B_6_Ind013_0014",
                                                   "B_6_Ind013_0016", "B_6_Ind013_0020",
                                                   "B_6_Ind013_0029", "B_6_Ind013_0031", "B_6_Ind013_0065",
                                                   "B_6_Ind013_0070", "B_6_Ind013_0085", "B_6_Ind013_0097"))
)
rownames(meirmans_polyploid) <- meirmans_polyploid$Individual
polygid <- df2genind(
    meirmans_polyploid[1:2],
    ploidy = 6,
    ncode = 3
  )
strata(polygid) <- meirmans_polyploid["pop"]
splitStrata(polygid) <- ~group/population/sample/genome
}

context("Published AMOVA results")

test_that("Amova returns published values", {

  skip_on_cran()
	expect_equivalent(res$componentsofcovariance[, 2], resper)
	expect_equivalent(res$componentsofcovariance[, 1], ressig)

	expect_equivalent(rescc$componentsofcovariance[, 2], resccper)
	expect_equivalent(rescc$componentsofcovariance[, 1], resccsig)
	expect_output(print(res), "ade4::amova")

})

test_that("pegas implemenation returns published values", {
  
  skip_on_cran()
  suppressWarnings({
  pres    <- poppr.amova(Aeut, ~Pop/Subpop, quiet = TRUE, method = "pegas")
  prescc  <- poppr.amova(Aeut, ~Pop/Subpop, quiet = TRUE, clonecorrect = TRUE, 
                         method = "pegas")
  })
  expect_equivalent(pres$varcomp, ressig[-4])
  expect_equivalent(pres$varcomp/sum(pres$varcomp), resper[-4]/100)
  expect_output(print(pres), "Variance components:")
  
  expect_equivalent(prescc$varcomp, resccsig[-4])
  expect_equivalent(prescc$varcomp/sum(prescc$varcomp), resccper[-4]/100)
})

context("AMOVA subsetting")

test_that("AMOVA handles subsetted genclone objects", {
	Athena.mlg   <- mlg.vector(Athena)
	agc.mlg      <- mlg.vector(agc)
	Athena.AMOVA <- poppr.amova(Athena, ~Subpop, quiet = TRUE)

	# All MLGs are represented
	expect_false(anyNA(match(1:max(agc.mlg), agc.mlg)))

	# Some MLGS are missing
	expect_true(anyNA(match(1:max(Athena.mlg), Athena.mlg)))

	#AMOVA will still work on the data set with missing MLGs
	expect_equal(class(Athena.AMOVA), "amova")
	})


test_that("AMOVA can take a distance matrix as input", {
  skip_on_cran()
  adist       <- diss.dist(clonecorrect(Aeut, strata = NA), percent = FALSE)
  root_adist  <- sqrt(adist)
  full_dist   <- dist(Aeut)
  
  # UNSQUARED
  usqres <- poppr.amova(Aeut, ~Pop/Subpop, quiet = TRUE, dist = root_adist, 
                        squared = FALSE)
  # SQUARED
  sqres <- poppr.amova(Aeut, ~Pop/Subpop, quiet = TRUE, dist = adist, 
                       squared = TRUE)
  
  # FULL provided distance will be trimmed
  fres <- poppr.amova(Aeut, ~Pop/Subpop, quiet = TRUE, dist = full_dist, 
                       squared = FALSE)
  expect_equivalent(usqres, res)
  expect_equivalent(usqres, sqres)
  expect_equivalent(sqres, fres)
  err <- c("Uncorrected.+?97.+?Clone.+?70.+?provided.+?119")
  expect_error(poppr.amova(Athena, ~Subpop, quiet = TRUE, dist = adist), err)
})

context("AMOVA correction and filtering")

test_that("AMOVA can work on clone correction from mlg.filter and do filtering", {
  skip_on_cran()
  data("monpop", package = "poppr")
  splitStrata(monpop) <- ~Tree/Year/Symptom
  expect_warning(poppr.amova(monpop, ~Symptom/Year), "Zero")
  mondist <- dist(monpop)
  THRESHOLD <- 1.75
  monfilt <- mlg.filter(monpop, distance = mondist, threshold = THRESHOLD, 
                        stats = "DIST")
  mlg.filter(monpop, distance = mondist) <- THRESHOLD
  monfiltdist <- as.dist(monfilt)
  monres <- poppr.amova(monpop, ~Symptom/Year, dist = monfiltdist, squared = FALSE)
  msg <- "Original.+?264"
  monpop@mlg <- mll(monpop, "original")
  expect_message(res <- poppr.amova(monpop, ~Symptom/Year, filter = TRUE, threshold = THRESHOLD), msg)
  expect_equivalent(monres$componentsofcovariance, res$componentsofcovariance)
})

test_that("AMOVA can calculate within individual variance for diploids", {
  skip_on_cran()
  data("microbov", package = "adegenet")
  strata(microbov) <- data.frame(other(microbov))
  mics <- microbov[pop = 1:2]
  expect_message(micwithin  <- poppr.amova(mics, ~breed), "Removing")
  expect_output(poppr.amova(mics, ~breed, correction = "cai"), "Cailliez constant")
  expect_output(poppr.amova(mics, ~breed, correction = "lin"), "Lingoes constant")
  expect_error(poppr.amova(mics, ~breed, correction = "wha", quiet = TRUE), "cailliez")
  expect_message(micwithout <- poppr.amova(mics, ~breed, within = FALSE), "Removing")
  expect_equal(dim(micwithin$componentsofcovariance), c(4, 2))
  expect_equal(dim(micwithout$componentsofcovariance), c(3, 2))
})

context("Polyploid AMOVA tests")


test_that("AMOVA can accurately calculate within-individual variance", {
  polyPhi  <- c(0.300757099877674, 0.0189408581368382, 0.127819751662941, 0.182803253215195)
  polyPerc <- c(69.9242900122326, 1.34999614286817, 10.4453885233797, 18.2803253215195)
  polyin   <- poppr.amova(polygid, ~group/population, within = TRUE)
  
  # Phi statistics are equal
  expect_equal(polyin$statphi$Phi, polyPhi)
  # Variance percentages are equal
  expect_equal(polyin$componentsofcovariance[4:1, "%", drop = TRUE], polyPerc)
})


test_that("AMOVA can accurately calculate rho", {
  rho      <- c(0.688374795330904, 0.445443116797669, 0.438064490572017)
  polyPerc <- c(31.1625204669096, 25.0310304758887, 43.8064490572017)
  polyno   <- poppr.amova(polygid, ~group/population, within = FALSE)
  
  # Phi statistics are equal
  expect_equal(polyno$statphi$Phi, rho)
  # Variance percentages are equal
  expect_equal(polyno$componentsofcovariance[3:1, "%", drop = TRUE], polyPerc)
})

test_that("Rho works with frequencies and counts for full data", {
  rho      <- c(0.688374795330904, 0.445443116797669, 0.438064490572017)
  polyPerc <- c(31.1625204669096, 25.0310304758887, 43.8064490572017)
  polyno   <- poppr.amova(polygid, ~group/population, within = FALSE, freq = FALSE)
  
  # Phi statistics are equal
  expect_equal(polyno$statphi$Phi, rho)
  # Variance percentages are equal
  expect_equal(polyno$componentsofcovariance[3:1, "%", drop = TRUE], polyPerc)
})

pg <- polygid
pg@tab[] <- as.integer(tab(pg) > 0)

test_that("AMOVA will not calculate within-individual variance for mixed ploids", {
  wrn      <- "ambiguous allele dosage"
  rho      <- c(0.688374795330904, 0.445443116797669, 0.438064490572017)
  expect_warning(res <- poppr.amova(pg, ~group/population, within = TRUE), wrn)
  expect_lt(sum(abs(res$statphi$Phi - rho)), 0.2)
})

test_that("AMOVA will give an extra warning for polyploids with zeroes", {
  skip_on_cran()
  wrn <- "zeroes encoded"
  expect_warning(res <- poppr.amova(recode_polyploids(pg, addzero = TRUE), ~group/population, within = TRUE), wrn)
})

context("AMOVA for vcfR-exported snp genind data")

m <- matrix(nrow = 10, ncol = 20)
set.seed(99)
strat     <- data.frame(POP = rep(c("A", "B"), 5))
m[1:5, ]  <- sample(c("00", "01", "11"), 5*20, replace = TRUE, prob = c(0.6, 0.3, 0.1))
m[6:10, ] <- sample(c("00", "01", "11"), 5*20, replace = TRUE, prob = c(0.2, 0.1, 0.7))
mm[mm == "00"] <- "AA"
mm[mm == "01"] <- "AC"
mm[mm == "11"] <- "CC"
x <- df2genind(m, strata = strat, ncode = 1)
y <- df2genind(mm, strata = strat, ncode = 1)

test_that("AMOVA will interpret the SNP data as having ambiguous allele dosage by default", {

  expect_warning({
    ax <- poppr.amova(x, ~POP, quiet = TRUE)
  }, "Data with mixed ploidy or ambiguous allele dosage")

  # We expect this to be equivalent to the true snps without asessing
  # within-sample differences
  ay <- poppr.amova(y, ~POP, quiet = TRUE, within = FALSE)
  expect_identical(ax, ay)
})


test_that("AMOVA will account for numeric SNP data if the user specifies so", {

  skip('2019-08-04: Still writing the functionality')
  
  op <- getOption('poppr.numeric.snp')
  options(poppr.numeric.snp = TRUE)

  ax <- poppr.amova(x, ~POP, quiet = TRUE)
  ay <- poppr.amova(y, ~POP, quiet = TRUE)
  expect_identical(ax, ay)

  options(poppr.numeric.snp = op)

})


context("AMOVA on genlight objects")

set.seed(99)
glite <- glSim(20, 10, 10, pop.freq = c(0.5, 0.5), ploidy = 2, parallel = FALSE, n.cores = 1L)
strata(glite) <- as.data.frame(other(glite))

test_that("AMOVA can be run on genlight objects", {
  skip_on_cran()
  # Both ade4 and pegas will run
  suppressWarnings({
    resa <- poppr.amova(glite, ~ancestral.pops, method = "ade4")
    resp <- poppr.amova(glite, ~ancestral.pops, method = "pegas")
  })
  # Both are (irritatingly) amova class
  expect_is(resa, "amova")
  expect_is(resp, "amova")
  # Both will have the same Mean Squared Error
  expect_equal(resp$tab$MSD, resa$results$`Mean Sq`)
})

test_that("AMOVA can work on filtered data", {
  skip_on_cran()
  noW <- poppr.amova(glite, ~ancestral.pops, within = FALSE)
  expect_warning(fil <- poppr.amova(glite, ~ancestral.pops, 
                                    filter = TRUE, 
                                    threshold = 0L, 
                                    quiet = TRUE), "within = FALSE")
  expect_message(fil2 <- poppr.amova(glite, ~ancestral.pops,
                                     filter = TRUE,
                                     within = FALSE,
                                     threshold = 2.25),
                 "Contracted multilocus genotypes ... 18")
  expect_identical(noW, fil)
  expect_failure(expect_identical(fil, fil2))
})
