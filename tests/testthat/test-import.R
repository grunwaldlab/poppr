context("Data import tests")

data(monpop, package = "poppr")
data(Pinf, package = "poppr")
pr <- recode_polyploids(Pinf, newploidy = TRUE)

test_that("basic text connections work", {
	y <- "13	6	1	6											
				7_09_BB											
	Ind	Pop	CHMFc4	CHMFc5	CHMFc12	SEA	SED	SEE	SEG	SEI	SEL	SEN	SEP	SEQ	SER
	A004	7_09_BB	224	85	163	132	133	156	144	116	143	227	257	142	145
	A002	7_09_BB	224	97	159	156	129	156	144	113	143	231	261	136	153
	A011	7_09_BB	224	97	159	160	133	156	126	119	147	227	257	134	149
	A009	7_09_BB	224	97	159	160	133	156	126	119	147	227	261	134	149
	A006	7_09_BB	224	97	159	160	133	156	126	119	147	235	261	134	149
	A013	7_09_BB	224	97	163	160	133	156	126	119	147	235	257	134	149"
	
	gen <- read.genalex(textConnection(y), sep = "\t")
	expect_equivalent(gen@tab, monpop[1:6, drop = TRUE]@tab)
})

test_that("genclone objects can be saved and restored", {
	mp <- file()
	genind2genalex(monpop, filename = mp, quiet = TRUE)
	gen <- read.genalex(mp)
	close(mp)
	
	expect_equal(gen@tab, monpop@tab)
})

test_that("polyploids can be saved", {
	skip_on_cran()
	file1 <- tempfile()
	file2 <- tempfile()
	genind2genalex(Pinf, filename = file1, quiet = TRUE)
	genind2genalex(pr, filename = file2, quiet = TRUE)
	Pinf2 <- read.genalex(file1, ploidy = 4)
	pr2   <- read.genalex(file2, ploidy = 3)
	expect_equal(summary(Pinf, verbose = FALSE)$He, summary(Pinf2, verbose = FALSE)$He)
	expect_equal(summary(pr2, verbose = FALSE)$NA.perc, summary(Pinf, verbose = FALSE)$NA.perc)
	expect_true(all(ploidy(pr2) == 3))
	expect_true(all(ploidy(Pinf2) == 4))
})

test_that("errors are reported", {
	skip_on_cran()
	file1 <- tempfile()
	file2 <- tempfile()
	genind2genalex(Pinf, filename = file1, quiet = TRUE)
	genind2genalex(pr, filename = file2, quiet = TRUE)
	expect_that(Pinf2 <- read.genalex(file1, ploidy = 4), not(throws_error()))

	expect_error(Pinf2 <- read.genalex(file1), "set the flag?")
	expect_error(Pinf2 <- read.genalex(file1, geo = TRUE), "geo = TRUE")
	expect_error(Pinf2 <- read.genalex(file1, region = TRUE), "region = TRUE")
})


test_that("genalex data can be imported with region data to genind and genclone", {
	skip_on_cran()
	rr1 <- system.file("files/rootrot.csv", package = "poppr")
	rr2 <- system.file("files/rootrot2.csv", package = "poppr")

	root1gc <- read.genalex(rr1)
	root1gd <- read.genalex(rr1, genclone = FALSE)
	
	root2gc <- read.genalex(rr2)
	root2gd <- read.genalex(rr2, genclone = FALSE)
	
	root2re <- read.genalex(rr2, region = TRUE)

	expect_is(root1gc, "genclone")
	expect_is(root1gd, "genind")

	expect_is(root2gc, "genclone")
	expect_is(root2gd, "genind")
	expect_equal(length(nameStrata(root2gc)), 1L)
	expect_equal(length(nameStrata(root2re)), 2L)
	expect_identical(nameStrata(root2re), c("Pop", "Region"))
})

test_that("genalex can import geographic information", {
	skip_on_cran()
	data("Pram", package = "poppr")
	filepram <- tempfile()
	sourpram <- tempfile()
	staypram <- tempfile()
	custpram <- tempfile()
	custpop  <- sample(.genlab("p", 10), nInd(Pram), replace = TRUE)

	expect_output(genind2genalex(Pram, filename = filepram, geo = TRUE), filepram)
	expect_output(genind2genalex(Pram, pop = ~SOURCE, allstrata = FALSE, filename = sourpram, geo = TRUE), sourpram)
	expect_output(genind2genalex(Pram, pop = ~STATE/YEAR, allstrata = FALSE, filename = staypram, geo = TRUE), staypram)
	expect_output(genind2genalex(Pram, pop = custpop, allstrata = FALSE, filename = custpram, geo = TRUE), custpram)

	expect_error(read.genalex(filepram))
	pall <- read.genalex(filepram, geo = TRUE)
	sour <- read.genalex(sourpram, geo = TRUE)
	stay <- read.genalex(staypram, geo = TRUE)
	cust <- read.genalex(custpram, geo = TRUE)

	splitStrata(pall) <- ~SOURCE/YEAR/STATE
	nameStrata(sour)  <- ~SOURCE
	splitStrata(stay) <- ~STATE/YEAR

	expect_equal(nInd(pall), nInd(Pram))
	expect_equal(nInd(sour), nInd(Pram))
	expect_equal(nInd(stay), nInd(Pram))
	expect_equal(nInd(cust), nInd(Pram))

	expect_equivalent(strata(pall), strata(Pram))
	expect_equivalent(strata(sour), strata(Pram, ~SOURCE, combine = FALSE))
	expect_equivalent(strata(stay), strata(Pram, ~STATE/YEAR, combine = FALSE))

	expect_equivalent(other(pall)$xy[1:513, ], other(Pram)$xy)
	expect_equivalent(other(sour)$xy[1:513, ], other(Pram)$xy)
	expect_equivalent(other(stay)$xy[1:513, ], other(Pram)$xy)
	expect_equivalent(other(cust)$xy[1:513, ], other(Pram)$xy)

})