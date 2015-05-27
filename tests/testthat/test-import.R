context("Data import tests")

data(monpop, package = "poppr")
data(Pinf, package = "poppr")
pr <- recode_polyploids(Pinf, newploidy = TRUE)

mysummary <- function(...){
	tmp <- tempfile()
	sink(tmp)
	x <- summary(...)
	sink()
	return(x)
}

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
	expect_that(mysummary(Pinf)$He, equals(mysummary(Pinf2)$He))
	expect_that(mysummary(pr2)$NA.perc, equals(mysummary(Pinf)$NA.perc))
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
