context("Data import tests")

data(monpop, package = "poppr")
data(Pinf, package = "poppr")
data(H3N2, package = "adegenet")
pr <- recode_polyploids(Pinf, newploidy = TRUE)
y <- "13	6	1	6											
			7_09_BB											
Ind	Pop	CHMFc4	CHMFc5	CHMFc12	SEA	SED	SEE	SEG	SEI	SEL	SEN	SEP	SEQ	SER
A004	7_09_BB	224	85	163	132	133	156	144	116	143	227	257	142	145
A002	7_09_BB	224	97	159	156	129	156	144	113	143	231	261	136	153
A011	7_09_BB	224	97	159	160	133	156	126	119	147	227	257	134	149
A009	7_09_BB	224	97	159	160	133	156	126	119	147	227	261	134	149
A006	7_09_BB	224	97	159	160	133	156	126	119	147	235	261	134	149
A013	7_09_BB	224	97	163	160	133	156	126	119	147	235	257	134	149"

yd <- "13	6	1	6											
			7_09_BB											
Ind	Pop	CHMFc4	CHMFc5	CHMFc12	SEA	SED	SEE	SEG	SEI	SEL	SEN	SEP	SEQ	SER
4	7_09_BB	224	85	163	132	133	156	144	116	143	227	257	142	145
2	7_09_BB	224	97	159	156	129	156	144	113	143	231	261	136	153
2	7_09_BB	224	97	159	160	133	156	126	119	147	227	257	134	149
9	7_09_BB	224	97	159	160	133	156	126	119	147	227	261	134	149
6	7_09_BB	224	97	159	160	133	156	126	119	147	235	261	134	149
3	7_09_BB	224	97	163	160	133	156	126	119	147	235	257	134	149"

zz <- "1	6	1	6
7_09_BB			
Ind	Pop	CHMFc4	CHMFc5
A004	7_09_BB	224	85
A002	7_09_BB	224	97
A011	7_09_BB	224	97
A009	7_09_BB	224	97
A006	7_09_BB	224	97
A013	7_09_BB	224	97"


zzna <- "13	6	1	6											
			7_09_BB											
Ind	Pop	CHMFc4	CHMFc5	CHMFc12	SEA	SED	SEE	SEG	SEI	SEL	SEN	SEP	SEQ	SER
A004	7_09_BB	224	85	163	132	133	156	 	116	143	227	257	142	145
A002	7_09_BB	224	97	159	156	129	156	144	113	143	231	261	136	153
A011	7_09_BB	224	97	159	160	133	156	126	119	147	227	257	134	149
A009	7_09_BB	224	97	159	160	133	156	126	119	147	227	261	134	149
A006	7_09_BB	224	97	159	160	133	156	126	119	147	235	261	134	149
A013	7_09_BB	224	97	163	160	133	156	126	119	147	235	257	134	149"

z <- "1	6	1	6
7_09_BB			
Ind	Pop	CHMFc4	
A004	7_09_BB	224	 
A002	7_09_BB	224	 
A011	7_09_BB	224	 
A009	7_09_BB	224	 
A006	7_09_BB	224	 
A013	7_09_BB	224	 "

zna <- "3	6	1	6											
			7_09_BB											
Ind	Pop	CHMFc4	CHMFc5	CHMFc12	SEA	SED	SEE	SEG	SEI	SEL	SEN	SEP	SEQ
A002	7_09_BB	224	97	159	156	129	156	144	113	143	231	261	136
A004	7_09_BB	224	85	163	0	133	156	144		143	227	257	0
A011	7_09_BB	224	97	159	160	133	156	126	119	147	227	257	134
A009	7_09_BB	224	97	159	160	133	156	126	119	147	227	261	134
A006	7_09_BB	224	97	159	160	133	156	126	119	147	235	261	134
A013	7_09_BB	224	97	163	160	133	156	126	119	147	235	257	134"

hapdip <- "6	4	1	4											
			7_09_BB											
Ind	Pop	CHMFc4	CHMFc5	CHMFc12	SEA	SED	SEE	SEG	SEI	SEL	SEN	SEP	SEQ
A011	7_09_BB	224	0	159	0	133	0	126	0	147		257	0
A009	7_09_BB	224	97	159	160	133	156	126	119	147	227	261	134
A006	7_09_BB	224	97	159	160	133	156	126	119	147	235	261	134
A013	7_09_BB	224	97	163	160	133	156	126	119	147	235	257	134"

bad_genalex <- "1	6	1	6
Ind	Pop	CHMFc4	
A004	7_09_BB	224	 
A002	7_09_BB	224	 
A011	7_09_BB	224	 
A009	7_09_BB	224	 
A006	7_09_BB	224	 
A013	7_09_BB	224	 "

missing_single <-  "13	6	1	6											
			7_09_BB											
Ind	Pop	CHMFc4	CHMFc5	CHMFc12	SEA	SED	SEE	SEG	SEI	SEL	SEN	SEP	SEQ	SER
A004	7_09_BB_A004	224	85	163	132	133	156	144	116	143	227	257	142	145
A002	7_09_BB_A002	0	0	0	0	0	0	0	0	0	0	0	0	0
A011	7_09_BB_A011	224	97	159	160	133	156	126	119	147	227	257	134	149
A009	7_09_BB_A009	224	97	159	160	133	156	126	119	147	227	261	134	149
A006	7_09_BB_A006	224	97	159	160	133	156	126	119	147	235	261	134	149
A013	7_09_BB_A013	224	97	163	160	133	156	126	119	147	235	257	134	149"

test_that("basic text connections work", {
	gen <- read.genalex(textConnection(y), sep = "\t")
	expect_equivalent(tab(gen), tab(monpop[1:6, drop = TRUE]))
})

test_that("names are corrected properly", {
	expect_warning(gen <- read.genalex(textConnection(yd), sep = "\t"),
                 "duplicate labels detected")
  expect_false(anyNA(strata(gen)))
  expect_named(other(gen), "original_names")
  expect_identical(indNames(gen), as.character(1:6))
  indNames(gen) <- sprintf("A%03d", c(4, 2, 11, 9, 6, 13))
	expect_equivalent(tab(gen), tab(monpop[1:6, drop = TRUE]))
})

test_that("missing rows and columns are eliminated", {
  gen <- read.genalex(textConnection(zzna), sep = "\t")
  expect_true(any(is.na(tab(gen))))
  expect_equal(nInd(gen), 6L)
  expect_equal(nLoc(gen), 13L)
})

test_that("single locus diploids can be imported", {
  gen <- read.genalex(textConnection(zz), sep = "\t")
  expect_equivalent(nLoc(gen), 1L)
  expect_output(show(gen), "diploid")
})

test_that("single locus haploids can be imported", {
  gen <- read.genalex(textConnection(z), sep = "\t")
  expect_equivalent(nLoc(gen), 1L)
  expect_output(show(gen), "haploid")
})

test_that("missing cells are converted to zeroes for polyploids", {
  skip_on_cran()
  gen <- read.genalex(textConnection(zna), sep = "\t", ploidy = 4L)
  expect_equivalent(nLoc(gen), 3L)
  expect_output(show(gen), "tetraploid")
  expect_output(show(recode_polyploids(gen, newploidy = TRUE)), "triploid \\(1\\) and tetraploid \\(5\\)")
})

test_that("haplodiploids can be imported correctly", {
  skip_on_cran()
  gen <- read.genalex(textConnection(hapdip), sep = "\t")
  expect_equivalent(nLoc(gen), 6L)
  expect_output(show(gen), "diploid")
  expect_output(show(recode_polyploids(gen, newploidy = TRUE)), "haploid \\(1\\) and diploid \\(3\\)")
})

test_that("duplicate columns are flagged and fixed", {
  skip_on_cran()
  f <- "4,5,1,5,,,,,,
,,,Admix,,,,,,
Ind,Pop,RM127, ,RM22, ,RM22, ,RM127, 
1,Admix,210,210,200,200,195,195,130,110
2,Admix,230,230,185,185,200,200,110,120
3,Admix,210,210,200,200,195,195,130,130
4,Admix,230,230,200,200,195,195,130,130
5,Admix,210,230,200,200,200,200,120,120"
  expect_warning(read.genalex(textConnection(f)), "col 7: RM22 -> RM22_1")
})

test_that("missing samples do not shift strata", {
  skip_on_cran()
  expect_warning(ms <- read.genalex(textConnection(missing_single), sep = "\t"), 
    "[Ii]ndividual[s(][ s][)]?(deleted|with no scored loci have been removed)")
  expect_equal(as.character(strata(ms)$Pop), as.character(pop(ms)))
  expect_equal(rownames(strata(ms)), indNames(ms))
})

test_that("missing samples do not shift strata, even with duplicated names", {
  skip_on_cran()
  missing_single2 <- gsub("A004\t", "A011\t", missing_single)
  expect_warning(ms <- read.genalex(textConnection(missing_single2), sep = "\t"), "duplicate labels detected")
  expect_equal(as.character(strata(ms)$Pop), as.character(pop(ms)))
  expect_equal(rownames(strata(ms)), indNames(ms))
})

test_that("improperly-formatted data causes an error", {
  skip_on_cran()
  msg <- "^.+?6 individuals.+?5 rows.+?Please inspect "
  tcmsg  <- paste0(msg, "textConnection\\(bad_genalex\\).+?$")
  expect_error(read.genalex(textConnection(bad_genalex), sep = "\t"), tcmsg)
  f <- tempfile()
  writeLines(bad_genalex, f)
  fmsg <- paste0(msg, f, ".+?$")
  expect_error(read.genalex(f, sep = "\t"), fmsg)
})

test_that("sample names with apostrophes can be imported", {
  skip_on_cran()
  better_than_yar <- "1,5,1,5
,,,7_09_BB,
Ind,Pop,CHMFc4,CHMFc5
phaser,7_09_BB,224,85
rock,7_09_BB,224,97
bat'leth,7_09_BB,224,97
paper,7_09_BB,224,97
scissors,7_09_BB,224,97"
  
  res <- poppr::read.genalex(textConnection(better_than_yar))
  expect_is(res, "genclone")
  expect_equal(nInd(res), 5L)
  expect_equal(nLoc(res), 1L)
  expect_equal(indNames(res), c("phaser", "rock", "bat'leth", "paper", "scissors"))
})

test_that("loci with entirely T loci are not converted to TRUE", {
  # https://github.com/grunwaldlab/poppr/issues/214
  tea <- read.genalex(test_path("genalex", "test.txt"))
  expected <- list(
    `605-4471` = c("T", "C"), 
    `681-4471` = c("G", "T", "T"), 
    `682-4471` = c("G", "T", "T")
  )
  expect_setequal(alleles(tea), expected)

})

context("Data export tests")

test_that("not specifying a file for genind2genalex will generate a tempfile", {
  skip_on_cran()
  expect_warning(f <- genind2genalex(monpop, quiet = TRUE), "temporary file")
  expect_match(f, "^/.+?file.+\\.csv$")
  expect_is(read.genalex(f), "genclone")
})

test_that("genind2genalex will prevent a file from being overwritten", {
  skip_on_cran()
  f <- tempfile()
  writeLines("hey!\n", f)
  expect_error(genind2genalex(monpop, filename = f, quiet = TRUE), "exists and will not be overwritten")
  expect_match(readLines(f)[1], "hey")
})

test_that("genclone objects can be saved and restored", {
	mp <- file()
	genind2genalex(monpop, filename = mp, quiet = TRUE)
	gen <- read.genalex(mp)
	close(mp)
	
	expect_equal(gen@tab, monpop@tab)
})

test_that("genalex will give a warning if user asks for geo data when there is none", {
  skip_on_cran()
  mp <- file()
  expect_warning(genind2genalex(monpop, filename = mp, quiet = TRUE, geo = TRUE), 
                 "monpop@other")
  close(mp)
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

test_that("diploid missing data is handled correctly", {
  skip_on_cran()
  data("nancycats", package = "adegenet")
  file1 <- tempfile()
  genind2genalex(nancycats, file1, quiet = TRUE)
  nan <- read.genalex(file1, genclone = FALSE)
  expect_identical(summary(nancycats, verbose = FALSE), summary(nan, verbose = FALSE))
})

test_that("sequence data is handled correctly", {
  skip_on_cran()
  tmp <- tempfile()
  htab <- tab(H3N2[1:10, loc = 1:10, drop = TRUE])
  genind2genalex(H3N2[1:10, loc = 1:10], filename = tmp, quiet = TRUE, overwrite = TRUE)
  h3n2 <- read.genalex(tmp)
  # The alleles are imported in a different order, so I have to resort with the
  # column names.
  expect_equivalent(htab, tab(h3n2)[, colnames(htab)])
  genind2genalex(H3N2[1:10, loc = 1:10], filename = tmp, sequence = TRUE, quiet = TRUE, overwrite = TRUE)
  h3n2numbers <- read.genalex(tmp)
  # The sequence option converts letters to numbers. If the last test worked,
  # then this test should work, too.
  expect_equivalent(tab(h3n2), tab(h3n2numbers))
})

test_that("errors are reported", {
	skip_on_cran()
	file1 <- tempfile()
	file2 <- tempfile()
	genind2genalex(Pinf, filename = file1, quiet = TRUE)
	genind2genalex(pr, filename = file2, quiet = TRUE)
	expect_error(Pinf2 <- read.genalex(file1, ploidy = 4), NA)

	expect_error(Pinf2 <- read.genalex(file1), "set the flag?")
	expect_error(Pinf2 <- read.genalex(file1, geo = TRUE), "geo = TRUE")
	expect_error(Pinf2 <- read.genalex(file1, region = TRUE), "region = TRUE")
})

context("Extra info data import tests")

test_that("genalex data can be imported with region data to genind and genclone", {
	skip_on_cran()
	rr1 <- system.file("files/rootrot.csv", package = "poppr")
	rr2 <- system.file("files/rootrot2.csv", package = "poppr")
	Xcoord <- rnorm(187)
	Ycoord <- rnorm(187)
	rrg <- read.csv(rr2, header = FALSE) 
	blank <- rep("", nrow(rrg))
	rrg <- cbind(rrg, blank, data.frame(X = c(NA, NA, "x", Xcoord), 
	                                    Y = c(NA, NA, "y", Ycoord)))
	rrfile <- tempfile()
	write.table(rrg, file = rrfile, quote = FALSE, sep = ",", row.names = FALSE,
	            col.names = FALSE, na = "")
	root1gc <- read.genalex(rr1)
	root1gd <- read.genalex(rr1, genclone = FALSE)
	
	root2gc <- read.genalex(rr2)
	root2gd <- read.genalex(rr2, genclone = FALSE)
	
	root2re <- read.genalex(rr2, region = TRUE)

	root2reg <- read.genalex(rrfile, region = TRUE, geo = TRUE)
	
	expect_is(root1gc, "genclone")
	expect_is(root1gd, "genind")

	expect_is(root2gc, "genclone")
	expect_is(root2gd, "genind")

	expect_is(root2reg, "genclone")
	
	expect_equal(length(nameStrata(root2gc)), 1L)
	expect_equal(length(nameStrata(root2re)), 2L)
	expect_identical(nameStrata(root2re), c("Pop", "Region"))
	expect_identical(nameStrata(root2reg), c("Pop", "Region"))
	expect_equivalent(other(root2reg)$xy, 
	                  data.frame(x = Xcoord, y = Ycoord))
})

test_that("genalex data can be imported with a region column", {
  skip_on_cran()
  yr <- read.table(textConnection(y), sep = "\t", header = FALSE, row.names = NULL)
  yr <- as.matrix(yr)
  yrnums <- yr[1, 1:4]
  region <- c("", "", "Region", rep(c("one", "two"), 3))
  blank <- rep("", 9)
  Xcoords <- rnorm(6)
  X       <- c("", "", "x", Xcoords)
  Ycoords <- rnorm(6)
  Y       <- c("", "", "y", Ycoords)
  yrg         <- yr
  yr          <- cbind(yr[, 2], region, yr[, -c(1:2)], blank, yr[, 1])
  yr[1, ]     <- c(yrnums, rep("", ncol(yr) - 4))
  yr[1, 5:7]  <- c("2", "3", "3")
  yr[2, 6:7]  <- c("one", "two")
  yrg         <- cbind(yr, X, Y)
  yrfile  <- tempfile()
  yrgfile <- tempfile()
  write.table(yr, file = yrfile, quote = FALSE, sep = ",", row.names = FALSE,
              col.names = FALSE)
  write.table(yrg, file = yrgfile, quote = FALSE, sep = ",", row.names = FALSE,
              col.names = FALSE)
  
  genind_region     <- read.genalex(yrfile, region = TRUE)
  genind_region_geo <- read.genalex(yrgfile, region = TRUE, geo = TRUE)
  
  expect_is(genind_region, "genclone")
  expect_is(genind_region_geo, "genclone")
  
  expect_equal(nameStrata(genind_region), c("Pop", "Region"))
  expect_equal(nameStrata(genind_region_geo), c("Pop", "Region"))
  
  setPop(genind_region) <- ~Region
  setPop(genind_region_geo) <- ~Region
  
  expect_equal(popNames(genind_region), c("one", "two"))
  expect_equal(popNames(genind_region_geo), c("one", "two"))
  
  expect_equivalent(other(genind_region_geo)$xy, 
                    data.frame(x = Xcoords, y = Ycoords))
})


test_that("genalex can import geographic information", {
	skip_on_cran()
	data("Pram", package = "poppr")
	filepram <- tempfile()
	sourpram <- tempfile()
	staypram <- tempfile()
	custpram <- tempfile()
	custpop  <- sample(.genlab("p", 10), nInd(Pram), replace = TRUE)

	expect_output(show(genind2genalex(Pram, filename = filepram, geo = TRUE)), filepram)
	expect_output(show(genind2genalex(Pram, pop = ~SOURCE, allstrata = FALSE, filename = sourpram, geo = TRUE)), sourpram)
	expect_output(show(genind2genalex(Pram, pop = ~STATE/YEAR, allstrata = FALSE, filename = staypram, geo = TRUE)), staypram)
	expect_output(show(genind2genalex(Pram, pop = custpop, allstrata = FALSE, filename = custpram, geo = TRUE)), custpram)

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
