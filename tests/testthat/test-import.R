context("Data import tests")

test_that("basic text connections work", {
	data(monpop, package = "poppr")
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
	expect_equal(gen@tab, monpop[1:6, drop = TRUE]@tab)
})

test_that("genclone objects can be saved and restored", {
	mp <- file()
	genind2genalex(monpop, filename = mp, quiet = TRUE)
	gen <- read.genalex(mp)
	close(mp)
	
	expect_equal(gen@tab, monpop@tab)
})