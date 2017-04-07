context("Analytical value tests")

snps <- c(1587L, 1451L, 910L, 1899L, 1474L, 986L, 539L, 44L, 1035L, 1054L, 
1396L, 183L, 1693L, 55L, 770L, 1410L, 1845L, 328L, 607L, 1559L, 
1637L, 530L, 722L, 966L, 244L, 28L, 533L, 1281L, 77L, 532L, 1805L, 
564L, 79L, 1507L, 1984L, 1020L, 1289L, 1658L, 460L, 1792L, 478L, 
924L, 1491L, 1386L, 975L, 1321L, 1393L, 1959L, 981L, 531L, 1760L, 
836L, 516L, 1268L, 370L, 1892L, 225L, 1624L, 1765L, 1144L, 860L, 
714L, 1234L, 492L, 1356L, 223L, 1996L, 609L, 1030L, 1930L, 709L, 
375L, 867L, 598L, 1987L, 1493L, 60L, 1463L, 98L, 1852L, 577L, 
168L, 1179L, 919L, 1677L, 899L, 138L, 686L, 927L, 473L, 1265L, 
1975L, 1201L, 921L, 1788L, 1182L, 1217L, 566L, 865L, 1004L, 1849L, 
989L, 1239L, 1002L, 1232L, 443L, 1271L, 604L, 1155L, 1319L, 422L, 
1800L, 556L, 1751L, 1567L, 1783L, 288L, 738L, 1840L, 1391L, 813L, 
1172L, 628L, 1378L, 1181L, 736L, 988L, 584L, 503L, 1242L, 1283L, 
1432L, 1741L, 1654L, 1101L, 998L, 1034L, 1490L, 590L, 1154L, 
1940L, 1L, 721L, 879L, 866L, 818L, 1186L, 1235L, 1350L, 578L, 
458L, 1258L, 970L, 1949L, 1300L, 959L, 450L, 1339L, 812L, 437L, 
1965L, 1146L, 944L, 1879L, 135L, 1371L, 56L, 1093L, 451L, 538L, 
1084L, 1448L, 1661L, 1672L, 1250L, 890L, 1422L, 685L, 1660L, 
1580L, 32L, 260L, 1261L, 1469L, 1614L, 493L, 362L, 316L, 1534L, 
1540L, 345L, 1538L, 1470L, 181L, 1369L, 434L, 296L, 291L, 1696L, 
1213L, 61L, 264L, 218L, 1088L, 23L, 717L, 776L, 352L, 546L, 599L, 
909L, 815L, 1500L, 1827L, 745L, 66L, 937L, 11L, 404L, 1238L, 
1280L, 549L, 1058L, 1704L, 1720L, 1698L, 1455L, 1535L, 41L, 1483L, 
287L, 180L, 956L, 337L, 497L, 1556L, 512L, 1767L, 824L, 1961L, 
554L, 33L, 85L, 339L, 713L, 1670L, 523L, 20L, 1665L, 432L, 1543L, 
1091L, 271L, 691L, 1332L, 1132L, 254L, 1357L, 1295L, 490L, 1240L, 
1963L, 374L, 38L, 1583L, 688L, 1440L, 863L, 1903L, 506L, 766L, 
632L, 1163L, 1531L, 101L, 878L, 82L, 1197L, 1013L, 1545L, 1502L, 
1752L, 958L, 786L, 1785L, 1428L, 1007L, 740L, 1511L, 954L, 1376L, 
281L, 1995L, 1377L, 36L, 1623L, 1969L, 286L, 1630L, 6L, 1033L, 
1573L, 471L, 1257L, 1367L, 859L, 614L, 480L, 1810L, 1218L, 1627L, 
528L, 561L, 1143L, 1916L, 1514L, 644L, 1173L, 1537L, 1098L, 1429L, 
49L, 1059L, 829L, 1781L, 543L, 505L, 911L, 429L, 1778L, 62L, 
1225L, 838L, 280L, 1641L, 907L, 1488L, 1499L, 1065L, 558L, 703L, 
359L, 1884L, 1921L, 129L, 994L, 400L, 1263L, 1216L, 1625L, 995L, 
1971L, 1411L, 125L, 1400L, 1078L, 1528L, 1249L, 118L, 667L, 1945L, 
886L, 1109L, 1900L, 552L, 948L, 1317L, 773L, 876L, 626L, 692L, 
1175L, 771L, 1170L, 664L, 1292L, 1691L, 332L, 190L, 1737L, 734L, 
672L, 1145L, 283L, 90L, 971L, 1418L, 1530L, 915L, 1669L, 892L, 
617L, 1140L, 1564L, 206L, 502L, 1365L, 1750L, 1209L, 527L, 263L, 
1923L, 1569L, 1403L, 1729L, 1318L, 600L, 819L, 239L, 1408L, 1160L, 
1640L, 1503L, 1977L, 785L, 1893L, 1690L, 754L, 1207L, 1320L, 
1826L, 1496L, 681L, 1842L, 1725L, 1354L, 1210L, 1557L, 1551L, 
1464L, 1610L, 1956L, 586L, 1306L, 1044L, 220L, 1438L, 896L, 1643L, 
1565L, 1520L, 1814L, 1843L, 491L, 279L, 849L, 1073L, 760L, 1359L, 
1113L, 992L, 1412L, 1632L, 1379L, 1868L, 247L, 1711L, 1860L, 
1085L, 874L, 1885L, 378L, 660L, 1703L, 1424L, 42L, 1913L, 583L, 
1326L, 1313L, 485L, 925L, 796L, 834L, 1026L, 272L, 1402L, 290L, 
633L, 1738L, 1414L, 1120L, 1070L, 1248L, 1323L, 515L, 96L, 187L, 
1364L, 1069L, 624L, 1329L, 1125L, 1363L, 1264L, 514L, 1187L, 
1766L, 1023L, 916L)

test_that("Bruvo's distance works as expected.", {
  testdf  <- data.frame(test = c("00/20/23/24", "20/24/26/43"))
  testgid <- df2genind(testdf, ploidy = 4, sep = "/")
  addloss <- as.vector(bruvo.dist(testgid, add = FALSE, loss = FALSE))
  ADDloss <- as.vector(bruvo.dist(testgid, add = TRUE, loss = FALSE))
  addLOSS <- as.vector(bruvo.dist(testgid, add = FALSE, loss = TRUE))
  ADDLOSS <- as.vector(bruvo.dist(testgid, add = TRUE, loss = TRUE))
  # Values from Bruvo et. al. (2004)
  expect_equal(addloss, 0.46875000000000)
  expect_equal(addLOSS, 0.34374987334013)
  expect_equal(ADDloss, 0.458333164453506)
  expect_equal(ADDLOSS, 0.401041518896818)
})

test_that("Bruvo's distance will trim extra zeroes.", {
  testdf  <- data.frame(test = c("00/20/24/26/43", "00/00/20/23/24"))
  testgid <- df2genind(testdf, ploidy = 5, sep = "/")
  addloss <- as.vector(bruvo.dist(testgid, add = FALSE, loss = FALSE))
  ADDloss <- as.vector(bruvo.dist(testgid, add = TRUE, loss = FALSE))
  addLOSS <- as.vector(bruvo.dist(testgid, add = FALSE, loss = TRUE))
  ADDLOSS <- as.vector(bruvo.dist(testgid, add = TRUE, loss = TRUE))
  # Values from Bruvo et. al. (2004)
  expect_equal(addloss, 0.46875000000000)
  expect_equal(addLOSS, 0.34374987334013)
  expect_equal(ADDloss, 0.458333164453506)
  expect_equal(ADDLOSS, 0.401041518896818)
})

test_that("Bruvo's distance will go through the recusion", {
  skip_on_cran()
  testdf  <- data.frame(test = c("00/00/00/00/24", "00/20/24/26/43"))
  testgid <- df2genind(testdf, ploidy = 5, sep = "/")
  ADDloss <- as.vector(bruvo.dist(testgid, add = TRUE, loss = FALSE))
  addloss <- as.vector(bruvo.dist(testgid, add = FALSE, loss = FALSE))
  addLOSS <- as.vector(bruvo.dist(testgid, add = FALSE, loss = TRUE))
  ADDLOSS <- as.vector(bruvo.dist(testgid, add = TRUE, loss = TRUE))
  expect_equal(ADDloss, 0.671874523162842)
  expect_equal(addLOSS, 0.351561918854713)
  expect_equal(addloss, 0.75)
  expect_equal(ADDLOSS, 0.511718221008778)
})

test_that("Repeat lengths can be in any order and length if named", {
  skip_on_cran()
  data("Pram")
  p10 <- Pram[sample(nInd(Pram), 10)]
  # Repeat length as normal
  pbruvo <- bruvo.dist(p10, replen = other(p10)$REPLEN)
  # Repeat lengths mixed up
  pbruvo_sampled <- bruvo.dist(p10, replen = sample(other(p10)$REPLEN))
  # Extra value in repeat lengths
  pbruvo_long <- bruvo.dist(p10, replen = c(other(p10)$REPLEN, NOPE = 23))
  expect_equivalent(pbruvo, pbruvo_sampled)
  expect_equivalent(pbruvo, pbruvo_long)
})

test_that("Bruvo's distance can be calculated per locus", {
  skip_on_cran()
  data("Pram")
  p10 <- Pram[sample(nInd(Pram), 10)]
  pbruvo <- bruvo.dist(p10, replen = other(p10)$REPLEN, by_locus = TRUE)
  expect_is(pbruvo, "list")
  expect_is(pbruvo[[1]], "dist")
  expect_equal(length(pbruvo), nLoc(p10))
})

test_that("Infinite Alleles Model works.",{
  x <- structure(list(V3 = c("228/236/242", "000/211/226"), 
                      V6 = c("190/210/214", "000/190/203")), 
                 .Names = c("V3", "V6"), row.names = c("1", "2"), 
                 class = "data.frame")
  gid <- df2genind(x, sep = "/", ploidy = 3)
  res <- bruvo.dist(gid, replen = c(2/10000,2/10000), add = FALSE, loss = FALSE)
  expect_equal(as.numeric(res), 0.833333333333333)
})

test_that("Dissimilarity distance works as expected.", {
  data(nancycats, package = "adegenet")
  nan1 <- popsub(nancycats, 1)
  nanmat <- diss.dist(nan1, mat = TRUE)
  expect_that(diss.dist(nan1), is_a("dist"))
  expect_that(nanmat, is_a("matrix"))
  expect_equal(diss.dist(nan1, mat = TRUE, percent = TRUE), (nanmat/2)/9)
  expect_equal(nanmat[2, 1], 4)
})

test_that("Index of association works as expected.", {
  data(Aeut, package = "poppr")
  # Values from Grünwald and Hoheisel (2006)
  res <- c(Ia = 14.3707995986407, rbarD = 0.270617053778004)
  expect_equal(ia(Aeut), res)
})

test_that("Internal function fix_negative_branch works as expected.", {
  the_tree <- structure(list(edge = structure(c(9L, 10L, 11L, 12L, 13L, 14L,
      14L, 13L, 12L, 11L, 10L, 9L, 9L, 10L, 11L, 12L, 13L, 14L, 2L, 3L, 6L, 7L,
      8L, 1L, 5L, 4L), .Dim = c(13L, 2L)), 
      edge.length = c(0, 0.0625, 0.0625, 0.09375, 0.15, -0.25, 0.25, -0.15,
      -0.09375, -0.0625, 0, 0, 0),
      tip.label = c("2340_50156.ab1 ", "2340_50149.ab1 ", "2340_50674.ab1 ",
      "2370_45312.ab1 ", "2340_50406.ab1 ", "2370_45424.ab1 ", "2370_45311.ab1
      ", "2370_45521.ab1 "), 
      Nnode = 6L
    ),
    .Names = c("edge", "edge.length", "tip.label", "Nnode"), 
    class = "phylo",
    order = "cladewise")
  fix_tree <- poppr:::fix_negative_branch(the_tree)
  # Not all branch lengths are positive
  expect_false(min(the_tree$edge.length) >= 0)
  # After fix, all branch lengths are positive
  expect_true(min(fix_tree$edge.length) >= 0)
  # The difference from fixed and unfixed is unfixed. This indicates that the
  # clones were set to zero and the fix set the branch lengths in the correct 
  # order.
  expect_equivalent(min(fix_tree$edge.length - the_tree$edge.length), min(the_tree$edge.length))
})

test_that("mlg.matrix returns a matrix and not table", {
  data(partial_clone)
  mat4row <- poppr:::mlg.matrix(partial_clone)
  pop(partial_clone) <- NULL
  mat1row <- poppr:::mlg.matrix(partial_clone)
  expect_that(mat4row, is_a("matrix"))
  expect_that(mat1row, is_a("matrix"))
  expect_false(inherits(mat1row, "table"))
  expect_false(inherits(mat4row, "table"))
})

test_that("diversity_stats returns expected values", {
  skip_on_cran()
  data("Aeut", package = "poppr")
  expected <- structure(c(4.06272002528149, 3.66843094399907, 4.55798828426928,
  42.1928251121076, 28.7234042553191, 68.9723865877712, 0.976299287915825,
  0.965185185185185, 0.985501444136235, 0.721008688944842, 0.725926650260449,
  0.720112175857993), .Dim = 3:4, .Dimnames = structure(list(Pop = c("Athena",
  "Mt. Vernon", "Total"), Index = c("H", "G", "simp", "E.5")), .Names = c("Pop",
  "Index")))
  
  atab       <- mlg.table(Aeut, plot = FALSE, total = TRUE)
  res        <- diversity_stats(atab)
  pop(Aeut)  <- NULL
  res_single <- diversity_stats(atab["Total", , drop = FALSE])
  
  expect_equivalent(res, expected)
  expect_false(is.matrix(res_single))
  expect_equivalent(res_single, res["Total", ])
})

test_that("get_boot_x works with one pop or one stat", {
  skip_on_cran()
  data(Pinf)
  Ptab <- mlg.table(Pinf, plot = FALSE)
  pop(Pinf) <- NULL
  ptab <- mlg.table(Pinf, plot = FALSE)
  
  Pboot_all <- diversity_ci(Ptab, 20L, plot = FALSE)
  Pboot_E   <- diversity_ci(Ptab, 20L, G = FALSE, H = FALSE, lambda = FALSE, plot = FALSE)
  pboot_all <- diversity_ci(ptab, 20L, plot = FALSE)
  pboot_E   <- diversity_ci(ptab, 20L, G = FALSE, H = FALSE, lambda = FALSE, plot = FALSE)
  
  Past <- poppr:::get_boot_stats(Pboot_all$boot)
  PEst <- poppr:::get_boot_stats(Pboot_E$boot)
  past <- poppr:::get_boot_stats(pboot_all$boot)
  pEst <- poppr:::get_boot_stats(pboot_E$boot)
  
  
  expect_equivalent(dim(Past), c(2, 4))
  expect_equivalent(dim(PEst), c(2, 1))
  expect_equivalent(dim(past), c(1, 4))
  expect_equivalent(dim(pEst), c(1, 1))
  
})

test_that("ia returns NA with less than three samples", {
  skip_on_cran()
  data(partial_clone)
  pc <- partial_clone[sample(50, 2)]
  expect_equivalent(rep(NA_real_, 2), ia(pc))
  expect_equivalent(rep(NA_real_, 4), ia(pc, sample = 9))
})

test_that("ia and pair.ia return same values", {
  skip_on_cran()
  data(partial_clone)
  pc_pair <- pair.ia(partial_clone, plot = FALSE, quiet = TRUE)
  expect_output(print(pc_pair), "Locus_1:Locus_2")
  
  # Randomly sample two loci
  set.seed(9001)
  locpair <- sample(locNames(partial_clone), 2)
  
  # Calculate ia for those
  pc_ia   <- ia(partial_clone[loc = locpair])
  
  # Find the pair in the table
  pair_posi <- grepl(locpair[1], rownames(pc_pair)) & grepl(locpair[2], rownames(pc_pair))

  expect_equivalent(pc_pair[pair_posi], pc_ia)
})


test_that("bitwise.ia can handle large samples", {
  skip_on_cran()
  set.seed(999)
  x <- glSim(n.ind = 200, n.snp.nonstruc = 2e3, ploidy = 2, parallel = FALSE)
  position(x) <- sort(sample(1e4, 2e3))
  res <- bitwise.ia(x[, snps], thread = 1L)
  expect_equal(res, 8.6296328853274e-06)
})

test_that("win.ia produces expected results", {
  skip_on_cran()
  set.seed(999)
  x <- glSim(n.ind = 10, n.snp.nonstruc = 5e2, n.snp.struc = 5e2, ploidy = 2)
  position(x) <- sort(sample(1e4, 1e3))
  pos.res <- c(0.0455840455840456, -0.121653107452204, -0.00595292101094314,
               0.0080791417110234, -0.0168361601753453, -0.00529454073927957,
               -0.00897633972053186, -0.0127370947405975, -0.0449102196309011,
               -0.00973722044247706, -0.00395345516884293, -0.00365271629340326,
               -0.00974233012233703, -0.00245476433207628, 0.00284396024347222,
               -0.00961958529835491, 0.0131920863329584, 0.153481745744672,
               0.267842562306642, 0.241714085779086, 0.218062029786597,
               0.160460027459485, 0.351036788453414, 0.165485011482322,
               0.275776270167356, 0.292113177590906, 0.201538950412372,
               0.146277858717017, 0.171243142808176, 0.213214941672605,
               0.200245399275377, 0.243007310007378, 0.12396191060666,
               0.0350644169627312 )

  nopos.res <- c(0.00124186533451839, 0.0334854151956646, 0.226632983481153, 
                 0.173916930910343)

  win.pos     <- win.ia(x, window = 300L, quiet = TRUE)
  position(x) <- NULL
  win.nopos   <- win.ia(x, window = 300L, quiet = TRUE)
  expect_equivalent(win.pos, pos.res)
  expect_equivalent(win.nopos, nopos.res)
  
  # Making sure that the genind version does the same thing.
  xmat <- as.matrix(x)
  xmat[xmat == 0]   <- "1/1"
  xmat[xmat == "1"] <- "1/2"
  xmat[xmat == "2"] <- "2/2"
  xgid <- df2genind(xmat[, 1:300], sep = "/", ind.names = .genlab("", 10), 
                    loc.names = .genlab("", 300L))
  gid.res <- ia(xgid)["rbarD"]
  expect_equivalent(gid.res, nopos.res[1])
})

test_that("win.ia can handle missing data for haploids", {
  skip_on_cran()
  set.seed(999)
  dat <- sample(c(0, 1, NA), 500, replace = TRUE, prob = c(0.47, 0.47, 0.06))
  mat <- matrix(dat, nrow = 5, ncol = 100)
  x   <- new("genlight", mat, parallel = FALSE)
  gid <- df2genind(mat + 1, ind.names = .genlab("", 5), ploidy = 1,
                   loc.names = .genlab("", 100))
  bit.res <- win.ia(x, quiet = TRUE)
  gid.res <- ia(gid)["rbarD"]
  expect_equivalent(bit.res, gid.res)
})

test_that("win.ia can handle missing data for diploids", {
  skip_on_cran()
  set.seed(999)
  dat <- sample(c(0, 1, 2, NA), 500, replace = TRUE, prob = c(0.22, 0.5, 0.22, 0.06))
  mat <- matrix(dat, nrow = 5, ncol = 100)
  x   <- new("genlight", mat, parallel = FALSE)
  mat[mat == 0]   <- "1/1"
  mat[mat == "1"] <- "1/2"
  mat[mat == "2"] <- "2/2"
  gid <- df2genind(mat, ind.names = .genlab("", 5), ploidy = 2, sep = "/",
                   loc.names = .genlab("", 100))
  bit.res <- win.ia(x, quiet = TRUE)
  gid.res <- ia(gid)["rbarD"]
  expect_equivalent(bit.res, gid.res)
})

test_that("win.ia knows the batman theme", {
  skip_on_cran()
  set.seed(999)
  x <- glSim(n.ind = 10, n.snp.nonstruc = 2.5e2, n.snp.struc = 2.5e2, ploidy = 2)
  position(x) <- sort(sample(5e3, 5e2))
  res <- win.ia(x, window = 325L, min.snps = 100L, quiet = TRUE)
  expect_true(all(is.na(res)))
  # NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA
  # BATMAN!
})

test_that("samp.ia works",{
  skip_on_cran()
  set.seed(999)
  x <- glSim(n.ind = 10, n.snp.nonstruc = 5e2, n.snp.struc = 5e2, ploidy = 2)
  position(x) <- sort(sample(1e4, 1e3))
  set.seed(900)
  pos.res <- samp.ia(x, n.snp = 20, reps = 2, quiet = TRUE)
  expect_equivalent(pos.res, c(0.0754028380744677, 0.0178846504949451))

  # Sampling is agnostic to position
  position(x) <- NULL
  set.seed(900)
  nopos.res <- samp.ia(x, n.snp = 20, reps = 2, quiet = TRUE)
  expect_equivalent(nopos.res, pos.res)
})

test_that("fix_replen works as expected", {
  data(nancycats)
  nanrep  <- rep(2, 9)
  nantest <- test_replen(nancycats, nanrep)
  nanfix  <- fix_replen(nancycats, nanrep)
  expect_equal(sum(nantest), 5)
  expect_true(all(floor(nanfix)[!nantest] == 1))
})

test_that("fix_replen throws errors for weird replens", {
  skip_on_cran()
  data(partial_clone)
  expect_warning(fix_replen(partial_clone, rep(10, 10)), paste(locNames(partial_clone), collapse = ", "))
  expect_warning(fix_replen(partial_clone, rep(10, 10), fix_some = FALSE), "Original repeat lengths are being returned")
  expect_warning(fix_replen(partial_clone, rep(2, 10)), "The repeat lengths for Locus_2, Locus_7, Locus_9 are not consistent.")
  expect_warning(fix_replen(partial_clone, rep(2, 10)), "Repeat lengths with some modification are being returned: Locus_3")
})

test_that("poppr_has_parallel returns something logical", {
  expect_is(poppr_has_parallel(), "logical")
})