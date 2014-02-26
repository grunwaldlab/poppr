context("Poppr table tests")

test_that("poppr returns expected PA values", {
  data(Aeut, package = "poppr")
  A.tab <- poppr(Aeut, quiet = TRUE)
  comparison <- structure(list(Pop = structure(1:3, 
    .Label = c("Athena", "Mt. Vernon", "Total"), 
    class = "factor"), 
    N = c(97, 90, 187), 
    MLG = c(70, 50, 119), 
    eMLG = c(65.9808393377379, 50, 68.4525682280634), 
    SE = c(1.24567688244012, 0, 2.98857840353058), 
    H = c(4.06272002528149, 3.66843094399907, 4.55798828426928), 
    G = c(42.1928251121076, 28.7234042553191, 68.9723865877712), 
    Hexp = c(0.986469072164949, 0.976029962546817, 0.990799838997182), 
    E.5 = c(0.721008688944842, 0.725926650260449, 0.720112175857993), 
    Ia = c(2.90602921191748, 13.3024309367662, 14.3707995986407), 
    rbarD = c(0.0723700801886747, 0.281642324983496, 0.270617053778004), 
    File = structure(c(1L, 1L, 1L), class = "factor", .Label = "rootrot.csv")), 
    .Names = c("Pop", "N", "MLG", "eMLG", "SE", "H", "G", "Hexp", "E.5", "Ia", "rbarD", "File"), 
    row.names = c(NA, -3L), 
    class = c("popprtable", "data.frame"))

  expect_that(A.tab$Pop, is_equivalent_to(comparison$Pop))
  expect_that(A.tab$N, equals(comparison$N))
  expect_that(A.tab$MLG, equals(comparison$MLG))
  expect_that(A.tab$eMLG, equals(comparison$eMLG))
  expect_that(A.tab$SE, equals(comparison$SE))
  expect_that(A.tab$H, equals(comparison$H))
  expect_that(A.tab$G, equals(comparison$G))
  expect_that(A.tab$Hexp, equals(comparison$Hexp))
  expect_that(A.tab$E.5, equals(comparison$E.5))
  expect_that(A.tab$Ia, equals(comparison$Ia))
  expect_that(A.tab$rbarD, equals(comparison$rbarD))
})

test_that("poppr returns expected codominant values", {
  data(partial_clone, package = "poppr")
  p.tab <- poppr(partial_clone, quiet = TRUE)
  comparison <- structure(list(Pop = structure(1:5, .Label = c("1", "2", "3", 
"4", "Total"), class = "factor"), N = c(13, 13, 12, 12, 50), 
    MLG = c(10, 12, 11, 9, 26), eMLG = c(9.46153846153846, 11.1538461538462, 
    11, 9, 9.93701353621932), SE = c(0.498518515262143, 0.36080121229411, 
    0, 0, 1.12881059579593), H = c(2.24503527412618, 2.45831132968308, 
    2.36938211969468, 2.09472904752765, 3.07152395656842), G = c(8.89473684210526, 
    11.2666666666667, 10.2857142857143, 7.2, 17.8571428571429
    ), Hexp = c(0.961538461538461, 0.987179487179487, 0.984848484848485, 
    0.939393939393939, 0.963265306122449), E.5 = c(0.935312405238733, 
    0.960842907662783, 0.958200460752105, 0.870390481833875, 
    0.819311895784525), Ia = c(2.1580763424628, 1.87492360969648, 
    1.15572679509632, 1.157153633392, 1.93513179817012), rbarD = c(0.243225877705591, 
    0.212786561587854, 0.132460530412697, 0.13328661732193, 0.217470007471919
    ), File = structure(c(1L, 1L, 1L, 1L, 1L), class = "factor", .Label = "partial_clone.dat")), .Names = c("Pop", 
"N", "MLG", "eMLG", "SE", "H", "G", "Hexp", "E.5", "Ia", "rbarD", 
"File"), row.names = c(NA, -5L), class = c("popprtable", "data.frame"
))

  expect_that(p.tab$Pop, is_equivalent_to(comparison$Pop))
  expect_that(p.tab$N, equals(comparison$N))
  expect_that(p.tab$MLG, equals(comparison$MLG))
  expect_that(p.tab$eMLG, equals(comparison$eMLG))
  expect_that(p.tab$SE, equals(comparison$SE))
  expect_that(p.tab$H, equals(comparison$H))
  expect_that(p.tab$G, equals(comparison$G))
  expect_that(p.tab$Hexp, equals(comparison$Hexp))
  expect_that(p.tab$E.5, equals(comparison$E.5))
  expect_that(p.tab$Ia, equals(comparison$Ia))
  expect_that(p.tab$rbarD, equals(comparison$rbarD))
})