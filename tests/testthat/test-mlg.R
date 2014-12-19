context("Multilocus genotype tests")

test_that("multilocus genotype vector is same length as samples", {
  data(Aeut, package = "poppr")
  data(partial_clone, package = "poppr")
  data(nancycats, package = "adegenet")
  amlg <- mlg.vector(Aeut)
  pmlg <- mlg.vector(partial_clone)
  nmlg <- mlg.vector(nancycats)
  expect_that(length(amlg), equals(nInd(Aeut)))
  expect_that(length(pmlg), equals(nInd(partial_clone)))
  expect_that(length(nmlg), equals(nInd(nancycats)))
  expect_that(length(unique(amlg)), equals(mlg(Aeut, quiet = TRUE)))
  expect_that(length(unique(pmlg)), equals(mlg(partial_clone, quiet = TRUE)))
  expect_that(length(unique(nmlg)), equals(mlg(nancycats, quiet = TRUE)))
})

test_that("multilocus genotype matrix matches mlg.vector and data", {
  data(Aeut, package = "poppr")
  data(partial_clone, package = "poppr")
  data(nancycats, package = "adegenet")
  aclone <- as.genclone(Aeut)
  atab   <- mlg.table(Aeut, bar = FALSE)
  ptab   <- mlg.table(partial_clone, bar = FALSE)
  ntab   <- mlg.table(nancycats, bar = FALSE)
  expect_that(nrow(atab), equals(length(Aeut@pop.names)))
  expect_that(nrow(ptab), equals(length(partial_clone@pop.names)))
  expect_that(nrow(ntab), equals(length(nancycats@pop.names)))
  expect_that(ncol(atab), equals(mlg(Aeut, quiet = TRUE)))
  expect_that(ncol(ptab), equals(mlg(partial_clone, quiet = TRUE)))
  expect_that(ncol(ntab), equals(mlg(nancycats, quiet = TRUE)))
  expect_that(sum(atab), equals(nInd(Aeut)))
  expect_that(sum(ptab), equals(nInd(partial_clone)))
  expect_that(sum(ntab), equals(nInd(nancycats)))
})

test_that("mlg.crosspop will work with subsetted genclone objects", {
  data(Aeut, package = "poppr")
  agc <- as.genclone(Aeut)
  Athena <- popsub(agc, "Athena")
  setpop(Athena) <- ~Subpop
  expected_output <- structure(list(MLG.13 = structure(c(1L, 1L), .Names = c("8", 
"9")), MLG.23 = structure(c(1L, 1L), .Names = c("4", "6")), MLG.24 = structure(c(1L, 
1L), .Names = c("9", "10")), MLG.32 = structure(c(1L, 1L), .Names = c("7", 
"9")), MLG.52 = structure(c(1L, 1L), .Names = c("5", "9")), MLG.63 = structure(c(1L, 
1L), .Names = c("1", "5"))), .Names = c("MLG.13", "MLG.23", "MLG.24", 
"MLG.32", "MLG.52", "MLG.63"))
  expected_mlgout <- c(13, 23, 24, 32, 52, 63)

  expect_that(x <- mlg.crosspop(Athena, quiet = TRUE), equals(expected_output))
  expect_that(y <- mlg.crosspop(Athena, indexreturn = TRUE), equals(expected_mlgout))
  expect_warning(z <- mlg.crosspop(Athena, mlgsub = c(14, 2:5)), "The following multilocus genotypes are not defined in this dataset: 2, 3, 4, 5")
})

test_that("mlg.id Aeut works", {
  data(Aeut, package = "poppr")
  expected_output <- structure(list(`1` = "55", `2` = c("101", "103"), `3` = "111", 
                                    `4` = "112", `5` = "110", `6` = "102", `7` = "20", `8` = "7", 
                                    `9` = "68", `10` = "69", `11` = "73", `12` = "75", `13` = c("72", 
                                                                                                "80"), `14` = c("74", "76", "77"), `15` = "79", `16` = c("4", 
                                                                                                                                                         "9"), `17` = c("3", "8"), `18` = "95", `19` = "94", `20` = c("22", 
                                                                                                                                                                                                                      "23", "24", "25", "27", "28", "29", "30", "31"), `21` = "60", 
                                    `22` = "43", `23` = c("38", "59"), `24` = c("84", "90"), 
                                    `25` = "63", `26` = "5", `27` = "71", `28` = "32", `29` = "78", 
                                    `30` = "26", `31` = c("89", "92"), `32` = c("65", "81"), 
                                    `33` = "53", `34` = "51", `35` = c("46", "48", "50"), `36` = c("45", 
                                                                                                   "47"), `37` = "88", `38` = "87", `39` = "56", `40` = "91", 
                                    `41` = "82", `42` = "6", `43` = "83", `44` = "13", `45` = "17", 
                                    `46` = "85", `47` = "61", `48` = "62", `49` = "66", `50` = "64", 
                                    `51` = "15", `52` = c("52", "86"), `53` = "2", `54` = "115", 
                                    `55` = "151", `56` = "113", `57` = "42", `58` = "109", `59` = c("159", 
                                                                                                    "57"), `60` = c("67", "70"), `61` = "58", `62` = "49", `63` = c("1", 
                                                                                                                                                                    "54"), `64` = "96", `65` = "40", `66` = c("33", "34", "36", 
                                                                                                                                                                                                              "39", "41"), `67` = "37", `68` = "35", `69` = c("145", "146", 
                                                                                                                                                                                                                                                              "148", "149"), `70` = c("124", "126", "127", "131", "133"
                                                                                                                                                                                                                                                              ), `71` = "156", `72` = c("152", "154"), `73` = "116", `74` = c("139", 
                                                                                                                                                                                                                                                                                                                              "140", "141"), `75` = c("134", "135", "137", "142", "147"
                                                                                                                                                                                                                                                                                                                              ), `76` = c("125", "162"), `77` = c("160", "168", "170"), 
                                    `78` = c("169", "177"), `79` = "175", `80` = c("107", "108", 
                                                                                   "117", "120", "121", "122", "164", "167", "172", "183"), 
                                    `81` = c("130", "182"), `82` = "99", `83` = "100", `84` = "114", 
                                    `85` = "157", `86` = "98", `87` = c("158", "171"), `88` = c("123", 
                                                                                                "166"), `89` = "118", `90` = c("128", "163"), `91` = c("104", 
                                                                                                                                                       "173"), `92` = "132", `93` = "10", `94` = "11", `95` = "180", 
                                    `96` = c("138", "144"), `97` = c("181", "184", "185", "186"
                                    ), `98` = "143", `99` = c("136", "165"), `100` = "150", `101` = c("174", 
                                                                                                      "187"), `102` = "176", `103` = c("178", "179"), `104` = "129", 
                                    `105` = "153", `106` = "119", `107` = "161", `108` = "97", 
                                    `109` = "93", `110` = "18", `111` = "21", `112` = "12", `113` = "16", 
                                    `114` = "19", `115` = "155", `116` = "106", `117` = "105", 
                                    `118` = "14", `119` = "44"), .Names = c("1", "2", "3", "4", 
                                                                            "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", 
                                                                            "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", 
                                                                            "27", "28", "29", "30", "31", "32", "33", "34", "35", "36", "37", 
                                                                            "38", "39", "40", "41", "42", "43", "44", "45", "46", "47", "48", 
                                                                            "49", "50", "51", "52", "53", "54", "55", "56", "57", "58", "59", 
                                                                            "60", "61", "62", "63", "64", "65", "66", "67", "68", "69", "70", 
                                                                            "71", "72", "73", "74", "75", "76", "77", "78", "79", "80", "81", 
                                                                            "82", "83", "84", "85", "86", "87", "88", "89", "90", "91", "92", 
                                                                            "93", "94", "95", "96", "97", "98", "99", "100", "101", "102", 
                                                                            "103", "104", "105", "106", "107", "108", "109", "110", "111", 
                                                                            "112", "113", "114", "115", "116", "117", "118", "119"))
  x    <- mlg.id(Aeut)
  Avec <- mlg.vector(Aeut)
  expect_that(x, equals(expected_output))
  expect_that(length(x), equals(length(unique(Avec))))
  expect_that(sapply(x, length), is_equivalent_to(as.vector(table(Avec))))
  expect_that(names(x[1]), equals("1"))
  })

test_that("mlg.id Pinf works", {
  data(Pinf, package = "poppr")
  expected_output <- structure(list(`1` = "PiEC06", `4` = "PiMX03", `5` = "PiMX04", 
                                    `6` = "PiMXT01", `7` = "PiPE03", `8` = "PiPE01", `10` = "PiPE07", 
                                    `11` = "PiPE06", `12` = c("PiPE10", "PiPE26"), `13` = "PiMX01", 
                                    `14` = "PiEC02", `15` = "PiPE04", `17` = c("PiCO01", "PiCO03", 
                                                                               "PiCO04"), `19` = "PiEC03", `21` = "PiMX42", `22` = "PiEC01", 
                                    `23` = "PiPE09", `24` = "PiCO02", `25` = "PiPE05", `30` = "PiMX07", 
                                    `33` = "PiMX20", `34` = c("PiMX48", "PiMX49", "PiMX50"), 
                                    `35` = "PiPE13", `36` = c("PiPE11", "PiPE12", "PiPE14"), 
                                    `37` = "PiMX06", `38` = "PiMX02", `39` = "PiMX12", `40` = "PiMXT06", 
                                    `41` = "PiMX19", `42` = "PiMX17", `45` = "PiMX13", `46` = "PiMX24", 
                                    `47` = c("PiPE02", "PiPE08"), `50` = "PiMX23", `51` = "PiMX10", 
                                    `52` = "PiMX29", `53` = "PiMX05", `54` = "PiCO05", `55` = "PiMXT07", 
                                    `56` = "PiMX11", `57` = "PiMX26", `58` = "PiMX22", `59` = "PiMX14", 
                                    `61` = "PiMX18", `62` = "PiMX15", `63` = c("PiPE22", "PiPE24", 
                                                                               "PiPE25"), `68` = "PiPE23", `69` = "PiEC10", `71` = "PiPE21", 
                                    `72` = "PiPE20", `74` = "PiEC12", `75` = c("PiEC13", "PiEC14"
                                    ), `77` = "PiMX28", `79` = "PiEC11", `80` = "PiMX16", `83` = "PiEC08", 
                                    `84` = "PiEC07", `93` = "PiMX30", `94` = "PiMX41", `95` = "PiMX27", 
                                    `96` = "PiMX43", `97` = c("PiMX44", "PiMX45", "PiMX46", "PiMX47"
                                    ), `98` = "PiMX25", `99` = "PiMX40", `104` = "PiMXT02", `105` = "PiMXT05", 
                                    `106` = "PiPE27", `109` = "PiMXT03", `110` = "PiMX21", `115` = "PiMXT04", 
                                    `116` = "PiMXt48", `117` = "PiMXt68"), .Names = c("1", "4", 
                                                                                      "5", "6", "7", "8", "10", "11", "12", "13", "14", "15", "17", 
                                                                                      "19", "21", "22", "23", "24", "25", "30", "33", "34", "35", "36", 
                                                                                      "37", "38", "39", "40", "41", "42", "45", "46", "47", "50", "51", 
                                                                                      "52", "53", "54", "55", "56", "57", "58", "59", "61", "62", "63", 
                                                                                      "68", "69", "71", "72", "74", "75", "77", "79", "80", "83", "84", 
                                                                                      "93", "94", "95", "96", "97", "98", "99", "104", "105", "106", 
                                                                                      "109", "110", "115", "116", "117"))
  x    <- mlg.id(Pinf)
  Pvec <- mlg.vector(Pinf)
  expect_that(x, equals(expected_output))
  expect_that(length(x), equals(length(unique(Pvec))))
  expect_that(sapply(x, length), is_equivalent_to(as.vector(table(Pvec))))
  expect_that(names(x[1]), equals("1"))
})

test_that("multilocus genotype filtering functions correctly", {
  data(Aeut, package = "poppr")
  data(partial_clone, package = "poppr")
  data(nancycats, package = "adegenet")
  amlg <- mlg.vector(Aeut)
  pmlg <- mlg.vector(partial_clone)
  nmlg <- mlg.vector(nancycats)
 
  # No clustering should happen if the threshold is set to 0
  expect_that(length(unique(amlg)), equals(length(unique(mlg.filter(Aeut,0,distance="diss.dist")@mll))))
  expect_that(length(unique(pmlg)), equals(length(unique(mlg.filter(partial_clone,0,distance="diss.dist")@mll))))
  expect_that(length(unique(nmlg)), equals(length(unique(mlg.filter(nancycats,0,distance="diss.dist")@mll))))
  # All clusters should be merged for an arbitrarily large threshold
  expect_that(1, equals(length(unique(mlg.filter(Aeut,10000,distance="diss.dist")@mll))))
  expect_that(1, equals(length(unique(mlg.filter(partial_clone,10000,distance="diss.dist")@mll))))
  expect_that(1, equals(length(unique(mlg.filter(nancycats,10000,distance="diss.dist")@mll))))
  # The different methods of passing distance should produce the same results
  adis <- diss.dist(missingno(Aeut,"mean",quiet=TRUE))
  pdis <- diss.dist(missingno(partial_clone,"mean",quiet=TRUE))
  ndis <- diss.dist(missingno(nancycats,"mean",quiet=TRUE))
  expect_that(mlg.filter(Aeut,0.3,missing="mean",distance=adis)@mll, equals(mlg.filter(Aeut,0.3,missing="mean",distance=diss.dist)@mll))
  expect_that(mlg.filter(Aeut,0.3,missing="mean",distance=adis)@mll, equals(mlg.filter(Aeut,0.3,missing="mean",distance="diss.dist")@mll))
  expect_that(mlg.filter(nancycats,0.3,missing="mean",distance=ndis)@mll, equals(mlg.filter(nancycats,0.3,missing="mean",distance=diss.dist)@mll))
  expect_that(mlg.filter(nancycats,0.3,missing="mean",distance=ndis)@mll, equals(mlg.filter(nancycats,0.3,missing="mean",distance="diss.dist")@mll))
  expect_that(mlg.filter(partial_clone,0.3,missing="mean",distance=pdis)@mll, equals(mlg.filter(partial_clone,0.3,missing="mean",distance=diss.dist)@mll))
  expect_that(mlg.filter(partial_clone,0.3,missing="mean",distance=pdis)@mll, equals(mlg.filter(partial_clone,0.3,missing="mean",distance="diss.dist")@mll))
})
