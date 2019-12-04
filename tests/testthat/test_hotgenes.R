#test_hotgenes.R
# Author: Yi Fei Huang
#
context("hotgenes")

# ==== BEGIN SETUP AND PREPARE =================================================

#load(system.file("data", "musCh1fc.Rda", package = "hotgenes"))
#data(musCh1fc)

# ==== END SETUP AND PREPARE ===================================================

test_that("corrupt input generates errors",  {
  expect_error(generateStrandModel())
  expect_error(generateLocationModel())
  expect_error(rebuildStrandModel(NULL, newStartBase=3000000, newEndBase=3217000, newScaling=1000000,
                                  startBase=1, endBase=195471971, scaling=42))
  expect_error(plotHeatedMap())
})

test_that("a sample input produces the expected output",  {
  temp <- hotgenes::generateStrandModel(startBase=1, endBase=195471971, fc=musCh1fc, chr="chr1", strand="-", scaling=1000)
  temp1 <- hotgenes::rebuildStrandModel(strandModel=temp, newStartBase=3000000, newEndBase=3217000, newScaling=100000,
                            startBase=1, endBase=195471971, scaling=1000)
  expect_gt(max(temp), 0.1)
  expect_gt(max(temp1), 0.1)
  expect_lt(max(hotgenes::simulateHeatSpread(t(as.matrix(c(0,1))), 1, 1))-0.5, 0.1)
  expect_equal(nrow(temp1), nrow(generateLocationModel(3000000, 3217000, 100000)))
})

# ==== BEGIN TEARDOWN AND RESTORE ==============================================
# Remove every persistent construct that the test has created, except for
# stuff in tempdir().
#
#remove(musCh1fc)

# ==== END  TEARDOWN AND RESTORE ===============================================

# [END]
