#test_hotgenes.R
# Author: Yi Fei Huang
#
context("hotgenes")

# ==== BEGIN SETUP AND PREPARE =================================================

load(system.file("data", "musCh1fc.Rda", package = "hotgenes"))

# ==== END SETUP AND PREPARE ===================================================

test_that("corrupt input generates errors",  {
  expect_error(generateStrandModel())
  expect_error(generateLocationModel())
  expect_error(plotHeatedMap())
})

test_that("a sample input produces the expected output",  {
  expect_gt(max(generateStrandModel(startBase=1, endBase=195471971, fc=musCh1fc, chr="chr1", strand="-", scaling=100000)), 0.1)
  expect_lt(max(simulateHeatSpread(t(as.matrix(c(0,1))), 1, 1))-0.5, 0.1)
})


# ==== BEGIN TEARDOWN AND RESTORE ==============================================
# Remove every persistent construct that the test has created, except for
# stuff in tempdir().
#

remove(musCh1fc)

# ==== END  TEARDOWN AND RESTORE ===============================================

# [END]
