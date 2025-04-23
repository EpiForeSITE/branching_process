test_that("pFinalSize() matches brute-force calculation", {
  pfs <- pFinalSize(n=1, j=4, R=1.2, k=0.33)
  
  d <- dnbinom(0:3, mu=1.2, size=0.33)
  dnb <- d[4] * d[1]^3 +
    d[3] * d[2] * d[1]^2 * 2 +
    d[2] * d[3] * d[1]^2 +
    d[2]^3 * d[1]
  
  expect_equal(pfs, dnb, tolerance = sqrt(.Machine$double.eps))
})
