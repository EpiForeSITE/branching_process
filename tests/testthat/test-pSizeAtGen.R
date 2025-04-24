test_that("pSizeAtGen() matches brute-force calculation", {
  pfs1 <- pSizeAtGen(g=1, n=1, j=4, R=1.2, k=0.33)
  pfs2 <- pSizeAtGen(g=2, n=1, j=4, R=1.2, k=0.33)
  pfs3 <- pSizeAtGen(g=3, n=1, j=4, R=1.2, k=0.33)
  
  d <- dnbinom(0:3, mu=1.2, size=0.33)
  dnb1 <- d[4]
  dnb2 <- d[3] * d[2] * d[1] * 2 + d[2] * d[3]
  dnb3 <- d[2]^3
  
  expect_equal(pfs1, dnb1, tolerance = sqrt(.Machine$double.eps))
  expect_equal(pfs2, dnb2, tolerance = sqrt(.Machine$double.eps))
  expect_equal(pfs3, dnb3, tolerance = sqrt(.Machine$double.eps))
})
