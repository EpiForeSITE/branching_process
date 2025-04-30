test_that("pFinalSizeAndGenSwitch2() matches brute-force calculation", {
  pfs1 <- pFinalSizeAndGenSwitch2(g=1, n=1, j=5, R0=1.2, k0=0.33, Rc=0.9, k=0.22)
  pfs2 <- pFinalSizeAndGenSwitch2(g=2, n=1, j=5, R0=1.2, k0=0.33, Rc=0.9, k=0.22)
  pfs3 <- pFinalSizeAndGenSwitch2(g=3, n=1, j=5, R0=1.2, k0=0.33, Rc=0.9, k=0.22)
  pfs4 <- pFinalSizeAndGenSwitch2(g=4, n=1, j=5, R0=1.2, k0=0.33, Rc=0.9, k=0.22)
  
  d0 <- dnbinom(0:4, mu=1.2, size=0.33)
  dc <- dnbinom(0:2, mu=0.9, size=0.22)
  dnb1 <- d0[5] * d0[1]^4
  dnb2 <- d0[4] * d0[2] * d0[1]^2 * dc[1] * 3 +
    d0[3] * d0[3] * d0[1] * dc[1]^2 * 2 +
    d0[3] * d0[2]^2 * dc[1]^2 +
    d0[2] * d0[4] * dc[1]^3
  dnb3 <- d0[3] * d0[2] * d0[1] * dc[2] * dc[1] * 2 +
    d0[2] * d0[3] * dc[2] * dc[1]^2 * 2 + 
    d0[2]^2 * dc[3] * dc[1]^2
  dnb4 <- d0[2]^2 * dc[2]^2 * dc[1]
  
  expect_equal(pfs1, dnb1, tolerance = sqrt(.Machine$double.eps))
  expect_equal(pfs2, dnb2, tolerance = sqrt(.Machine$double.eps))
  expect_equal(pfs3, dnb3, tolerance = sqrt(.Machine$double.eps))
  expect_equal(pfs4, dnb4, tolerance = sqrt(.Machine$double.eps))
})
