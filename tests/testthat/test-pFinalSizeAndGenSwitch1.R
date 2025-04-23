test_that("pFinalSizeAndGenSwitch1() matches brute-force calculation", {
  pfs1 <- pFinalSizeAndGenSwitch1(g=1, n=1, j=4, R0=1.2, k0=0.33, Rc=0.9, k=0.22)
  pfs2 <- pFinalSizeAndGenSwitch1(g=2, n=1, j=4, R0=1.2, k0=0.33, Rc=0.9, k=0.22)
  pfs3 <- pFinalSizeAndGenSwitch1(g=3, n=1, j=4, R0=1.2, k0=0.33, Rc=0.9, k=0.22)

  d0 <- dnbinom(0:3, mu=1.2, size=0.33)
  dc <- dnbinom(0:2, mu=0.9, size=0.22)
  dnb1 <- d0[4] * dc[1]^3
  dnb2 <- d0[3] * dc[2] * dc[1]^2 * 2 +
    d0[2] * dc[3] * dc[1]^2
  dnb3 <- d0[2] * dc[2]^2 * dc[1]
  
  expect_equal(pfs1, dnb1, tolerance = sqrt(.Machine$double.eps))
  expect_equal(pfs2, dnb2, tolerance = sqrt(.Machine$double.eps))
  expect_equal(pfs3, dnb3, tolerance = sqrt(.Machine$double.eps))
})
