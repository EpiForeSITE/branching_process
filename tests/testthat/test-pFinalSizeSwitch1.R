test_that("pFinalSizeSwitch1() matches brute-force calculation: n=1, j=4", {
  pfs <- pFinalSizeSwitch1(n=1, j=4, R0=1.2, k0=0.33, Rc=0.9, kc=1.1)
  
  d0 <- dnbinom(0:3, mu=1.2, size=0.33)
  dc <- dnbinom(0:2, mu=0.9, size=1.1)
  dnb <- d0[4] * dc[1]^3 +
    d0[3] * dc[2] * dc[1]^2 * 2 +
    d0[2] * dc[3] * dc[1]^2 +
    d0[2] * dc[2]^2 * dc[1]
  
  expect_equal(pfs, dnb, tolerance = sqrt(.Machine$double.eps))
})

test_that("pFinalSizeSwitch1() matches brute-force calculation: n=1, j=5", {
  pfs <- pFinalSizeSwitch1(n=1, j=5, R0=1.2, k0=0.33, Rc=0.9, kc=1.1)
  
  d0 <- dnbinom(0:4, mu=1.2, size=0.33)
  dc <- dnbinom(0:3, mu=0.9, size=1.1)
  dnb <- d0[5] * dc[1]^4 +
    d0[4] * dc[2] * dc[1]^3 * 3 +
    d0[3] * dc[3] * dc[1]^3 * 2 +
    d0[3] * dc[2]^2 * dc[1]^2 +
    d0[3] * dc[2]^2 * dc[1]^2 * 2 +
    d0[2] * dc[4] * dc[1]^3 +
    d0[2] * dc[3] * dc[2] * dc[1]^2 * 2 +
    d0[2] * dc[2] * dc[3] * dc[1]^2 +
    d0[2] * dc[2]^3 * dc[1]
  
  expect_equal(pfs, dnb, tolerance = sqrt(.Machine$double.eps))
})
