test_that("pNextGenSize() with one transmitter matches dnbinom", {
  expect_equal(pNextGenSize(x=1, y=3, R=2, k=0.5), dnbinom(3,mu=2,size=0.5), tolerance = sqrt(.Machine$double.eps))
})

test_that("pNextGenSize() with multiple transmitters matches brute-force dnbinom", {
  pngs <- pNextGenSize(x=3, y=2, R=1.1, k=0.44)
  
  d <- dnbinom(0:2, mu=1.1, size=0.44)
  dnb <- 3 * d[1]^2 * d[3] + 3 * d[1] * d[2]^2
  
  expect_equal(pngs, dnb, tolerance = sqrt(.Machine$double.eps))
})
