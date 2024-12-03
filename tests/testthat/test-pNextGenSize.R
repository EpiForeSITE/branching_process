test_that("pNextGenSize() matches dnbinom", {
  expect_equal(pNextGenSize(x=1, y=3, R=2, k=0.5), dnbinom(3,mu=2,size=0.5), tolerance = sqrt(.Machine$double.eps))
})
