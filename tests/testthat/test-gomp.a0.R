test_that("gomp.a0() produces correct output with default values", {
  gomp.a0_output <- gomp.a0()
  expect_identical(gomp.a0_output, c(-66.76844784, -2.32502545, 0.0823))
  })

test_that("gomp.a0() produces correct output with non-default values", {
  gomp.a0_output <- withr::with_seed(
    seed = 123,
    gomp.a0(sampling = 1000,
                            b_min = 0.02,
                            b_max = 0.1,
                            minimum_age = 20))
  expect_equal(gomp.a0_output, c(-61.9776483, -2.3090956, 0.0826092))
})
