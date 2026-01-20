data_min <- lapply(miheart, \(x) x[, 4:6, 1:3])

test_that("the MSTCAR model works", {
  expect_no_error(
    mstcar(
      "test",
      data_min,
      miadj,
      show_plots = FALSE,
      verbose = FALSE,
      seed = 1234
    )
  )
})
