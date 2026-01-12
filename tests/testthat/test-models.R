data_list <- list(
  lapply(miheart, \(x) x[, 4, 1]),
  lapply(miheart, \(x) x[, 4:6, 1]),
  lapply(miheart, \(x) x[, 4:6, 1:3])
)

test_that("all CAR models work", {
  expect_warning(
    car(
      "test",
      data_list[[1]],
      miadj,
      show_plots = FALSE,
      verbose = FALSE,
      seed = 1234
    )
  )
})

test_that("all RCAR models work", {
  expect_no_error(
    rcar(
      "test",
      data_list[[1]],
      miadj,
      show_plots = FALSE,
      verbose = FALSE,
      seed = 1234
    )
  )
})

test_that("all MCAR models work", {
  expect_no_error(
    mcar(
      "test",
      data_list[[2]],
      miadj,
      show_plots = FALSE,
      verbose = FALSE,
      seed = 1234
    )
  )
})

test_that("all MSTCAR models work", {
  expect_no_error(
    mstcar(
      "test",
      data_list[[3]],
      miadj,
      show_plots = FALSE,
      verbose = FALSE,
      seed = 1234
    )
  )
})

test_that("all MSTCAR models work with update_rho = TRUE", {
  expect_no_error(
    mstcar(
      "test",
      data_list[[3]],
      miadj,
      update_rho = TRUE,
      show_plots = FALSE,
      verbose = FALSE,
      seed = 1234
    )
  )
})
