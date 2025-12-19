data_list <- list(
  lapply(miheart, \(x) x[, 4, 1]),
  lapply(miheart, \(x) x[, 4:6, 1]),
  lapply(miheart, \(x) x[, 4, 1:3]),
  lapply(miheart, \(x) x[, 4:6, 1:3])
)

test_that("all CAR models work", {
  for (x in data_list) {
    expect_no_error(car("test", x, miadj, show_plots = FALSE, seed = 1234))
  }
})

test_that("all RCAR models work", {
  for (x in data_list) {
    expect_no_error(rcar("test", x, miadj, show_plots = FALSE, seed = 1234))
  }
})

test_that("all MCAR models work", {
  for (x in data_list) {
    expect_no_error(mcar("test", x, miadj, show_plots = FALSE, seed = 1234))
  }
})

test_that("all MSTCAR models work", {
  for (x in data_list) {
    expect_no_error(mstcar("test", x, miadj, show_plots = FALSE, seed = 1234))
  }
})

test_that("all MSTCAR models work with update_rho = TRUE", {
  for (x in data_list) {
    expect_no_error(
      mstcar(
        "test",
        x,
        miadj,
        update_rho = TRUE,
        show_plots = FALSE,
        seed = 1234
      )
    )
  }
})
