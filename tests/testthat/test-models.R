data_list <- list(
  lapply(miheart, \(x) x[, 4, 1]),
  lapply(miheart, \(x) x[, 4:6, 1]),
  lapply(miheart, \(x) x[, 4, 1:3]),
  lapply(miheart, \(x) x[, 4:6, 1:3])
)

test_that("all UCAR models work", {
  lapply(data_list, \(x) {
    expect_no_error(
      ucar("test", x, miadj, show_plots = FALSE)
    )
  })
})

test_that("all EUCAR models work", {
  lapply(data_list, \(x) {
    expect_no_error(
      eucar("test", x, miadj, show_plots = FALSE)
    )
  })
})

test_that("all MCAR models work", {
  lapply(data_list, \(x) {
    expect_no_error(
      mcar("test", x, miadj, show_plots = FALSE)
    )
  })
})

test_that("all MSTCAR models work", {
  lapply(data_list, \(x) {
    expect_no_error(
      mstcar("test", x, miadj, show_plots = FALSE)
    )
  })
})

test_that("all MSTCAR models work with rho updates", {
  lapply(data_list, \(x) {
    expect_no_error(
      mstcar("test", x, miadj, update_rho = TRUE, show_plots = FALSE)
    )
  })
})
