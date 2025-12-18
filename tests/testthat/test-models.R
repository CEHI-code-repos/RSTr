run_all_models <- function(data) {
  ucar("test", data, miadj, show_plots = FALSE)
  eucar("test", data, miadj, show_plots = FALSE)
  mcar("test", data, miadj, show_plots = FALSE)
  mstcar("test", data, miadj, show_plots = FALSE)
}

test_that("models run with regions only", {
  data_min <- lapply(miheart, \(x) x[, 4, 1])
  expect_no_error(run_all_models(data_min))
})

test_that("models run with regions + groups", {
  data_min <- lapply(miheart, \(x) x[, 4:6, 1])
  expect_no_error(run_all_models(data_min))
})

test_that("models run with regions + time", {
  data_min <- lapply(miheart, \(x) x[, 4, 1:3])
  expect_no_error(run_all_models(data_min))
})

test_that("models run with regions + groups + time", {
  data_min <- lapply(miheart, \(x) x[, 4:6, 1:3])
  expect_no_error(run_all_models(data_min))
})
