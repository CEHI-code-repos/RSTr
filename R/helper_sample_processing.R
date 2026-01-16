create_array_new <- function(sample_new, sample, margin, bind_new, new_name) {
  if (bind_new) {
    new_dim <- dim(sample)
    new_dim[margin] <- 1
    new_dimnames <- dimnames(sample)
    new_dimnames[[margin]] <- c(new_dimnames[[margin]], new_name)
    agg_sample <- array(sample_new, dim = new_dim)
    array_new <- abind::abind(sample, agg_sample, along = margin)
    dimnames(array_new) <- new_dimnames
  } else {
    new_dim <- dim(sample)[-margin]
    new_dimnames <- dimnames(sample)[-margin]
    array_new <- array(sample_new, dim = new_dim, dimnames = new_dimnames)
  }
  array_new
}

arr_to_matrix <- function(array, perm, nrow, ncol) {
  array |> aperm(perm) |> matrix(nrow = nrow, ncol = ncol)
}
