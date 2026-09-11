test_that("current exact and sparse scientific baselines remain unchanged", {
  # Captured from authoritative 156c60d, before correctness edits.
  # These intentionally differ: legacy L1/Ward.D2 vs Euclidean sparse Ward.
  set.seed(17)
  x <- array(10+rnorm(4*5*9),c(4,5,9))
  exact <- c(1L,1L,1L,1L,2L,1L,1L,2L,1L,1L,1L,2L,3L,3L,2L,1L,2L,1L,1L,2L)
  sparse <- c(1L,1L,1L,1L,1L,2L,2L,1L,2L,2L,2L,1L,3L,3L,1L,1L,2L,1L,2L,1L)
  partition <- function(x)outer(as.vector(x),as.vector(x),"==")
  expect_equal(partition(segment(x,Ncomp=3)$cluster_map),partition(exact))
  z <- segment_large(x,Ncomp=3,knn_k=8,verbose=FALSE)
  expect_equal(partition(z$cluster_map),partition(sparse))
  expect_equal(z$original_cube$imDat,x)
  reconstructed <- reconstruct_flux_preserving_cube(z)
  expect_lt(reconstructed$flux_check$max_abs_difference,1e-10)
})
