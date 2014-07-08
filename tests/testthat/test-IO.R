context("input-output operations")
test_that("output is correct for a given input", {
  Q = sparseMatrix(i=c(1,1,2,2),j=c(1,2,1,2),x=c(0.1,0.2,0.2,1))
  X <- cholPermute(Q)
  S_partial = Takahashi_Davis(Q,cholQp = X$Qpermchol,P=X$P)
  expect_that(S_partial,is_a("dgCMatrix"))
  expect_that(diag(S_partial),equals(diag(solve(Q))))
  expect_that(dim(S_partial),equals(c(2,2)))
  expect_that(logdet(chol(Q)),equals(log(det(Q))))
  
  Q2 <- as(Q,"matrix") # now make it dense
  pivot_chol1 <- chol(Q2,pivot=1)
  pivot_chol2 <- find_pivot_Bastos(Q2)
  expect_that(attr(pivot_chol1,"pivot"),equals(pivot_chol2$piv))
  expect_that(matrix(pivot_chol1),equals(matrix(pivot_chol2$R)))
})