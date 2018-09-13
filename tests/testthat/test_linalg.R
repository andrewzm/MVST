context("linalg")
set.seed(1)
a <- matrix(rnorm(10))
A <- matrix(rnorm(25), 5, 5)
B <- tcrossprod(A, A)
C <- matrix(rnorm(100), 10, 10)
D <- t(A)
e <- matrix(rnorm(5))

test_that("logdet works as expected", {

  logdetB_raw <- as.numeric(determinant(B)$modulus)
  logdetB <- logdet(chol(B))
  expect_equal(logdetB_raw, logdetB)

})

test_that("chol works", {

  Q <- GMRF_RW(n = 6)@Q + Diagonal(6)
  L1 <- chol(Q)
  X1 <- cholPermute(Q)
  X2 <- cholPermute_old(Q)

  L1 <- X1$Qpermchol
  P1 <- X1$P

  L2 <- X2$Qpermchol
  P2 <- X2$P

  expect_equal(as.matrix(Q), as.matrix(P1 %*% L1 %*% t(L1) %*% t(P1)))
  expect_equal(as.matrix(Q), as.matrix(P2 %*% L2 %*% t(L2) %*% t(P2)))
      
})
