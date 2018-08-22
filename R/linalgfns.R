#' @title Pivoted Cholesky factorisation
#'
#' @description Pivoted Cholesky (same as pivot=T when carrying out Cholesky with dense matrices)
#'
#' @param A matrix (sparse or dense), the Cholesky factor of which needs to be found
#' @details This function should just be used to verify the pivot=T with dense matrices. Since it is an R implementation it is much slower than that in the base package.
#' @return A list with two elements, R (the Cholesky factor) and piv (the pivoting order)
#' @keywords Cholesky factor
#' @export
#' @examples
#' find_pivot_Bastos(matrix(c(0.1,0.2,0.2,1),2,2))
#' @references Leonardo S. Bastos and A. O'Hagan (2007). Diagnostics for Gaussian Process Emulators. \url{www.tonyohagan.co.uk/academic/pdf/diagtech.pdf}
find_pivot_Bastos <- function(A) {
 n <- nrow(A)
 R <- matrix(0,n,n)
 piv=1:n
 for (k in 1:(n-1)) {
  q <- which.max(diag(A[k:n,k:n])) + k - 1
  A <- swapcol(A,k,q)
  R <- swapcol(R,k,q)
  A <- swaprow(A,k,q)
  piv <- swap(piv,k,q)
  R[k,k] <- sqrt(A[k,k])
  R[k,((k+1):n)] <- (1/R[k,k]) * A[k,((k+1):n)]
  A[((k+1):n),((k+1):n)] <- A[((k+1):n),((k+1):n)] - outer(R[k,((k+1):n)],R[k,((k+1):n)])
 }
 R[n,n] <- sqrt(A[n,n])
 return(list(R=R,piv=piv))
}

#' @title Create an empty matrix
#'
#' @description Creates an empty sparse matrix of size 0 x 0
#' @export
#' @examples
#' require(Matrix)
#' Q <- emptySp()
emptySp <- function() {
  as(matrix(0,0,0),"dgCMatrix")
}

#' @title Create a sparse identity matrix
#'
#' @description Creates a sparse identity matrix of size n x n
#' @param n size of matrix
#' @export
#' @examples
#' require(Matrix)
#' Q <- Imat(4)
Imat <- function(n) {
  sparseMatrix(i=1:n, j=1:n,x=1)
}

#' @title Create an empty sparse matrix
#'
#' @description Creates an empty sparse matrix of size ni x nj
#' @param ni number of rows
#' @param nj number of columns. If NULL a square matrix is produced
#' @export
#' @examples
#' require(Matrix)
#' Q <- Zeromat(2,5)
Zeromat <- function(ni,nj=NULL) {
  if(is.null(nj)) nj <- ni
   return(as(sparseMatrix(i={},j={},dims=c(ni,nj)),"dgCMatrix"))
}


#' @title Find the log determinant
#'
#' @description Find the log determinant of a matrix Q from its Cholesky factor L (which could be permutated or not)
#' @param L the Cholesky factor of Q
#' @examples
#' require(Matrix)
#' Q <- sparseMatrix(i=c(1,1,2,2),j=c(1,2,1,2),x=c(0.1,0.2,0.2,1))
#' logdet(chol(Q))
logdet <- function(L)  {
  ## Find the log-determinant of Q from its Cholesky L
  diagL <- diag(L)
  return(2*sum(log(diagL)))
}


#' @title Create a sparse diagonal matrix
#'
#' @description Creates a sparse diagonal matrix of length(xx) where xx is the vector containing the elements on the diagonal
#' @param xx diagonal vector
#' @export
#' @examples
#' require(Matrix)
#' Q <- sparsediag(c(1,2,3,4))
sparsediag <- function(xx) {
  n <- length(xx)
  return(sparseMatrix(i=1:n,j=1:n,x=xx))
}


# Internal functions
swapcol <- function(A,p,q) {
 temp <- A[,p]
 A[,p] <- A[,q]
 A[,q] <- temp
 return(A)
}

swap <- function(v,p,q) {
  temp <- v[p]
  v[p] <- v[q]
  v[q] <- temp
  return(v)
}

swaprow <- function(A,p,q) {
  temp <- A[p,]
  A[p,] <- A[q,]
  A[q,] <- temp
  return(A)
}

# Find cholPermute using the spam package
.spam_chol <- function(Q,amd=T) {
  Qspam <- as.spam.dgCMatrix(Q)
  if(amd) {
    P <- linalg:::amd_Davis(Q)
    X  <- spam::chol(Qspam,pivot=P)
  } else {
    X  <- spam::chol(Qspam)
  }
  P <- sparseMatrix(i=X@pivot,j=1:nrow(X),x=1)
  Qpermchol <- as(as.dgCMatrix.spam(t(X)),"dtCMatrix")
  return(list(Qpermchol = Qpermchol,P=P))
}


cholMATLAB <- function(Q,matlab_server) {
  Q <- as(Q,"dgTMatrix")
  i = Q@i+1
  j = Q@j+1
  x = Q@x
  setVariable(matlab_server, i=i,j=j,x=x)
  cat("Doing Cholesky in MATLAB",sep="\n")
  evaluate(matlab_server, paste("Q = sparse(i,j,x);",
                                "L = (chol(Q));"))
  L <- getVariable(matlab_server,"L")$L
  L <- as(L,"dtCMatrix")
  cat("Finished Cholesky in MATLAB",sep="\n")
  evaluate(matlab_server, "clearvars Q L i j x")
  return(L)

}

# deprecated Takahashi trials
#----------------------------
Takahashidiag <- function(L) {
  ## Takahashi diag: Find marginal variance from Cholesky factor
  # For now convert to full matrices for indexing. We need to do this intelligently in the future for matrices > 10000x10000 which would fill up memory
  n <- dim(L)[1]
  X <- which(L!=0,arr.ind=T)
  i_ind <- X[,1]
  j_ind <- X[,2]
  Sigma <- Sigma + t(Sigma) - sparseMatrix(i=1:n,j=1:n,x=diag(Sigma))
  Sigma <- as(Sigma,"dgTMatrix") # coerce to i,j format
  Sigma[n,n] <- 1/L[n,n]^2
  numcomp = 0

  # Sigma <- as.matrix(Sigma)
  # L <- as.matrix(L)

  tic()
  for (i in seq(n-1,1,-1)) {
    Lii <- L[i,i]
    nz_indices <- intersect(i_ind[j_ind == i],(i+1):n)
    Lip1n_i<-L[nz_indices,i]
    for (j in intersect(seq(n,i,-1),which(abs(L[,i])>0))) {
      Sigma[i,j] <-  Sigma[j,i] <- (i==j)/(Lii^2) - 1/Lii*sum(Lip1n_i*Sigma[j,nz_indices ])
      numcomp = numcomp + 1
    }
  }
  toc()



  return(diag(Sigma))
}
Takahashidiag_Cseke <- function(L) {

  n <- dim(L)[1]
  invdiagL2 <- 1/diag(L)^2
  S <- L + t(L)
  S[S >0] = 1
  S[n,n] =invdiagL2[n]

  for (i in seq(n-1,1,-1)) {
    I   = i+which(abs(L[seq(i+1,n,1),i]) > 0)
    S[I,i] = -(S[I,I]%*%L[I,i])/L[i,i]
    S[i,I] = t(S[I,i])
    S[i,i] = invdiagL2[i] - (S[i,I]%*%L[I,i])/L[i,i];

  }
  return(diag(S))
}

Takahashi <- function(L, diag = TRUE, method = c("R", "C")) {
  #### R and C implementation of Takahashi diagonal equations by Jonty

  method <- match.arg(method)
  stopifnot(inherits(L, "CsparseMatrix")) # has @i and @p


  n <- nrow(L)
  ii <- L@i + 1 # in {1,...,n}
  dp <- diff(L@p)
  jj <- rep(seq_along(dp), dp) # in {1,...,n}, non-decreasing
  N = length(ii)

  stopifnot(ii >= jj,              # lower triangular
            1:n %in% ii[ii == jj]) # full diagonal


  if (method == "C") {
    tic();
    dyn.load("Test.so")
    X <- .C("TakahashiC",as.integer(n),as.integer(N),as.integer(ii),as.integer(jj),as.double(L@x),results = double(N))
    toc();
    return(X$results[ii==jj])
  } else if (method == "R") {

    if (diag) { # this to speed up the calculation
      tic()
      S <- L; S@x[] <- -999

      S@x[ii == n & jj == n] <- 1 / L[n, n]^2

      if (n > 1)
        for (i in (n-1):1) {

          k <- ii[ii > i & jj == i]      # Find row numbers with non zero indices at this column
          if (length(k) == 0) {
            S@x[ii == i & jj == i] <- 1 / L[i, i]^2
          } else {
            Lii <- L[i, i]
            Lki <- L[k, i]

            js <- rev(jj[ii %in% k & jj >= i]) # going backwards
            #for (j in js) {
            #  skj <- S@x[ii == pmax(k, j) & jj == pmin(k, j)] # select from lower triangle
            #  S@x[ii == j & jj == i] <- ((i==j) / Lii - sum(Lki * skj)) / Lii
            #}
            js = unique(js)
            js <- c(k,i)
            for(j in js) {
              skj <- apply(matrix(k),1,function(ind){ S@x[ii == max(ind,j) & jj == min(j,ind)] } )
              S@x[ii == j & jj == i] <- ((i==j) / Lii - sum(Lki * skj)) / Lii
            }

          }
        }
      toc()
      return(diag(S))

    } else { # diag = FALSE : use full-size S but only fill lower triangle

      S <- matrix(NA, n, n)

      S[n, n] <- 1 / L[n, n]^2

      if (n > 1)
        for (i in (n-1):1) {

          k <- ii[ii > i & jj == i]
          if (length(k) == 0) {
            S@x[ii == i & jj == i] <- 1 / L[i, i]^2
          } else {
            Lii <- L[i, i]
            Lki <- L[k, i]

            js <- n:i # going backwards

            for (j in js) {
              skj <- S[pmax(k, j), pmin(k, j)] # select from lower triangle
              S[j, i] <- ((i==j) / Lii - sum(Lki * skj)) / Lii
            }
          }
        }

      return(ifelse(is.na(S), t(S), S))
    }

  } else stop("Never get here!")
}
