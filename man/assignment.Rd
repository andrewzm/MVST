\docType{methods}
\name{[<--methods}
\alias{[<-,Basis,ANY,ANY,ANY-method}
\alias{[<-,Basis,ANY,ANY-method}
\alias{[<-,GMRF,ANY,ANY,ANY-method}
\alias{[<-,GMRF,ANY,ANY-method}
\alias{[<-,GMRF_basis,ANY,ANY,ANY-method}
\alias{[<-,GMRF_basis,ANY,ANY-method}
\alias{[<-,Obs,ANY,ANY,ANY-method}
\alias{[<-,Obs,ANY,ANY-method}
\alias{[<--methods}
\title{Assignment methods}
\usage{
\S4method{[}{GMRF_basis,ANY,ANY,ANY}(x, i, j) <- value

\S4method{[}{GMRF,ANY,ANY,ANY}(x, i, j) <- value

\S4method{[}{Obs,ANY,ANY,ANY}(x, i, j) <- value

\S4method{[}{Basis,ANY,ANY,ANY}(x, i, j) <- value
}
\arguments{
  \item{x}{ object from which to extract element(s) or in
  which to replace element(s).  }

  \item{i}{ indices specifying elements to extract or
  replace.  Indices are \code{numeric} or \code{character}
  vectors or empty (missing) or \code{NULL}.  Numeric
  values are coerced to integer as by
  \code{\link{as.integer}} (and hence truncated towards
  zero).  Character vectors will be matched to the
  \code{\link{names}} of the object (or for
  matrices/arrays, the \code{\link{dimnames}}): see
  \sQuote{Character indices} below for further details.

  For \code{[}-indexing only: \code{i}, \code{j},
  \code{\dots} can be logical vectors, indicating
  elements/slices to select.  Such vectors are recycled if
  necessary to match the corresponding extent.  \code{i},
  \code{j}, \code{\dots} can also be negative integers,
  indicating elements/slices to leave out of the selection.

  When indexing arrays by \code{[} a single argument
  \code{i} can be a matrix with as many columns as there
  are dimensions of \code{x}; the result is then a vector
  with elements corresponding to the sets of indices in
  each row of \code{i}.

  An index value of \code{NULL} is treated as if it were
  \code{integer(0)}.  }

  \item{value}{typically an array-like \R object of a
  similar class as \code{x}.}
}
\description{
Methods for \code{"[<-"}, i.e., extraction or subsetting of
elements in the data frame of the MVST object
}
\examples{
data(surf_fe)
Mesh <- initFEbasis(p=surf_fe$p, t = surf_fe$t, M = surf_fe$M, K = surf_fe$K)
Mesh["z"] <- sin(Mesh["x"]/1000)*cos(Mesh["y"]/1000)
}

