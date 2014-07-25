\docType{methods}
\name{plot_interp}
\alias{plot,FEBasis,ANY-method}
\alias{plot_interp}
\alias{plot_interp,FEBasis,character-method}
\alias{plot_interp,FEBasis-method}
\title{Plots an interpolated field from a mesh}
\usage{
plot_interp(x, y, ds, ...)

\S4method{plot_interp}{FEBasis,character}(x, y, ds = 40, ...)

\S4method{plot}{FEBasis,ANY}(x, y, ...)
}
\arguments{
  \item{x}{a Mesh object}

  \item{y}{a character indicated which column to plot}

  \item{ds}{the resolution of the plot (defaulted to 40 x
  40)}

  \item{max}{upperbound on colour scale}

  \item{min}{lowerbound on colour scale}

  \item{...}{Arguments to be passed to methods, such as
  \link{graphical parameters} (see \code{\link{par}}).
  Many methods will accept the following arguments:
  \describe{ \item{\code{type}}{what type of plot should be
  drawn.  Possible types are \itemize{ \item \code{"p"} for
  \bold{p}oints, \item \code{"l"} for \bold{l}ines, \item
  \code{"b"} for \bold{b}oth, \item \code{"c"} for the
  lines part alone of \code{"b"}, \item \code{"o"} for both
  \sQuote{\bold{o}verplotted}, \item \code{"h"} for
  \sQuote{\bold{h}istogram} like (or \sQuote{high-density})
  vertical lines, \item \code{"s"} for stair \bold{s}teps,
  \item \code{"S"} for other \bold{s}teps, see
  \sQuote{Details} below, \item \code{"n"} for no plotting.
  } All other \code{type}s give a warning or an error;
  using, e.g., \code{type = "punkte"} being equivalent to
  \code{type = "p"} for S compatibility.  Note that some
  methods, e.g. \code{\link{plot.factor}}, do not accept
  this.  } \item{\code{main}}{an overall title for the
  plot: see \code{\link{title}}.} \item{\code{sub}}{a sub
  title for the plot: see \code{\link{title}}.}
  \item{\code{xlab}}{a title for the x axis: see
  \code{\link{title}}.} \item{\code{ylab}}{a title for the
  y axis: see \code{\link{title}}.} \item{\code{asp}}{the
  \eqn{y/x} aspect ratio, see \code{\link{plot.window}}.} }
  }
}
\value{
a ggplot2 object
}
\description{
This function takes a Mesh and a column name in the mesh,
to generate an interpolated field of the indicated column
}
\examples{
\dontrun{
data(surf_fe)
Mesh <- initFEbasis(p=surf_fe$p, t = surf_fe$t, M = surf_fe$M, K = surf_fe$K)
Mesh["z"] <- sin(Mesh["x"]/1000)*cos(Mesh["y"]/1000)
g <- plot_interp(Mesh,"z",ds=200)}
}

