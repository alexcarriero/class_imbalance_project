\name{tryCatch.W.E}
\alias{tryCatch.W.E}
\title{Catching and Storing Warnings and Errors Simultaneously}
\description{
   Catches and saves both warnings (\code{\link{warning}}) and errors
   (\code{\link{stop}}) and in the case of a warning, also the
   computed result.
}
\usage{
tryCatch.W.E(expr)
}
\arguments{
  \item{expr}{expression to be evaluated, typically a function call.}
}
\value{
  \code{\link{list}} with components
  \item{\code{value}}{value of \code{expr} \emph{or} error message
    (see \code{\link{simpleError}} or \code{\link{stop}()}).}
  \item{\code{warning}}{warning  message (see \code{\link{simpleWarning}} or
    \code{\link{warning}()}) or \code{\link{NULL}}.}
}
\details{
  This function is particularly useful in large(r) simulation studies to
  check all computations and guarantee their correctness.
}
\author{
  Marius Hofert and Martin Maechler, based on hints from
  Luke Tierney and Bill Dunlap, see
  \url{https://stat.ethz.ch/pipermail/r-help/2010-December/262626.html}.
}
\references{see \code{\link{simsalapar-package}}.}
\seealso{
  the base function \code{\link{tryCatch}()} and
  \code{demo(error.catching)}. Also, \code{\link{doCallWE}()}, of which
  \code{\link{tryCatch.W.E}()} is the \dQuote{workhorse}.
}
\examples{
## Adapted from demo(error.catching) :
str(r1 <- tryCatch.W.E( log( 2 ) ) )
str(r2 <- tryCatch.W.E( log(-1 ) ) )
str(r3 <- tryCatch.W.E( log("a") ) )
stopifnot(is.null(r1$warning),
          is.na  (r2$value),   inherits(r2$warning, "warning"),
          is.null(r3$warning), inherits(r3$value, "error"))
}
\keyword{programming}
