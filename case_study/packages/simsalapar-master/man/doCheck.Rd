\name{doCheck}
\title{Checking a User's doOne}
\alias{doCheck}
\description{
  \code{doCheck()} checks, if possible, a user's doOne() function for
  return objects of correct sizes.
}
\usage{
doCheck(doOne, vList, nChks = ng, verbose = TRUE)
}
\arguments{
  \item{doOne}{a user-supplied \code{\link{function}} for computing one
    row of the (physical) grid.}
  \item{vList}{a \code{\link{list}} of variable specifications. Each
    variable spec is itself a named list which must contain a
    \code{"value"} component.}
  \item{nChks}{number of rows randomly picked from the (physical) grid
    which are used for a basic test of the evaluation and return value
    of \code{doOne()}.}
  \item{verbose}{\code{\link{logical}} indicating whether check output
    is displayed.}
}
\value{
  None.
}
\author{Marius Hofert and Martin Maechler.}
\examples{
## definition
doCheck
}
\keyword{utilities}
