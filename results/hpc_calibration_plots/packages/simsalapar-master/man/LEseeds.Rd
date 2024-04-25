\name{LEseeds}
\alias{LEseeds}
\title{Advancing .Random.seed for "L'Ecuyer-CMRG"}
\description{
  \code{LEseeds()} creates a list of advanced .Random.seed's for "L'Ecuyer-CMRG".
}
\usage{
LEseeds(n)
}
\arguments{
  \item{n}{number of steps to advance \code{\link{.Random.seed}}.}
}
\value{
  A list of length \code{n} containing the advanced \code{\link{.Random.seed}}'s.
}
\details{
  See, for example, Hofert and \enc{Mächler}{Maechler} (2014) for how
  to use these functions.
}
\author{Marius Hofert and Martin Maechler.}
\seealso{
  \code{\link{.Random.seed}} for information about random number
  generators and seeds.
}
% \examples{
% }
\keyword{utilities}
