\name{doCallWE}
\title{Innermost Computation: Error Catching Version of do.call()}
\alias{doCallWE}
\alias{mkTimer}
\description{
  \code{doCallWE()} performs the innermost computation of the simulation
  study at hand.  It is a version of \code{\link{do.call}(f, argl, *)},
  with care of catching and storing \emph{both} error and warnings (via
  \code{\link{tryCatch.W.E}()}) and measures user time.  This is useful
  in large(r) simulation studies.

  \code{mkTimer()} returns a \emph{function} to be passed as
  \code{timer} to \code{doCallWE()}.
}
\usage{
doCallWE(f, argl, % f.value = NA_real_; introduce, if required
         timer = mkTimer(gcFirst=FALSE))

mkTimer(gcFirst)
}
\arguments{
  \item{f}{a \code{\link{function}} which given data and parameters,
    computes the statistic we are simulating.}
  \item{argl}{list of arguments for \code{f()}.}
  % \item{f.value}{prototype return value of \code{f()}.}
  \item{timer}{a \code{\link{function}} similar to
    \code{\link{system.time}()}; by default, measure user time in
    milliseconds.}

  \item{gcFirst}{logical, passed to \code{\link{system.time}()}, as it
    is called from the resulting function \code{mkTimer()}.}
}
\value{
  \code{doCallWE()} returns a \code{\link{list}} with components
  \item{value}{\code{f}\eqn{(\langle \mathtt{argl} \rangle)}{(<argl>)},
    if there was no error, \code{\link{NULL}} otherwise.}
  \item{error}{error message (see \code{\link{simpleError}} or
    \code{\link{stop}()}), \code{\link{NULL}} otherwise.}
  \item{warning}{warning message (see \code{\link{simpleWarning}} or
    \code{\link{warning}()}), \code{\link{NULL}} otherwise.}
  \item{time}{time, as measured by \code{timer()}; defaults to
    milliseconds without garbage collection.}
}
\details{
  Note that \code{gcFirst=FALSE} is default for a good reason: if a call
  to \code{doOne()} is relatively fast, calling \code{\link{gc}()} every
  time is unnecessarily expensive and may completely dominate the
  overall simulation run time. For serious run time measurement,
  \code{gcFirst=TRUE} is preferable, as it ensures less
  variable timings, see \code{\link{system.time}}.
}
\author{Marius Hofert and Martin Maechler.}
\seealso{
  \code{\link{do.call}}, \code{\link{tryCatch.W.E}}.
}
\examples{
set.seed(61)
L <- log(abs(rt(n=100, df = 1.5)))
r <- doCallWE(quantile, list(L, probs= 0.95))
## set timer for  "no timing" :
u <- doCallWE(quantile, list(L, probs= 0.95), timer = function(E) { E; NULL })
stopifnot(is.null(r$error),
	  all.equal(r$value, quantile(L, 0.95)),
	  identical(r[1:3], u[1:3]), is.null(u[["time"]]))
}
\keyword{programming}
