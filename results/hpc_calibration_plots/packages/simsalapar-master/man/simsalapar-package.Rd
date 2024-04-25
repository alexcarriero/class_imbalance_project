\name{simsalapar-package}
\alias{simsalapar-package}
\alias{simsalapar}
\docType{package}
\title{Tools for Simulation Studies in Parallel with \R}
\description{
  Tools for setting up, conducting, and evaluating larger-scale
  simulation studies, including parallel computations, in \R.
}
\details{

  The DESCRIPTION file:
  \packageDESCRIPTION{simsalapar}
  \packageIndices{simsalapar}

  \subsection{Setting up a simulation}{
    \describe{
      \item{\code{\link{varlist}}()}{creates a variable specification list.}
      \item{\link{dimnames2varlist}()}{creates a variable specification list from
	given dimension names.}
      \item{\link{getEl}()}{extracts elements from a variable list.}
      \item{mkGrid()}{function for creating a grid of all variables of type
	\dQuote{grid}; see \code{\link{mkGrid}()}.}
      \item{mkNms()}{builds a list of names from a variable list; see
	\code{\link{mkNms}()}.}
      \item{get.n.sim()}{extracts \dQuote{n.sim}; see \code{\link{get.n.sim}()}.}
      \item{get.nonGrids()}{extracts all variables not of type
	\dQuote{grid}; see \code{\link{get.nonGrids}()}.}
    }
  }
  \subsection{Conducting a simulation}{
    \describe{
      \item{tryCatch.W.E()}{catching and storing warnings and errors
	simultaneously; see \code{\link{tryCatch.W.E}()}.}
      \item{doCallWE()}{innermost computation (return value of doOne()):
	returns value, error, warning, and run time; see \code{\link{doCallWE}()}.}
      \item{LEseeds()}{create a list of advanced \code{.Random.seed}'s for
	\dQuote{L'Ecuyer-CMRG}; see \code{\link{LEseeds}()}.}
      \item{printInfo()}{displays information about the sub-job just
	finished; see \code{\link{printInfo}()}.}
      \item{subjob()}{computes one row of the virtual grid in a
	simulation; see \code{\link{subjob}()}.}
      \item{mkTimer()}{creates a function to be passed to doCallWE() as
	timer; see \code{\link{mkTimer}()}.}
      \item{\code{\link{doLapply}()}}{ sequentially iterates over all
	subjobs via standard \code{\link{lapply}()}.}
      \item{\code{\link{doForeach}()}}{ iterates over all subjobs in parallel (via
	\code{foreach()}, package \CRANpkg{foreach}).}
      \item{\code{\link{doRmpi}()} }{iterates over all subjobs in parallel
	(via \CRANpkg{Rmpi}'s \code{mpi.apply()}).}
      \item{\code{\link{doMclapply}()} }{iterates over all subjobs in parallel (via
	\code{\link{mclapply}()}).}
      \item{\code{\link{doClusterApply}()} }{iterates over all subjobs
	in parallel (via \code{\link{clusterApply}()}).}
    }
  }
  \subsection{Analysis}{
    \describe{
      \item{doRes.equal()}{convenience wrapper for comparing two results
	of the \code{do*} lapply-like functions; see \code{\link{doRes.equal}()}.}
      \item{mkAL()}{converts a list of named 5-lists to an array of lists;
	see \code{\link{mkAL}()}.}
      \item{saveSim()}{(optionally) converts a result list to an
        array of lists using mkAL(); see \code{\link{saveSim}()}.}
      \item{maybeRead()}{(optionally) reads the provided \code{.rds}; see
	\code{\link{maybeRead}()}.}
      \item{getArray()}{gets an array of 4-lists and computes an
	array of values, errors, warnings, or run times; see
	\code{\link{getArray}()}.}
      \item{\code{\link{array2df}()} }{conveniently converts an array to
	a \code{\link{data.frame}}.}
      \item{\code{\link{toLatex}()}:}{an S3 method for varlist and ftable.}
      \item{\code{\link{fftable}()} }{essentially calls
	\code{format.ftable()} and adds attributes \code{ncv} and
	\code{nrv} to the return object.}
      \item{\code{\link{tablines}()} }{computes ingredients for converting a
	character matrix with attributes to a LaTeX table.}
      \item{\code{\link{wrapLaTable}()}}{wraps a table and tabular
	environment around the lines of the body of a LaTeX table.}
      \item{\code{\link{mayplot}()}:}{a matrix-like plot for arrays up
	to rank 5, with grid and gridBase.}
    }
  }
}
\author{
Marius Hofert and Martin Maechler <maechler@stat.math.ethz.ch>

Maintainer: Marius Hofert <marius.hofert@math.ethz.ch>
}
\references{
  \describe{
    \item{Publication}{
      Marius Hofert, Martin Maechler (2016).
      Parallel and Other Simulations in \R Made Easy: An End-to-End Study.
      \emph{Journal of Statistical Software}, \bold{69}(4), 1--44.
      \doi{10.18637/jss.v069.i04}}

    \item{Preprint (for simsalapar 1.0-0; including timing info):}{
      Hofert, M. and \enc{Mächler}{Maechler}, M. (2013).
      Parallel and other simulations in \R made easy: An end-to-end study.
      \url{https://arxiv.org/abs/1309.4402}}
  }
}
\examples{
\dontrun{ % we have a faster version as test
 demo(TGforecasts)
}
}
\keyword{package}
