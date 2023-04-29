\name{array-stuff}
\alias{ul}
\alias{mkAL}
\alias{saveSim}
\alias{maybeRead}
\alias{getArray}
\alias{array2df}
\title{Tools For Converting To and From Arrays, Lists, and Array of Lists}
\description{
  \describe{
    \item{ul()}{is a simple wrapper for \code{\link{unlist}()} with
      \code{recursive=FALSE}.}
    \item{mkAL()}{gets a list \code{x} with elements that are named
      lists of length five, see \code{x} below, and converts it to an
      array of \code{\link{list}}s.}
    \item{saveSim()}{(optionally) converts a result list to an array using
      \code{mkAL()} and (optionally) saves it to a file via \code{\link{saveRDS}()}.}
    \item{maybeRead()}{if the provided \file{.rds} file exists, this function
      reads it via \code{\link{readRDS}()}; otherwise, nothing is done.}
    \item{getArray()}{gets an array of 4-\code{\link{list}}s as returned
      by \code{mkAL()}, picks out the specified component \code{comp},
      applies the specified function \code{FUN} (with useful defaults),
      and builds an \code{\link{array}}.}
    \item{array2df()}{auxiliary function to convert an array to a
      \code{\link{data.frame}} (correctly dealing with \code{n.sim}).}
  }
}
\usage{
ul(x)

mkAL   (x, vList, repFirst, check=TRUE)
saveSim(x, vList, repFirst, sfile, check=TRUE, doAL=TRUE)

maybeRead(sfile, msg=TRUE)

getArray(x, comp=c("value", "error", "warning", "time"),
         FUN=NULL, err.value=NA)

array2df(x, responseName = "value")
}
\arguments{
  \item{x}{for
    \describe{
      \item{ul()}{a list.}
      \item{mkAL(), saveSim()}{a list (of length
        \code{n.sim * nrow(pGrid)})
	where each element is a list of length
	five, containing the named elements \code{"value"},
         \code{"error"}, \code{"warning"}, \code{"time"}, and
         \code{".Random.seed"}, the first four as returned by
         \code{\link{doCallWE}()}.}
      \item{getArray()}{an array of \code{\link{list}}s as
	returned by \code{mkAL()}.}
      \item{array2df()}{a numeric array as returned by
	\code{getArray(*, "value")}.}
    }
  }
  \item{vList}{a \code{\link{list}} of variable specifications.  Each
    variable specification is itself a named list which must contain a
    \code{"value"} component.}
  \item{repFirst}{logical; \bold{must} match the value of \code{repFirst}
    in the \code{x <- do*Apply()} call where \code{x} has been created.}
  \item{check}{\code{\link{logical}} activating consistency checks for
    \code{x}.}
  \item{sfile}{a file name, typically with extension \file{.rds} or \code{\link{NULL}}.}
  \item{doAL}{\code{\link{logical}} indicating if \code{mkAL()} should
    be called, or rather just \code{x} be saved.}
  \item{msg}{\code{\link{logical}} indicating whether a message is
    printed when an object is read from \code{sfile}.}
  \item{err.value}{\code{\link{numeric}} which is used to replace the
    value of the array entry in case of an error.}
  \item{comp}{character string denoting the component.}
  \item{FUN}{function to be applied right before the resulting array
    \code{\link{array}} is constructed.}
  \item{responseName}{(for \code{arrady2df()}:) a string specifying the
    name of the \dQuote{value} column of the resulting data frame.}
}
\value{
  For \describe{
    \item{ul()}{the unlisted list; see \code{\link{unlist}()}.}
    \item{mkAL()}{an \code{\link{array}} of \code{\link{list}}s.}
    \item{saveSim()}{the \code{\link{array}} returned by
      \code{mkAL()}.}
    \item{maybeRead()}{the object read by \code{\link{readRDS}()} from
      \code{sfile} or nothing (if \code{sfile} does not exist).}
    \item{getArray()}{an \code{\link{array}} containing the values
      of the specified component \code{comp} after applying \code{FUN}
      to them.  The default \code{FUN} produces an \code{\link{array}},
      depending on \code{comp}, of
      \describe{
	\item{"value":}{values or \code{err.value} (in case of an error)}
	\item{"error":}{logicals indicating whether there was an error}
	\item{"warning":}{logicals indicating whether there was a warning}
	\item{"time":}{timings as returned by \code{\link{doCallWE}()},
	  i.e., typically (from \code{\link{mkTimer}}'s
	  \code{\link{proc.time}()[1]}) the number of milliseconds
	  of \dQuote{"CPU user time"}.}
      }
    }
    \item{array2df(x)}{a \code{\link{data.frame}} with several columns
      built from the \code{\link{dimnames}(x)} and a column named
      \code{responseName} with the values of \code{x}.}
  }
}
\details{
  \code{mkAL()} is useful when creating arrays from result
  lists returned from large(r) simulation studies which use
  \code{\link{doCallWE}()}. To create a proper argument \code{x} for
  \code{mkAL()}, the function \code{ul()} turns out to
  be useful to (stepwise) unlist nested lists.

  \code{getArray()} converts arrays of lists as returned by
  \code{mkAL()} to an \code{\link{array}} of \code{\link{numeric}} (or
  \code{\link{logical}}, see below) after applying the specified \code{FUN}.
  \cr%%--------------------
  In case of an error, the corresponding entry in the resulting array is
  replaced by \code{err.value}.

  The default \code{FUN} converts possible errors and warnings to
  \code{\link{logical}} (indicating whether there was a error or
  warning, respectively) and run times to \code{\link{numeric}}.  For
  \code{comp="value"}, the situation is trickier.  First of all, the
  resulting array contains dimensions for variables of type
  \dQuote{inner} and, if greater than 1, for the variable of type
  \dQuote{N} (typically called \code{"n.sim"}); see the vignette for
  details.
  \cr%----------------
  Use \code{FUN = identity} to get at the full error or warning objects,
  for \code{comp = "error"} or for \code{comp = "warning"}, respectively.

  \code{saveSim()} and \code{maybeRead()} are
  useful for creating and (re)storing arrays from large(r) simulation
  studies (to avoid recomputation, to ease the data analysis etc.).
  \code{saveSim()} calls \code{mkAL()}, nowadays wrapped in
  \code{\link{tryCatch}(.)}, such that the simulation is not lost, even
  when the resulting format cannot correctly be treated by
  \code{mkAL()}.  Consequently, \code{doAL} is not much needed anymore.
  Note that both \code{saveSim()} and \code{maybeRead()} accept
  \code{sfile=\link{NULL}} in which case nothing is saved or read.
}
\author{Marius Hofert and Martin Maechler.}
\references{see \code{\link{simsalapar-package}}.}
\seealso{
  \code{\link{getEl}()} and \code{\link{mkNms}()} used by
  \code{mkAL()}. \code{\link{saveRDS}()} and
  \code{\link{readRDS}()}, the \dQuote{workhorses} of
  \code{saveSim()} and \code{maybeRead()}, respectively.
}
\examples{
\dontrun{
## Get at the full error objects, notably (message, call):
errObjs <- getArray(res, "error", FUN=identity)
}

if(FALSE) ## A longer, "interesting" example is in
 demo(robust.mean)
}
\keyword{utilities}
