\name{varlist}
\title{Variable Specification List - Generation and Class}
\alias{varlist}
\alias{dimnames2varlist}
\alias{varlist-class}
\alias{coerce,varlist,list-method}
\alias{show,varlist-method}
\description{
  Generate variable specification lists.  These are objects of the
  formal (aka \dQuote{S4}) class \code{"varlist"}.  This class
  simply extends \code{"\linkS4class{namedList}"} and has a validity
  method (see \code{\link{validObject}}).
}
\usage{
varlist(...)
dimnames2varlist(dmn)
%% \S4method{coerce}{varlist,list}(from, to)
%% as(object, "list")% <- needs \alias{as} and documentation of "list"
\S4method{show}{varlist}(object)
}
\arguments{
  \item{\dots}{of the form
    \describe{
      \item{\code{nam1} = list(....),}{}
      \item{\code{nam2} = list(....),}{}
      \item{........}{}
      \item{\code{namk} = list(....)}{}
    }
    i.e, a \dQuote{list} of variable specifications using \dQuote{sub
      lists} \code{list(....) = list(value = <vv>, type = <tp>, expr =
      <e>)}, see the details and the examples below.
  }
  \item{dmn}{\bold{named} \code{\link{dimnames}}, a \code{\link{list}}.}
  \item{object}{a \code{"varlist"} object.}
}
\value{
  an object of formal (aka \dQuote{S4}) class \code{"varlist"}.
}
\details{
  \code{value} is typically an atomic vector (\code{\link{is.atomic}})
    or a \code{\link{list}}, e.g., of \code{\link{function}}s; in the
    latter case, typically with \code{\link{names}}.

  \code{type} can be one of \code{"N"}, \code{"frozen"}, \code{"grid"}, or
  \code{"inner"}. In short:
  \describe{
    \item{\code{"N"}}{This type is reserved for a (single) variable
      named \code{n.sim} which provides the simulation replications;
      if it is not given, \code{n.sim} is implicitly treated as 1.}
    \item{\code{"frozen"}}{Variables of this type remain fixed (they do
      not vary) throughout the whole simultion study.  They affect the
      final result but do not appear as a \emph{dimension} in the result
      array of the simulation study.  This is the \bold{default}
      \code{type} (apart from \code{n.sim} which defaults to \code{"N"}).}
    \item{\code{"grid"}}{Variables of this type are used to build a
      (physical) grid (a \code{\link{data.frame}}) with number of rows
      equal to the product of the lengths of all variables of this
      type. The simulation will use this grid to iterate \code{n.sim}
      times over all of its rows for conducting the required
      computations. Conceptually, this corresponds to
      iterating over a \emph{virtual grid} seen
      as \code{n.sim} copies of the (physical) grid pasted together.
      The computations for one row in this virtual grid form one
      sub-job. One can use one of
      \code{\link{doLapply}()}, \code{\link{doForeach}()},
      \code{\link{doRmpi}()}, \code{\link{doMclapply}()}, or
      \code{\link{doClusterApply}()} to iterate over all sub-jobs.}
    \item{\code{"inner"}}{Variables of this type are all dealt with within
      a sub-job for reasons of convenience, speed, load balancing etc.}
  }

  The \code{dimnames2varlist()} functions creates a \code{varlist} from
  a \emph{named} list of \code{\link{character}} vectors, typically
  resulting from \code{\link{dimnames}(tt)} of a \code{\link{table}}
  \code{tt}, see the Titanic example below.

  For more details, see Hofert and \enc{Mächler}{Maechler} (2014), and
  also the examples in \code{ demo(package="simsalapar") }.
}
\author{Martin Maechler.}
\seealso{
  \code{\linkS4class{namedList}}; \code{\link{getEl}} for easy
  extraction of elements from a \code{"varlist"}.

  The \code{\link{toLatex}} method for \code{varlist}s,
  \code{\link{toLatex.varlist}}.

  \code{\link{doLapply}()}, \code{\link{doForeach}()},
  \code{\link{doRmpi}()}, \code{\link{doMclapply}()},
  \code{\link{doClusterApply}()} for the functions to iterate over the
  virtual grid.
}
\examples{
showClass("varlist")

vList <- varlist(
    n.sim = list(value = 1000, expr = quote(N[sim])), # type = N
    n     = list(type="grid", value = c(20, 100, 500)), # sample sizes
    meth  = list(type="grid", expr = quote(italic(method)),
                 value = c("classical", "robust")),
    alpha = list(value = 0.95)) # default type = "frozen"

str(vList)# note the default 'expr' for  n and alpha; and type of alpha

## For more extensive examples, see also
demo(package="simsalapar")

## coerce to simple list  .. and back :
lvl <- as(vList, "list")
stopifnot(identical(
    do.call(varlist, lvl),
    vList ))

## From a data.frame to a LaTeX table :
str(dimnames(Titanic))
vlTitan <- dimnames2varlist(dimnames(Titanic))
vlTitan  # default 'type = "grid"' here
toLatex(vlTitan)
}
\keyword{classes}
\keyword{utilities}
