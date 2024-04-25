\name{grid-stuff}
\alias{getEl}
\alias{mkGrid}
\alias{mkNms}
\alias{get.n.sim}
\alias{set.n.sim}
\alias{get.nonGrids}
\title{Tools For Working with Variable Specification Lists}
\description{
  From a variable specification list (\code{\link{varlist}}),
  \describe{
    \item{getEl()}{gets elements of a variable specification
      list that match the given variable \code{type}.}
    \item{mkGrid()}{builds a grid, e.g., for parallel evaluation,
      basically by calling \code{do.call(\link{expand.grid}, <list>)}.}
    \item{mkNms()}{builds a \code{\link{list}}
      of names, e.g., to be used as \code{\link{dimnames}} for a
      corresponding simulation result \code{\link{array}}.}
    \item{get.n.sim()}{extracts \code{n.sim} or returns 1 if it
      is not contained in \code{vList}.}
    \item{set.n.sim()}{modifies or sets \code{n.sim} in \code{vList}.}
    \item{get.nonGrids()}{extracts all variables not having
      \code{type="grid"} and returns \code{n.sim} the same as
      \code{get.n.sim()}.
    }
  }
}
\usage{
getEl       (vList, type = "ALL", comp = "value")
mkGrid      (vList)
mkNms       (vList, addNms = FALSE)
get.n.sim   (vList)
set.n.sim   (vList, n)
get.nonGrids(vList)
}
\arguments{
  \item{vList}{a \code{\link{list}} of variable specifications,
    typically resulting from \code{\link{varlist}()}.  Each
    variable spec is itself a named list which must contain a
    \code{"value"} component.}
  \item{type}{\code{\link{character}} vector of variable type or types
    to restrict the selection to.  The default, \code{"ALL"} implies no
    restriction and hence returns \emph{all} variables.}
  \item{comp}{either a \code{\link{character}} string containing the
    component name to pick out or \code{\link{NA}} (in which case all
    components are picked out).}
  \item{addNms}{logical, specifying if the resulting names should be of
    the form \code{<var>=<value>} instead of just \code{<value>}.}
  \item{n}{for \code{set.n.sim()}: the value \code{n.sim} should be set
    to; an integer or \code{\link{NULL}}.}
}
\value{
  For \describe{
    \item{getEl()}{a named \code{\link{list}} containing the selected
      components of those variables that match the provided \code{type}.}
    \item{mkGrid()}{a data frame (\code{\link{data.frame}}).}
    \item{mkNms()}{a named \code{\link{list}} of the same
      \code{\link{length}()} and with the same \code{\link{names}()} as
      \code{vList}.}
    \item{get.n.sim()}{\code{n.sim} if it is contained in \code{vList},
      1 otherwise.}
    \item{set.n.sim()}{the \code{\link{varlist}} \code{vList} with a
      \emph{modified} \code{n.sim}.}
    \item{get.nonGrids()}{\code{\link{list}} of length 2 containing the
      (possibly modified) \code{n.sim} and a \code{\link{list}}
      containing all variables not having \code{type="grid"}.}
  }
}
\details{
  These functions are useful when working with variable specification lists.
}
\author{Marius Hofert and Martin Maechler.}
\seealso{
  \code{\link{varlist}}, for construction of variable lists.
  \code{\link{expand.grid}}, the \dQuote{workhorse} of \code{mkGrid()}.
}
\examples{
vList <-
 varlist(n.sim = list(type="N", expr = quote(N[sim]), value = 64),
         n     = list(type="grid",
                      value = c(20, 100, 500)), # sample sizes
         p     = list(type="grid",
                      value = c(3, 7, 15, 25)), # dimensions
         meth  = list(type="grid", expr = quote(italic(method)),
                      value = c("classical", "robust")))

getEl(vList, type="grid") # for those of type "grid", get all values
## for those of type "grid", get all components :
str(getEl(vList, type="grid", comp=NA))
stopifnot(identical(as(vList, "list"),
                    getEl(vList, type=c("N","grid"), comp = NA)))

(grd <- mkGrid(vList))
stopifnot(nrow(grd) == 3*4*2, ncol(grd) == 3)

getEl(vList)# -> all "value"s: the same as lapply(., `[[`, "value") :
stopifnot(identical(lapply(vList, `[[`, "value"),
                    getEl(vList)))
mkNms(vList)
mkNms(vList, addNms=TRUE)

get.n.sim(vl. <- set.n.sim(vList, NULL)) # 1
vl.$n.sim # NULL
set.n.sim(vl., 12)
\dontshow{
stopifnot(is.null(vl.$n.sim),
	  identical(set.n.sim(vl., 12)$n.sim,
		    list(value = 12, type = "N", expr = quote(N[sim]))))
}%dont
}
\keyword{utilities}
