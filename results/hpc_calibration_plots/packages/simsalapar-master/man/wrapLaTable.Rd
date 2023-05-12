\name{wrapLaTable}
\alias{wrapLaTable}
\title{Wrapper for a floating LaTeX Table}
\description{
  \code{wrapLaTable()} wraps (a table and tabular environment) around
  the lines of the body of a LaTeX table and utilizes
  \code{\link{writeLines}()} to write the LaTeX table.
}
\usage{
wrapLaTable(x, align, do.table = TRUE, placement = "htbp", center = TRUE,
            fontsize = "normalsize", booktabs = TRUE,
            caption = NULL, label = NULL)
}
\arguments{
  \item{x}{a \code{\link{character}} \code{\link{vector}} containing the
    lines of the body of the table (for \dQuote{booktabs} tables,
    everything strictly between \code{\\midrule} and
    \code{\\bottomrule}). A table header can be passed via attributes of
    \code{x}.}
  \item{align}{table columns alignment string (e.g.,
    \code{"lcccS[table-format=1.2]"}, the notation of \code{"S[...]"}
    coming from the LaTeX package \code{siunitx}).}
  \item{do.table}{\code{\link{logical}} indicating whether a LaTeX
    \sQuote{table} environment should be used at all.}
  \item{placement}{(if \code{do.table}:) \code{\link{character}} string
    containing a LaTeX table placement string such as \code{"htbp"}.}
  \item{center}{\code{\link{logical}} indicating whether centering
    should happen.}
  \item{fontsize}{\code{\link{character}} string giving a fontsize
    (such as \code{"tiny"}, \code{"scriptsize"}, \code{"footnotesize"},
    \code{"small"}, \code{"normalsize"}, \code{"large"}, \code{"Large"},
    \code{"LARGE"}, \code{"huge"}, or \code{"Huge"}).}
  \item{booktabs}{\code{\link{logical}} indicating whether a LaTeX
    table in the format of the LaTeX \code{booktabs} package is
    created. }
  \item{caption}{(if \code{do.table}:) \code{\link{character}} string
    containing the table caption or \code{\link{NULL}} for no caption.}
  \item{label}{(if \code{do.table}:) \code{\link{character}} string
    containing the table label or \code{\link{NULL}} for no label.}
}
\value{
  a \dQuote{LaTeX table}, of class \code{"Latex"} (where the
  \code{\link{print}} method uses \code{\link{writeLines}()}).
}
\details{
  Note that necessary LaTeX packages (such as \code{tabularx}) have to
  be loaded in the preambel of the corresponding .tex or .Rnw file.
}
\author{Marius Hofert.}
\references{see \code{\link{simsalapar-package}}.}
\seealso{
  \code{\link[=toLatex.ftable]{toLatex}()} where it is used to create
  a LaTeX table.% and tablines() is there, too.
}
\examples{
ftab <- ftable(Titanic, row.vars = 1:2)
fftab <- fftable(ftab)
tlist <- tablines(fftab)

wrapLaTable(structure(tlist$body, head = tlist$head), align = tlist$align,
            caption="The Titanic data set.", label="tab:titanic")
}
\keyword{utilities}
