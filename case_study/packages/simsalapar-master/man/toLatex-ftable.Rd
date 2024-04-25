\name{toLatex-ftable}
\title{Convert Flat Contingency Table (ftable) and VarLists to LaTeX Table}
\alias{toLatex.ftable}
\alias{toLatex.varlist}
\alias{fftable}
\alias{tablines}
\alias{cattablines}
\description{
  The \code{\link{ftable}} method of \code{\link{toLatex}()} converts
  an \code{ftable} to a LaTeX table via \code{tablines()}.

  Analogously, the \code{\link{varlist}} method of
  \code{\link{toLatex}()} converts an \code{varlist} to a LaTeX table.

  \code{fftable()} essentially calls \code{format.ftable()} and adds
  attributes \code{ncv} and \code{nrv} to the return object.

  \code{tablines()} computes ingredients for converting a
  \code{\link{character}} \code{\link{matrix}} with attributes to a
  LaTeX table.

  \code{cattablines()} is a small auxiliary function which creates rows
  of a LaTeX table from a given matrix.
}
\usage{
\S3method{toLatex}{ftable}(object, vList = NULL,
        x.escape = FALSE, exprFUN = expr2latex, escapeFUN = escapeLatex,
        align = NULL, booktabs = TRUE, head = NULL,
	rsep = "\\\\\\\\", sp = if(booktabs) 3 else 1.25, rsep.sp = NULL,
	csep = " & ", quote = FALSE, lsep=" \\\\textbar\\\\ ",
	do.table = TRUE, placement = "htbp", center = TRUE,
	fontsize = "normalsize", caption = NULL, label = NULL, \dots)

\S3method{toLatex}{varlist}(object,
	col.vars = c("Variable", "expression", "type", "value"),
	exprFUN = expr2latex, escapeFUN = escapeLatex,
	align = NULL, booktabs = TRUE, head = NULL,
	rsep = "\\\\\\\\", sp = if(booktabs) 3 else 1.25, rsep.sp = NULL, csep = " & ",
	do.table = TRUE, placement = "htbp", center = TRUE,
	fontsize = "normalsize", caption = NULL, label = NULL, ...)

fftable(x, lsep = " | ", quote = FALSE, method = "compact", ...)

tablines(x, align = NULL, booktabs = TRUE, head = NULL,
	 rsep = "\\\\\\\\", sp = if(booktabs) 3 else 1.25, rsep.sp = NULL,
	 csep = " & ", quote = FALSE)

cattablines(x, rsep = "\\\\\\\\", csep = " & ", include.rownames = TRUE)
}
\arguments{
  \item{object}{an \code{\link{ftable}} to be converted to a LaTeX
    table. This is accomplished via \code{\link{format}}ting it.}
  \item{x}{for \code{fftable()} an \code{\link{ftable}} object;
    for \code{tablines()} a \code{\link{character}} \code{\link{matrix}} with
    attributes \code{nrv} and \code{ncv} (as returned by \code{fftable()})
    giving the number of row and column variables, respectively; for
    \code{cattablines()} a \code{\link{numeric}} or
    \code{\link{character}} \code{\link{matrix}}.}
  \item{vList}{a variable specification list see \code{\link{varlist}}.}
  \item{x.escape}{logical indicating if the \dQuote{body} entries of the
    table should be escaped by \code{espaceFUN()}; if false, as by
    default, only the column and row variables are escaped.}
  \item{exprFUN}{a function, by default \code{\link{expr2latex}}, for
    transforming plotmath expressions to equivalent LaTeX strings.}
  \item{escapeFUN}{a function, by default \code{\link{escapeLatex}} which
    \dQuote{escapes} each of its input character strings to valid LaTeX strings.}
  \item{align}{either a \code{\link{character}}
    (e.g., \code{"*{3}{c} S[table-format=1.2]"}) or \code{\link{character}}
    vector (e.g., \code{c("c", "c", "c", "S[table-format=1.2]")}), or
    \code{\link{NULL}} (default).}
  \item{booktabs}{\code{\link{logical}} indicating whether a LaTeX
    table in the format of the LaTeX \code{booktabs} package is
    created (requires the LaTeX \code{booktabs} package loaded in the preamble).}
  \item{head}{either
    \describe{
      \item{\code{\link{character}}}{a \code{\link{vector}} containing
	the lines of the header.}
      \item{\code{\link{NA}}}{do not construct a header.}
      \item{\code{\link{NULL}}}{construct a default header.}
    }
  }
  \item{rsep}{\code{\link{character}} to be inserted at the end of each row.}
  \item{sp}{\code{\link{numeric}} scaling factor for separating
    blocks of rows if \code{rsep.sp} is \code{\link{NULL}}.}
  \item{rsep.sp}{\code{\link{numeric}} of length equal to the number of
    different groups of rows minus one, giving the spaces (interpreted
    as pt) between different groups of rows. If \code{\link{NULL}}, a
    suitable default is constructed.}
  \item{csep}{\code{\link{character}} string for separating different
    cells in a row.}
  \item{quote, lsep, method}{see \code{\link{format.ftable}()} (\R-3.0.0
	 or later).}
  \item{col.vars}{character vector of length 3 or 4 (\code{"expression"}
	 can be omitted), specifying the column names.}
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
  \item{caption}{(if \code{do.table}:) \code{\link{character}} string
    containing the table caption or \code{\link{NULL}} for no caption.}
  \item{label}{(if \code{do.table}:) \code{\link{character}} string
    containing the table label or \code{\link{NULL}} for no label.}
  \item{include.rownames}{\code{\link{logical}} indicating whether row
    names are included in the first column.}
  \item{\dots}{additional arguments passed to \code{\link{format.ftable}()}.}
}
\value{
  \code{toLatex()} returns an object as from \code{\link{wrapLaTable}()}.

  \code{fftable()} returns a formatted flat contingency table as
  returned by \code{format.ftable()} with added attributes \code{ncv}
  (number of column variables) and \code{nrv} (number of row variables).

  \code{tablines()} a list with components
    \item{body}{\code{\link{character}} \code{\link{vector}} of lines of
      the table body.}
    \item{body.raw}{\code{\link{character}} \code{\link{matrix}} of
      cells of the table body.}
    \item{head}{\code{\link{character}} \code{\link{vector}} of lines of
      the table head.}
    \item{head.raw}{\code{\link{character}} \code{\link{matrix}} of
      cells of the table head..}
    \item{align}{alignment string.}
    \item{rsepcol}{\code{\link{character}} \code{\link{vector}}
      containing the row separators (last entries of each row).}

  \code{cattablines()} outputs the formatted lines for copy-and-paste
  into a LaTeX table.
}
\author{Marius Hofert and Martin Maechler.}
\references{see \code{\link{simsalapar-package}}.}
\seealso{
  \code{\link{wrapLaTable}()} for how to wrap the lines of a LaTeX table
  created by \code{tablines()} in a LaTeX table and tabular environment.
}
\examples{
## Different table layouts for the same content
(ft1 <- ftable(Titanic, col.vars = 1:4))
(ft2 <- ftable(Titanic, row.vars = 1))
(ft3 <- ftable(Titanic, row.vars = 1:2))
(ft4 <- ftable(Titanic, row.vars = 1:3))
(ft5 <- ftable(Titanic, row.vars = 1:4))

## What tablines() returns
tablines(fftable(ft2))

## LaTeX (booktabs/non-booktabs) versions
toLatex(ft1, do.table=FALSE)
toLatex(ft1, booktabs=FALSE)
toLatex(ft1, method="col.compact")
toLatex(ft1)
toLatex(ft2)
toLatex(ft3)
toLatex(ft4)
toLatex(ft5, booktabs=FALSE)
toLatex(ft5, method="col.compact")
toLatex(ft5)

## ``poor-man's approach'' for creating lines of a LaTeX table
set.seed(271)
tab <- matrix(runif(6), ncol=3)
ftab <- formatC(tab, digits=4, format="f")
cattablines(ftab)
rownames(ftab) <- LETTERS[1:nrow(ftab)]
cattablines(ftab)
}
\keyword{utilities}
