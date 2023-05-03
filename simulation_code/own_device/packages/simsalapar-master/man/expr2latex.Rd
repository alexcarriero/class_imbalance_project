\name{expr2latex}
\alias{expr2latex}
\alias{escapeLatex}
\title{Translate 'plotmath' expressions to LaTeX}
\description{
  \code{expr2latex()} translates a \dQuote{R graphics annotation}
  expression to the corresponding LaTeX one.

  \code{escapeLatex()}, very similar to its original,
  \code{escape_latex()} from \pkg{fortunes}, escapes certain character
  combinations, such that the result can be used in LaTeX.
}
\usage{
expr2latex(expr)
escapeLatex(x)
}
\arguments{
  \item{expr}{an \R object of \code{\link{class}} \code{expression} or
    \code{language}, typically as from \code{\link{quote}(...)}.}
  \item{x}{a \code{\link{character}} vector.}
}
\details{
  The \code{expr2latex()} function is recursively rendering (sub)
  expressions, until it uses the internal \code{renderAtom()} for simple
  symbols (\code{\link{is.symbol}}).

  We currently work with some tables of math annotation expressions,
  lifted from the corresponding C source of \R itself. (Hidden in
  \pkg{simsalpar}'s namespace, we have \code{AccentTable},
  \code{BinTable}, \code{RelTable}, \code{Lgreek} and \code{Ugreek},
  currently.)

  The current implementation is still incomplete.%, and not fully tested.
}
\value{
  a \code{\link{character}} string with the LaTeX %\eqn{\LaTeX}{LaTeX}
  expression corresponding to \dQuote{R graphics annotation} expression
  \code{expr}.
}
%% \references{
%% }
\author{Martin Maechler.}
\seealso{
  \code{\link{plotmath}} for mathematical expressions to annotate \R
  graphics.

  \code{\link{toLatex}()} and its \code{ftable} method,
  \code{\link[=toLatex.ftable]{toLatex}()}.
}
%% TODO MM:  demo(expr2latex)  should do what we currently have in
%% -------   ../../expr2latex-ex.R
\examples{
expr2latex( quote( N[sim] ) )
expr2latex( quote( N[sim] ~ O(n) ) )
expr2latex( quote(x \%notin\% N) )
expr2latex( quote(x \%+-\% epsilon) )
expr2latex( quote(N[s*m^2]) )
expr2latex( quote( 2^{N[sim] - 3} ~~~ O(n^{n^2}) ) )

escapeLatex(c("#{positives}", "A | B"))
}
\keyword{utilities}
