\name{mayplot}
\alias{mayplot}
\title{Matrix-like Plot for Arrays up to Rank 5}
\description{
  Produces a matrix-like plot for arrays up to rank 5,
  using \pkg{grid} and \pkg{gridBase} which allows traditional
  \pkg{graphics}, optionally via a user specified panel function
  \code{panel}.
}
\usage{
mayplot(x, vList, row.vars = NULL, col.vars = NULL,
        xvar, method = if(has.n.sim) "boxplot" else "lines",
       panel.first = NULL, panel.last = NULL,
       type = "l", pch = NULL, ylim = "global",
       log = "", do.legend = TRUE,
       spc = c(0.04/max(1,n.x-1), 0.04/max(1,n.y-1)),
       axlabspc=c(0.12, 0.08), labspc=c(0.04, 0.04),
       n.sim.spc = 0.06, auxcol = c("gray40", "gray78", "gray90", "white"),
       pcol = c("black", "blue", "red", "orange"), grid.lwd = 1.6, ax.lwd = 2,
       tx.cex = 1.2, leg.cex = 1, xlab = NULL, ylab = NA,
       do.n.sim = has.n.sim,
       verbose = getOption("verbose"), show.layout = verbose, \dots)
}
\arguments{
  \item{x}{numeric \code{\link{array}} of \sQuote{rank} 5, i.e.,
    \code{length(dim(x)) == 5}, with \emph{named}
    \code{\link{dimnames}}; typically resulting from a call like
    \code{\link{getArray}(\link{doMclapply}(..))}.}
  \item{vList}{a \code{\link{list}} of variable specifications, see
    \code{\link{varlist}} and \code{\link{mkGrid}}.}

  \item{row.vars}{a dimension name of \code{x}, a string; this variable
    is plotted in the plot rows.}
  \item{col.vars}{a dimension name of \code{x}, a string; this variable
    is plotted in the plot columns.}
  \item{xvar}{dimension name of \code{x}, a string; this variable is
    plotted on the x axis of each sub-plot.}
  \item{method}{\code{\link{character}} string indicating the plot
    method used.  Currently available are \code{"boxplot"} (the default
    if \code{vList} has \code{n.sim}) or \code{"lines"} (otherwise;
    \code{type} adjusts the type of lines used).}
  \item{panel.first, panel.last}{\code{\link{function}} or \code{NULL}
    (default).  If specified, \code{panel.first(x, y, col, ...)} is
    called before and \code{panel.last(x, y, col, ...)} is called after
    the main plotting functions (think \code{boxplot.matrix()} and
    \code{lines()}) are called in each panel.}
  \item{type}{character indicating the type of plotting in the
    non-boxplot case; actually any of the \code{type}s as in
    \code{\link{plot.default}}.}
  \item{pch}{\code{\link{logical}} indicating whether a plot symbol is
    to be used in the legend (default \code{NULL} determines this from
    \code{type}).}
  \item{ylim}{either string \code{"global"}, \code{"local"}, or a
    numeric vector, as for \code{\link{plot.default}}.}

  \item{log}{\code{\link{logical}} indicating if logarithmic scales should be used (in
    the individual plots).}
  \item{do.legend}{\code{\link{logical}} indicating if a legend should be added.}

  \item{spc}{dimensions (x, y) in \dQuote{npc} for the space between sub-plots.
    The default uses a simple adaption to the number of sub-plots in
    each direction.}
  \item{axlabspc}{vector of length two containing the width
    of the y axis label and the height of the x axis label in \dQuote{npc}.}
  \item{labspc}{vector of length two containing the width
    of the box of the row labels and the height of the box of
    the column labels in \dQuote{npc}.}
  \item{n.sim.spc}{space for \code{n.sim} on the bottom right of the
    plot in \dQuote{npc} (only if available).}

  \item{auxcol}{auxiliary colors; vector with four components:
    \enumerate{
      \item color of axes and ticks
      \item background color for the row and column labels
      \item background color for the plots
      \item color of grid lines
    }}
  \item{pcol}{plot base colors.  If more colors than the provided ones
    are required, \code{\link{colorRampPalette}()} is used.}

  \item{grid.lwd}{lwd for grid}
  \item{ax.lwd}{lwd for axes}
  \item{tx.cex}{cex for row and column labels}
  \item{leg.cex}{cex of legend text and \code{n.sim} if appropriate}

  \item{xlab}{x axis label (spanned over all plot columns); when
    \code{NULL}, the default is \code{vList[[xvar]]$expr}; to suppress,
    use \code{NA}.}
  \item{ylab}{y axis label (spanned over all plot rows): Typically a
    label for the \code{"value"} of the simulation.}

  %% \item{ylab.rot}{rotation of the y axis label}
  %% \item{ylas}{the \code{las} argument for the y-axis; the default means
  %%   horizontal tick labels.}
  \item{do.n.sim}{\code{\link{logical}} indicating whether
    \code{n.sim} is displayed on the bottom right of the plot
    (only if available).}
  %% \item{rot.n.sim}{rotation of \code{n.sim} on the bottom right of the
  %%   plot (only if available).}

  \item{verbose}{logical indicating whether more
    information is displayed during plotting.}
  \item{show.layout}{logical indicating whether the grid layout is
    displayed.}

  \item{\dots}{optional arguments passed to \code{panel()}.}
}
\value{
  the layout, invisibly.
}
\author{Marius Hofert and Martin Maechler.}
\seealso{
  \code{\link{matplot}}
  \code{\link[grid]{unit}} and \code{\link[grid]{grid.layout}} from
  package \pkg{grid}.
}
\examples{
vLis <-
 varlist(d = list(type="grid", value = c(10, 100, 1000)),
         family=list(type="grid", value = c("Clayton", "Gumbel")),
         tau = list(type="grid", value = c(0.25, 0.5)),
         alpha = list(type="inner", value = c(0.95, 0.99, 0.999)))

iP <- c(4, 1:3)# <- permutation, putting alpha first
dNms <- mkNms(vLis)[iP]
## an array as from  x <- getArray( doMclapply(vLis, ..) ) :
x <- array(
    c(6.1981, 8.0478, 8.4265, 46.883, 74.359, 86.4394, 432.585, 743.27, 859.35,
      4.8508, 6.0286, 6.3965, 26.380, 35.132, 47.1517, 243.113, 311.36, 342.84,
      7.8546, 8.9769, 9.2199, 78.235, 89.493, 92.2875, 785.674, 893.63, 923.62,
      7.7164, 8.2866, 8.8169, 75.959, 82.806, 88.0626, 756.786, 831.65, 874.70),
    dim = sapply(dNms, length), dimnames = dNms)

mayplot(x, vLis, row.vars="family", col.vars="tau", xvar="alpha", log="y",
       ylab=bquote(widehat(VaR)[alpha]))
## the same, but no xlab and no ylab :
mayplot(x, vLis, row.vars="family", col.vars="tau", xvar="alpha", log="y", xlab=NA)
%% ## and the same with a user specified panel:
%% mayplot(x, vLis, row.vars="family", col.vars="tau", xvar="alpha", log="y",
%%        panel = function(x,y,col,nv) points(x,y, type="l", col=col, lwd=4, lty=6))
\dontshow{
## Should also work, identically:
for(i in 1:5) ## use a random permutation of x's dimensions:
 mayplot(aperm(x, sample(4)),
        vLis, row.vars="family", col.vars="tau", xvar="alpha", log="y")
}% dont
}
\keyword{hplot}
