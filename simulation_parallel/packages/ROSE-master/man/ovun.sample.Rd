\name{ovun.sample}
\alias{ovun.sample}
\title{
Over-sampling, under-sampling, combination of over- and under-sampling.
}

\description{
Creates possibly balanced samples by random over-sampling minority examples, under-sampling majority examples or combination of over- and under-sampling.
}

\usage{
ovun.sample(formula, data, method="both", N, p=0.5, 
            subset=options("subset")$subset,
            na.action=options("na.action")$na.action, seed)
}

\arguments{
  \item{formula}{
An object of class \code{\link{formula}} (or one that can be coerced to that class). 
	See \code{\link{ROSE}} for information about interaction among predictors or 
	their transformations.
}
  \item{data}{
	An optional data frame, list or environment (or object
	coercible to a data frame by \code{as.data.frame}) in which 
   to preferentially interpret ``formula''. 
   If not specified, the variables are taken from ``environment(formula)''.
}
\item{method}{
  One among \code{c("over", "under", "both")} to perform over-sampling minority examples, under-sampling majority 
  examples or combination of over- and under-sampling, respectively. 
}
  \item{N}{
 The desired sample size of the resulting data set. 
 If missing and \code{method} is either \code{"over"} or \code{"under"} the sample size is determined by oversampling or, respectively, undersampling examples so that the minority class occurs approximately in proportion \code{p}. 
 When \code{method = "both"} the default value is given by the length of vectors 
 specified in \code{formula}. 
}
 \item{p}{
 The probability of resampling from the rare class.
 If missing and \code{method} is either \code{"over"} or \code{"under"} this proportion is determined by oversampling 
 or, respectively, undersampling examples so that the sample size is equal to \code{N}. 
 When \code{method ="both"} the default value given by 0.5.
}
   \item{subset}{
An optional vector specifying a subset of observations to be used in the sampling process.
The default is set by the \code{\link{subset}} setting of \code{\link{options}}.
 } 
 \item{na.action}{
 A function which indicates what should happen when the data contain 'NA's.  
 The default is set by the \code{\link{na.action}} setting of \code{\link{options}}.
}
  \item{seed}{
A single value, interpreted as an integer, recommended to specify seeds and keep trace of the  
  sample.
}
}


\value{
The value is an object of class \code{ovun.sample} which has components
	\item{Call}{The matched call.}
	\item{method}{The method used to balance the sample. Possible choices are \cr \code{c("over", "under", "both")}.}
	\item{data}{ The resulting new data set.} 
}



\seealso{
\code{\link{ROSE}}.
}

\examples{

# 2-dimensional example
# loading data
data(hacide)

# imbalance on training set
table(hacide.train$cls)

# balanced data set with both over and under sampling
data.balanced.ou <- ovun.sample(cls~., data=hacide.train,
                                N=nrow(hacide.train), p=0.5, 
                                seed=1, method="both")$data

table(data.balanced.ou$cls)

# balanced data set with over-sampling
data.balanced.over <- ovun.sample(cls~., data=hacide.train, 
                                  p=0.5, seed=1, 
                                  method="over")$data

table(data.balanced.over$cls)

}

