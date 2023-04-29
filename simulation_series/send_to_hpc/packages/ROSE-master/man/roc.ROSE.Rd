\name{roc.curve}
\alias{roc.curve}
\title{
ROC curve
}

\description{
This function returns the ROC curve and computes the area under the curve (AUC) for binary classifiers.
}

\usage{
roc.curve(response, predicted, plotit = TRUE, add.roc = FALSE, 
          n.thresholds=100, ...)
}

\arguments{
  \item{response}{A vector of responses containing two classes to be used to compute the ROC curve. It can be of class \code{"factor"}, \code{"numeric"} 
or \code{"character"}.}
  \item{predicted}{A vector containing a prediction for each observation. This can be of class \code{"factor"} or \code{"character"} if the predicted label classes are provided or \code{"numeric"} for the probabilities of the rare class (or a monotonic function of them).}
  \item{plotit}{Logical, if \code{TRUE} the ROC curve is plotted in a new window. Default value is set to \code{TRUE}.}
  \item{add.roc}{Logical, if \code{TRUE} the ROC curve is added to an existing window. Default value is set to \code{FALSE}.}
  \item{n.thresholds}{Number of \code{thresholds} at which the ROC curve is computed. 
			Default value is the minimum between 100 and the number of elements in \code{response}. 
	A value of \code{n.thresholds} greater than the length of \code{response} is ignored.}
  \item{\dots}{Further arguments to be passed either to \code{plot} or \code{lines}.}
}

\value{
The value is an object of class \code{roc.curve} which has components
  \item{Call}{The matched call.}
  \item{auc}{The value of the area under the ROC curve.}\cr
  \item{false positive rate}{The false positive rate (or equivalently the complement of sensitivity) of the classifier at the evaluated \code{thresholds}.}
  \item{true positive rate}{The true positive rate (or equivalently the specificity) of the classifier at the evaluated \code{thresholds}.}
  \item{thresholds}{Thresholds at which the ROC curve is evaluated.}
}

\references{
Fawcet T. (2006). An introduction to ROC analysis. \emph{Pattern Recognition Letters}, 27 (8), 861--875.
}


\seealso{
\code{\link{accuracy.meas}}, \code{\link[pROC:roc]{roc}}.

}

\examples{
# 2-dimensional example
# loading data
data(hacide)

# check imbalance on training set
table(hacide.train$cls)

# model estimation using logistic regression
fit.hacide  <- glm(cls~., data=hacide.train, family="binomial")

# prediction on training set
pred.hacide.train <- predict(fit.hacide, newdata=hacide.train)

# plot the ROC curve (training set)
roc.curve(hacide.train$cls, pred.hacide.train, 
          main="ROC curve \n (Half circle depleted data)")

# check imbalance on test set 
table(hacide.test$cls)

# prediction using test set
pred.hacide.test <- predict(fit.hacide, newdata=hacide.test)

# add the ROC curve (test set)
roc.curve(hacide.test$cls, pred.hacide.test, add=TRUE, col=2, 
          lwd=2, lty=2)
legend("topleft", c("Resubstitution estimate", "Holdout estimate"), 
        col=1:2, lty=1:2, lwd=2)
}

\keyword{ supervised classification }

