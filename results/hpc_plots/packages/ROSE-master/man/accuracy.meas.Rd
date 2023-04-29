\name{accuracy.meas}
\alias{accuracy.meas}

\title{
Metrics to evaluate a classifier accuracy in imbalanced learning
}

\description{
This function computes precision, recall and the F measure of a prediction. 
}

\usage{
accuracy.meas(response, predicted, threshold = 0.5) 
}

\arguments{
	\item{response}{
A vector of responses containing two classes to be used to evaluate prediction accuracy. 
It can be of class \code{"factor"}, \code{"numeric"} or \code{"character"}.
}
	\item{predicted}{
A vector containing a prediction for each observation. This can be of class \code{"factor"} or \code{"character"} if the predicted label classes are provided or \code{"numeric"} for the probabilities of the rare class (or a monotonic function of them).
}
	\item{threshold}{
When \code{predicted} is of class \code{numeric}, it defines the probability threshold to classify an example as positive. 
Default value is meant for predicted probabilities and is set to 0.5. See further details below. 
Ignored if \code{predicted} is of class \code{factor}}
}

\details{
Prediction of positive or negative labels depends on the classification threshold, 
here defined as the value such that observations with predicted value greater than the
threshold are assigned to the positive class. Some caution is due in setting the 
threshold as well as in using the default setting both because the default value is meant 
for predicted probabilities and because the default 0.5 is 
not necessarily the optimal 
choice for imbalanced learning. Smaller values set for the threshold correspond to assign a larger 
misclassification costs to the rare class, which is usually the case.

Precision is defined as follows:
\deqn{\frac{\mbox{true positives}}{\mbox{true positives + false positives}}}
Recall is defined as:
\deqn{\frac{\mbox{true positives}}{\mbox{true positives + false negative}}}
The F measure is the harmonic average between precision and recall:
\deqn{2 \cdot \frac{\mbox{precision} \cdot \mbox{recall}}{\mbox{precision+recall}}}

}

\value{
The value is an object of class \code{accuracy.meas} which has components
  \item{Call}{The matched call.}
  \item{threshold}{The selected threshold.}
  \item{precision}{A vector of length one giving the precision of the prediction}
  \item{recall}{A vector of length one giving the recall of the prediction}
  \item{F}{A vector of length one giving the F measure}
}

\references{
Fawcet T. (2006). An introduction to ROC analysis. \emph{Pattern Recognition Letters}, 27 (8), 861--875.
}

\seealso{
\code{\link{roc.curve}} 
}


\examples{
# 2-dimensional example
# loading data
data(hacide)

# imbalance on training set
table(hacide.train$cls)

# model estimation using logistic regression
fit.hacide  <- glm(cls~., data=hacide.train, family="binomial")

# prediction on training set
pred.hacide.train <- predict(fit.hacide, newdata=hacide.train,
                             type="response")

# compute accuracy measures (training set)
accuracy.meas(hacide.train$cls, pred.hacide.train, threshold = 0.02)

# imbalance on test set 
table(hacide.test$cls)

# prediction on test set
pred.hacide.test <- predict(fit.hacide, newdata=hacide.test,
                            type="response")

# compute accuracy measures (test set)
accuracy.meas(hacide.test$cls, pred.hacide.test, threshold = 0.02)

}

\keyword{ supervised classification }

