\name{predict.modelBst}
\alias{predict.modelBst}
\title{Predict Method for modelBst Object}
\usage{
  \method{predict}{modelBst}(object, newdata, type = "prob", ...)
}
\arguments{
  \item{object}{A object of \emph{modelBst} class.}
  \item{newdata}{A \emph{data frame} object containing new instances.}
  \item{type}{Types of output, which can be \bold{prob} (probability) and \bold{class} (predicted label). Default is prob.}
  \item{...}{Not used currently.}
}
\description{
Predicting instances in test set using modelBst object
}
\value{
Two type of output can be selected:
 \item{prob}{Estimated probability of being a minority instance (i.e. 1). The probability is averaged by using a majority vote by all weak learners, weighted by error estimation.}
 \item{class}{Predicted class of the instance. Instances of probability larger than 0.5 are predicted as 1, otherwise 0.}
}
\examples{
data("iris")
iris <- iris[1:70, ]
iris$Species <- factor(iris$Species, levels = c("setosa", "versicolor"), labels = c("0", "1"))
samp <- sample(nrow(iris), nrow(iris) * 0.7)
train <- iris[samp, ]
test <- iris[-samp, ]
model <- rus(Species ~ ., data = train, size = 10, alg = "c50") # Build RUSBoost model
prob <- predict(model, newdata = test, type = "prob") # return probability estimation
pred <- predict(model, newdata = test, type = "class") # return predicted class
}


