\name{measure}
\alias{measure}
\title{Calculating Performance Measurement in Class Imbalance Problem}
\usage{
measure(label, probability, metric, threshold = 0.5)
}
\arguments{
  \item{label}{A vector of actual labels of target variable in test set.}
  \item{probability}{A vector of probability estimated by the model.}
  \item{metric}{Measurement used for assessing model performance. \bold{auc}, \bold{gmean}, \bold{tpr}, \bold{tnr}, \bold{f}, and \bold{acc} are available. Please see Details for more information.}
  \item{threshold}{Probability threshold for determining the class of instances. A numerical value ranging from 0 to 1. Default is 0.5}
}
\description{
The function is an interation of multiple performance measurements that can be used to assess model performance in class imbalance problem. Totally six measurements are included.
}
\details{
This function integrates six common measurements. It uses pROC::roc() and pROC::auc() to calculate \bold{auc} (Area Under Curve), while calculates other measurements without dependency on other package: \bold{gmean} (Geometric Mean), \bold{tpr} (True Positive Rate), \bold{tnr} (True Negative Rate),and \bold{f} (F-Measure).

\bold{acc} (Accuracy) is also included for any possible use, although such measurement can be misleading when the classes of test set is highly imbalanced.

\emph{threshold} is the probability cutoff for determing the predicted class of instances. For AUC, users do not need to specify threshold because AUC is not affected by the probability cutoff. However, the threshold is required for other five measurements.
}
\examples{
data("iris")
iris <- iris[1:70, ]
iris$Species <- factor(iris$Species, levels = c("setosa", "versicolor"), labels = c("0", "1"))

# Creat training and test set
samp <- sample(nrow(iris), nrow(iris) * 0.7)
train <- iris[samp, ]
test <- iris[-samp, ]

# Model building and prediction
model <- rus(Species ~ ., data = train, size = 10, alg = "c50")
prob <- predict(model, newdata = test, type = "prob")

# Calculate measurements
auc <- measure(label = test$Species, probability = prob, metric = "auc")
gmean <- measure(label = test$Species, probability = prob, metric = "gmean", threshold = 0.5)
}


