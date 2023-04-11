\name{sbag}
\alias{sbag}
\title{Implementation of SMOTEBagging}
\usage{
sbag(formula, data, size, alg, smote.k = 5, rf.ntree = 50, svm.ker = "radial")
}
\arguments{
  \item{formula}{A formula specify predictors and target variable. Target variable should be a factor of 0 and 1. Predictors can be either numerical and categorical.}
  \item{data}{A data frame used for training the model, i.e. training set.}
  \item{size}{Ensemble size, i.e. number of weak learners in the ensemble model.}
  \item{alg}{The learning algorithm used to train weak learners in the ensemble model. \emph{cart}, \emph{c50}, \emph{rf}, \emph{nb}, and \emph{svm} are available. Please see Details for more information.}
  \item{smote.k}{Number of k applied in SMOTE algorithm. Default is 5.}
  \item{rf.ntree}{Number of decision trees in each forest of the ensemble model when using \emph{rf} (Random Forest) as base learner. Integer is required.}
  \item{svm.ker}{Specifying kernel function when using svm as base algorithm. Four options are available: \bold{linear}, \bold{polynomial}, \bold{radial}, and \bold{sigmoid}. Default is radial. Equivalent to that in e1071::svm().}
}
\description{
The function implements SMOTEBagging for binary classification. It returns a list of weak learners that are built on training-sets manipulated by SMOTE and random over-sampling. They together consist the ensemble model.
}
\details{
SMOTEBagging uses both SMOTE (Synthetic Minority Over-sampling TEchnique) and random over-sampling to increase minority instances in each bag of Bagging in order to rebalance class distribution. The manipulated training sets contain equal numbers of majority and minority instances, but the proportions of minority instances from SMOTE and random over-sampling vary for different bags, determined by an assigned re-sampling rate \emph{a}. The re-sampling rate \emph{a} is always the multiple of 10, and the function automatically generates a vector of \emph{a}, therefore users do not need to self-define.

The function requires the target varible to be a factor of 0 and 1, where 1 indicates minority while 0 indicates majority instances. Only binary classification is implemented in this version.

Argument \emph{alg} specifies the learning algorithm used to train weak learners within the ensemble model. Totally five algorithms are implemented: \bold{cart} (Classification and Regression Tree), \bold{c50} (C5.0 Decision Tree), \bold{rf} (Random Forest), \bold{nb} (Naive Bayes), and \bold{svm} (Support Vector Machine). When using Random Forest as base learner, the ensemble model is consisted of forests and each forest contains a number of trees.

The object class of returned list is defined as \emph{modelBag}, which can be directly passed to predict() for predicting test instances.
}
\examples{
data("iris")
iris <- iris[1:70, ]
iris$Species <- factor(iris$Species, levels = c("setosa", "versicolor"), labels = c("0", "1"))
model1 <- sbag(Species ~ ., data = iris, size = 10, alg = "c50")
model2 <- sbag(Species ~ ., data = iris, size = 20, alg = "rf", rf.ntree = 100)
model3 <- sbag(Species ~ ., data = iris, size = 40, alg = "svm", svm.ker = "sigmoid")
}
\references{
Wang, S. and Yao, X. 2009. Diversity Analysis on Imbalanced Data Sets by Using Ensemble Models. IEEE Symposium on Computational Intelligence and Data Mining, CIDM '09.

Galar, M., Fernandez, A., Barrenechea, E., Bustince,  H., and Herrera, F. 2012. A Review  on  Ensembles  for  the  Class  Imbalance  Problem:  Bagging-,  Boosting-,  and Hybrid-Based Approaches. IEEE Transactions on Systems, Man, and Cybernetics, Part C (Applications and Reviews). 42(4), pp. 463-484.
}
