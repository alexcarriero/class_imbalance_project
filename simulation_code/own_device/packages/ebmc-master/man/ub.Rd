\name{ub}
\alias{ub}
\title{Implementation of UnderBagging}
\usage{
ub(formula, data, size, alg, ir = 1, rf.ntree = 50, svm.ker = "radial")
}
\arguments{
  \item{formula}{A formula specify predictors and target variable. Target variable should be a factor of 0 and 1. Predictors can be either numerical and categorical.}
  \item{data}{A data frame used for training the model, i.e. training set.}
  \item{size}{Ensemble size, i.e. number of weak learners in the ensemble model.}
  \item{alg}{The learning algorithm used to train weak learners in the ensemble model. \emph{cart}, \emph{c50}, \emph{rf}, \emph{nb}, and \emph{svm} are available. Please see Details for more information.}
  \item{ir}{Imbalance ratio. Specifying how many times the under-sampled majority instances are over minority instances. Interger is not required and so such as ir = 1.5 is allowed.}
  \item{rf.ntree}{Number of decision trees in each forest of the ensemble model when using \emph{rf} (Random Forest) as base learner. Integer is required.}
  \item{svm.ker}{Specifying kernel function when using svm as base algorithm. Four options are available: \bold{linear}, \bold{polynomial}, \bold{radial}, and \bold{sigmoid}. Default is radial. Equivalent to that in e1071::svm().}
}
\description{
The function implements UnderBagging for binary classification. It returns a list of weak learners that are built on random under-sampled training-sets. They together consist the ensemble model.
}
\details{
UnderBagging uses random under-sampling to reduce majority instances in each bag of Bagging in order to rebalance class distribution. A 1:1 under-sampling ratio (i.e. equal numbers of majority and minority instances) is set as default.

The function requires the target varible to be a factor of 0 and 1, where 1 indicates minority while 0 indicates majority instances. Only binary classification is implemented in this version.

Argument \emph{alg} specifies the learning algorithm used to train weak learners within the ensemble model. Totally five algorithms are implemented: \bold{cart} (Classification and Regression Tree), \bold{c50} (C5.0 Decision Tree), \bold{rf} (Random Forest), \bold{nb} (Naive Bayes), and \bold{svm} (Support Vector Machine). When using Random Forest as base learner, the ensemble model is consisted of forests and each forest contains a number of trees.

\emph{ir} refers to the intended imbalance ratio of training sets for manipulation. With ir = 1 (default), the numbers of majority and minority instances are equal after class rebalancing. With ir = 2, the number of majority instances is twice of that of minority instances. Interger is not required and so such as ir = 1.5 is allowed.

The object class of returned list is defined as \emph{modelBag}, which can be directly passed to predict() for predicting test instances.
}
\examples{
data("iris")
iris <- iris[1:70, ]
iris$Species <- factor(iris$Species, levels = c("setosa", "versicolor"), labels = c("0", "1"))
model1 <- ub(Species ~ ., data = iris, size = 10, alg = "c50", ir = 1)
model2 <- ub(Species ~ ., data = iris, size = 20, alg = "rf", ir = 1, rf.ntree = 100)
model3 <- ub(Species ~ ., data = iris, size = 40, alg = "svm", ir = 1, svm.ker = "sigmoid")
}
\references{
Barandela,  R.,  Sanchez,  J.,  and  Valdovinos,  R.  2003.  New  Applications  of Ensembles of Classifiers. Pattern Analysis and Applications. 6(3), pp. 245-256.

Galar, M., Fernandez, A., Barrenechea, E., Bustince,  H., and Herrera, F. 2012. A Review  on  Ensembles  for  the  Class  Imbalance  Problem:  Bagging-,  Boosting-,  and Hybrid-Based Approaches. IEEE Transactions on Systems, Man, and Cybernetics, Part C (Applications and Reviews). 42(4), pp. 463-484.
}
