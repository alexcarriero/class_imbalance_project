\name{sbo}
\alias{sbo}
\title{Implementation of SMOTEBoost}
\usage{
sbo(formula, data, size, alg, over = 100, smote.k = 5, rf.ntree = 50, svm.ker = "radial")
}
\arguments{
  \item{formula}{A formula specify predictors and target variable. Target variable should be a factor of 0 and 1. Predictors can be either numerical and categorical.}
  \item{data}{A data frame used for training the model, i.e. training set.}
  \item{size}{Ensemble size, i.e. number of weak learners in the ensemble model.}
  \item{alg}{The learning algorithm used to train weak learners in the ensemble model. \emph{cart}, \emph{c50}, \emph{rf}, \emph{nb}, and \emph{svm} are available. Please see Details for more information.}
  \item{over}{Specifying over-sampling rate of SMOTE. Only multiple of 100 is acceptable.}
  \item{smote.k}{Number of k applied in SMOTE algorithm. Default is 5.}
  \item{rf.ntree}{Number of decision trees in each forest of the ensemble model when using \emph{rf} (Random Forest) as base learner. Integer is required.}
  \item{svm.ker}{Specifying kernel function when using svm as base algorithm. Four options are available: \bold{linear}, \bold{polynomial}, \bold{radial}, and \bold{sigmoid}. Default is radial. Equivalent to that in e1071::svm().}
}
\description{
The function implements SMOTEBoost for binary classification. It returns a list of weak learners that are built on SMOTE-manipulated training-sets, and a vector of error estimations of each weak learner. The weak learners altogether consist the ensemble model.
}
\details{
Based on AdaBoost.M2, SMOTEBoost uses SMOTE (Synthetic Minority Over-sampling TEchnique) to increase minority instances in each iteration of training weak learners. An over-sampling rate of SMOTE can be defined by users with argument \emph{over}.

The function requires the target varible to be a factor of 0 and 1, where 1 indicates minority while 0 indicates majority instances. Only binary classification is implemented in this version.

Argument \emph{alg} specifies the learning algorithm used to train weak learners within the ensemble model. Totally five algorithms are implemented: \bold{cart} (Classification and Regression Tree), \bold{c50} (C5.0 Decision Tree), \bold{rf} (Random Forest), \bold{nb} (Naive Bayes), and \bold{svm} (Support Vector Machine). When using Random Forest as base learner, the ensemble model is consisted of forests and each forest contains a number of trees.

The object class of returned list is defined as \emph{modelBst}, which can be directly passed to predict() for predicting test instances.
}
\value{
 The function returns a list containing two elements:
 \item{weakLearners}{A list of weak learners.}
 \item{errorEstimation}{Error estimation of each weak learner. Calculated by using (pseudo_loss + smooth) / (1 - pseudo_loss + smooth). \emph{smooth} helps prevent error rate = 0 resulted from perfect classfication during trainging iterations. For more information, please see Schapire et al. (1999) Section 4.2.}
}
\examples{
data("iris")
iris <- iris[1:70, ]
iris$Species <- factor(iris$Species, levels = c("setosa", "versicolor"), labels = c("0", "1"))
model1 <- sbo(Species ~ ., data = iris, size = 10, over = 100, alg = "c50")
model2 <- sbo(Species ~ ., data = iris, size = 20, over = 200, alg = "rf", rf.ntree = 100)
model3 <- sbo(Species ~ ., data = iris, size = 40, over = 300, alg = "svm", svm.ker = "sigmoid")
}
\references{
Chawla, N., Lazarevic, A., Hall, L., and Bowyer, K. 2003. SMOTEBoost: Improving Prediction of the Minority Class in Boosting. In Proceedings European Conference on Principles of Data Mining and Knowledge Discovery. pp. 107-119

Galar, M., Fernandez, A., Barrenechea, E., Bustince,  H., and Herrera, F. 2012. A Review  on  Ensembles  for  the  Class  Imbalance  Problem:  Bagging-,  Boosting-,  and Hybrid-Based Approaches. IEEE Transactions on Systems, Man, and Cybernetics, Part C (Applications and Reviews). 42(4), pp. 463-484.

Freund, Y. and Schapire, R. 1997. A Decision-Theoretic Generalization of On-Line Learning and an Application to Boosting. Journal of Computer and System Sciences. 55, pp. 119-139.

Freund, Y.  and  Schapire,  R.  1996.  Experiments  with  a  new  boosting  algorithm. Machine Learning: In Proceedings of the 13th International Conference. pp. 148-156

Schapire, R. and Singer, Y. 1999. Improved Boosting Algorithms Using Confidence-rated Predictions. Machine Learning. 37(3). pp. 297-336.
}
