\name{ROSE-package}
\alias{ROSE-package}
\alias{ROSEpack}
\docType{package}
\title{
ROSE: Random Over-Sampling Examples}
\description{
	Functions to deal with binary classification
  problems in the presence of imbalanced classes. Synthetic balanced samples are
  generated according to ROSE (Menardi and Torelli, 2014).
  Functions that implement more traditional remedies to the class imbalance
  are also provided, as well as different metrics to evaluate a learner accuracy.
  These are estimated by holdout, bootrstrap or cross-validation methods.
}
\details{
%\tabular{ll}{
%Package: \tab ROSE\cr
%Type: \tab Package\cr
%Version: \tab 0.0-1\cr
%Date: \tab 2013-01-28\cr
%License: \tab GPL 2\cr
%}
  The package pivots on function \code{\link{ROSE}} which generates synthetic balanced
  samples and thus allows to strenghten the subsequent estimation of any binary classifier.
  ROSE (Random Over-Sampling Examples) is a bootstrap-based technique
  which aids the task of binary classification in the presence of rare classes.
  It handles both continuous and categorical data by generating synthetic examples from
  a conditional density estimate of the two classes. %(with the caution of selecting a degenerate kernel to draw synthetic categorical data).
  Different metrics to evaluate a learner accuracy are supplied by
  functions \code{\link{roc.curve}} and \code{\link{accuracy.meas}}.
  Holdout, bootstrap or cross-validation estimators of these accuracy metrics are
  computed by means of ROSE and provided by function \code{\link{ROSE.eval}}, to be used in
  conjuction with virtually any binary classifier.
  Additionally, function \code{\link{ovun.sample}} implements more traditional remedies
  to the class imbalance, such as over-sampling the minority class, under-sampling the majority
  class, or a combination of over- and under- sampling.
 }
\author{
Nicola Lunardon, Giovanna Menardi, Nicola Torelli

Maintainer:
Nicola Lunardon <lunardon@stat.unipd.it>
}
\references{
Lunardon, N., Menardi, G., and Torelli, N. (2014). ROSE: a Package for Binary Imbalanced Learning. \emph{R Jorunal}, 6:82--92.

Menardi, G. and Torelli, N. (2014). Training and assessing classification rules with imbalanced data. \emph{Data Mining and Knowledge Discovery}, 28:92--122.
}
%~~ Optionally other standard keywords, one per line, from file KEYWORDS in ~~
%~~ the R documentation directory ~~
\keyword{ package }
\keyword{ machine learning }
\keyword{ imbalanced data }

\seealso{
\code{\link[nnet:nnet]{nnet}},
\code{\link[rpart:rpart]{rpart}}
}
\examples{
# loading data
data(hacide)

# check imbalance
table(hacide.train$cls)

# train logistic regression on imbalanced data
log.reg.imb <- glm(cls ~ ., data=hacide.train, family=binomial)

# use the trained model to predict test data
pred.log.reg.imb <- predict(log.reg.imb, newdata=hacide.test,
                            type="response")

# generate new balanced data by ROSE
hacide.rose <- ROSE(cls ~ ., data=hacide.train, seed=123)$data

# check (im)balance of new data
table(hacide.rose$cls)

# train logistic regression on balanced data
log.reg.bal <- glm(cls ~ ., data=hacide.rose, family=binomial)

# use the trained model to predict test data
pred.log.reg.bal <- predict(log.reg.bal, newdata=hacide.test,
                            type="response")

# check accuracy of the two learners by measuring auc
roc.curve(hacide.test$cls, pred.log.reg.imb)
roc.curve(hacide.test$cls, pred.log.reg.bal, add.roc=TRUE, col=2)

# determine bootstrap distribution of the AUC of logit models
# trained on ROSE balanced samples
# B has been reduced from 100 to 10 for time saving solely
boot.auc.bal <- ROSE.eval(cls ~ ., data=hacide.train, learner= glm,
                          method.assess = "BOOT",
                          control.learner=list(family=binomial),
                          trace=TRUE, B=10)

summary(boot.auc.bal)
}
