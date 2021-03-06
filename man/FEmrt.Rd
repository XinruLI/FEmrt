% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FEmrt.R
\name{FEmrt}
\alias{FEmrt}
\title{Fixed effects meta-regression tree}
\usage{
FEmrt(formula, vi, data, c = 1, control = rpart.control(xval = 10, minbucket
  = 5, minsplit = 10, cp = 1e-04))
}
\arguments{
\item{formula:}{A formula, with a response variable (usually the effect size) and the interested moderators but no interaction terms.}

\item{vi:}{The column name of the sampling variance in the data set.}

\item{data:}{A data frame of a meta-analytic data set, including the effect sizes, sampling variance, and the potential moderators.}

\item{c:}{A non-negative scalar.The pruning parameter to prune the initial regression tree by "c-standard-error" rule.}

\item{formula:}{A formula, with a response variable (usually the effect size) and the interested moderators but no interaction terms.}

\item{vi:}{The column name of the sampling variance in the data set.}

\item{data:}{A data frame of a meta-analytic data set, including the effect sizes, sampling variance, and the potential moderators.}

\item{c:}{A non-negative scalar.The pruning parameter to prune the initial regression tree by "c-standard-error" rule.}
}
\value{
If no moderator effect is detected, the function will return a list including the following objects:

If no moderator effect is detected, the function will return a list containing the following objects:

no.: The total number of the studies

Q: The Q-statistics of the heterogeneity test

df: The degree of freedoms of the heterogeneity test

pval.Q: The p-value of the heterogeneity test

g: The overall effect size for all studies, based on Hedges'g

se: The standard error of the overall effect size

zval: The test statistic of the overall effect sie

pval: The p-value for the test statistic of the overall effect size

ci.lb: The lower bound of the confidence interval for the overall effect size

ci.ub: The upper bound of the confidence interval for the overall effect size

call: The matched call

If  moderator effect(s) is(are) detected, the function will return a list including the following objects:

tree: An rpart object. The tree that represents the moderator effects and interaction effects between them.

labs: A data frame showing how the studies were split by the moderators into subgroups

no.: The number of the studies in each subgroup

Qb: The between-subgroups Q-statistic

df: The degree of freedoms of the between-subgroups Q test

pval.Qb: The p-value of the between-subgroups Q test

Qw: The within-subgroup Q-statistic in each subgroup

g: The subgroup overall effect sizes, based on Hedges'g

se: The standard errors of subgroup overall effect sizes

zval: The test statistics of the subgroup overall effect sizes

pval: The p-value for the test statistics of the subgroup overall effect sizes

ci.lb: The lower bounds of the confidence intervals

ci.ub: The upper bounds of the confidence intervals

call: The matched call

type: The type of the tree
}
\description{
A function to fit fixed effects meta-regression trees to meta-analytic data.
Note that the model is assuming fixed effects within and between subgroups.
}
\examples{
data(SimData)
test <- FEmrt(efk~m1+m2+m3+m4+m5, vark, data=SimData, c=1)
test
plot(test)
}

