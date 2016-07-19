#' Fixed effects meta-classification tree
#'
#' Function to fit fixed effects meta-classification trees to meta-analytic data. Note that
#' The model is assuming fixed effects within and between subgroups.
#' Note that meta-classifcation trees may result in spurious splits, and have lower performance comparing to
#' meta-regression trees. See in the documentation of the metacart-package for more details.
#'
#' @param formula: A formula, with a response variable (usually the effect size) and the interested moderators but no interaction terms.
#' @param vi: The column name of the sampling variance in the data set.
#' @param data: A data frame of a meta-analytic data set, including the effect sizes, sampling variance, and the potential moderators.
#' @param c: A non-negative scalar.The pruning parameter to prune the initial classification tree by "c-standard-error" rule.
#' @return If no moderator effect is detected, the function will return a list containing the following objects:
#' @return no.: The total number of the studies
#' @return Q: The Q-statistics of the heterogeneity test
#' @return df: The degree of freedoms of the heterogeneity test
#' @return pval.Q: The p-value of the heterogeneity test
#' @return g: The overall effect size for all studies, based on Hedges'g
#' @return se: The standard error of the overall effect size
#' @return zval: The test statistic of the overall effect sie
#' @return pval: The p-value for the test statistic of the overall effect size
#' @return ci.lb: The lower bound of the confidence interval for the overall effect size
#' @return ci.ub: The upper bound of the confidence interval for the overall effect size
#' @return call: The matched call
#' @return If  moderator effect(s) is(are) detected, the function will return a list including the following objects:
#' @return tree: An rpart object. The tree that represents the moderator effects and interaction effects between them.
#' @return labs: A data frame showing how the studies were split by the moderators into subgroups
#' @return no.: The number of the studies in each subgroup
#' @return Qb: The between-subgroups Q-statistic
#' @return df: The degree of freedoms of the between-subgroups Q test
#' @return pval.Qb: The p-value of the between-subgroups Q test
#' @return Qw: The within-subgroup Q-statistic in each subgroup
#' @return g: The subgroup overall effect sizes, based on Hedges'g
#' @return se: The standard errors of subgroup overall effect sizes
#' @return zval: The test statistics of the subgroup overall effect sizes
#' @return pval: The p-value for the test statistics of the subgroup overall effect sizes
#' @return ci.lb: The lower bounds of the confidence intervals
#' @return ci.ub: The upper bounds of the confidence intervals
#' @return call: The matched call
#' @return type: The type of the tree
#' @examples data(SimData)
#' test <- FEmct(efk~m1+m2+m3+m4+m5, vark, data=SimData, c=0)
#' test
#' plot(test)
#' @export
FEmct <- function(formula, vi, data, c = 0, dich.value = NULL, control = rpart.control(xval=10,minbucket=5,minsplit=10,cp=0.0001)) {
  # function that applies fixed effects meta-CART algorithms to data set
  #
  # Argument:
  #  formula: the formula to specify the response variable (i.e., effect size) and the moderators
  #       vi: the sampling variance of each study
  #     data: the data set
  #        c: the pruning parameter
  #
  # Returns:
  # a list contains the pruned tree, and the resutls of subgroup meta-analysis
  data <- as.data.frame(data)
  mf <- match.call()
  mf.vi <- mf[[match("vi", names(mf))]]
  vi <- eval(mf.vi, data, enclos = sys.frame(sys.parent()))
  assign("vi.femct.temp", vi, envir = .GlobalEnv)
  formula <- as.formula(formula)
  mf.y <- data[deparse(formula[[2]])][,1]
  if (is.null(dich.value)) {
    data$dich <- (mf.y >= (sum(mf.y/vi.femct.temp)/sum(1/vi.femct.temp)))*1
  } else {
    if (!is.numeric(dich.value)) {stop("dich.value needs to be a scalar")}
    else {data$dich <- ifelse(mf.y >= dich.value, 1, 0)}
  }

  formula[[2]] <- as.name("dich")
  tree <- rpart(formula, weights = (length(vi.femct.temp)/vi.femct.temp)/sum(1/vi.femct.temp), method = "class", data = data, control = control)
  prunedtree <- treepruner(tree, c)

  if (length(unique(prunedtree$where)) < 2) {
    warning("no moderator effect was detected")
    no. <- length(mf.y)
    g <- sum(mf.y/vi.femct.temp)/sum(1/vi.femct.temp)
    Q <- sum((mf.y-g)^2/vi.femct.temp)
    df <- no. - 1
    pval.Q <- pchisq(Q, df, lower.tail = FALSE)
    se <- 1/sqrt(sum(1/vi.femct.temp))
    zval <- g/se
    pval <- pnorm(abs(zval),lower.tail=FALSE)*2
    ci.lb <- g - qnorm(0.975)*se
    ci.ub <- g + qnorm(0.975)*se

    res <- list(no. = no. ,  Q = Q,
                df = df, pval.Q = pval.Q, g = g, se = se, zval = zval,
                pval = pval, ci.lb = ci.lb, ci.ub = ci.ub, call = mf, type="Classification")

  } else {
    no. <- table(prunedtree$where)
    g <- tapply(mf.y/vi.femct.temp,
                prunedtree$where, sum)/tapply(1/vi.femct.temp,
                                              prunedtree$where, sum)
    gt <- sum(mf.y/vi.femct.temp)/sum(1/vi.femct.temp)
    Qt <- sum((mf.y-gt)^2/vi.femct.temp)
    Qw <- tapply((mf.y - g[as.character(prunedtree$where)])^2/vi.femct.temp, prunedtree$where, sum )
    Qb <- Qt - sum(Qw)
    df <- length(unique(prunedtree$where))-1
    pval.Qb <- pchisq(Qb, df, lower.tail = FALSE)
    se <- tapply(vi.femct.temp, prunedtree$where, function(x) sqrt(1/sum(1/x)))
    zval <- g/se
    pval <- pnorm(abs(zval),lower.tail=FALSE)*2
    ci.lb <- g - qnorm(0.975)*se
    ci.ub <- g + qnorm(0.975)*se
    mod.names <- levels(prunedtree$frame$var)[levels(prunedtree$frame$var) != "<leaf>"]
    labs <- sapply(1:length(mod.names), function(x)
      tapply(t(data[mod.names[x]]), prunedtree$where, function(x)(paste(unique(x), collapse = "/")) ))
    colnames(labs) = mod.names  # only for categorical variables
    res <- list(tree =  prunedtree,  labs = labs, no. = no., Qb = Qb, df = df, pval.Qb = pval.Qb,
                Qw = Qw, g = g, se = se, zval =zval, pval = pval, ci.lb = ci.lb,
                ci.ub = ci.ub, call = mf, type = "Classification")

  }
  class(res) <- "metacart"
  remove(vi.femct.temp, envir = .GlobalEnv)
  res
}

