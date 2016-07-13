
#' Prune a tree
#'
#' Prune an initial rpart tree by "c-standard-error" rule.
#' @param tree: A initial tree fitted by rpart, needs to an rpart object.
#' @param c: A scalar to prune the  tree by selecting the tree with minum cross-validation error plus the standard error multiplied by c.
#' @return The pruned tree, also an rpart object.
treepruner <- function(tree, c){
  # function that prunes a CART with c-SE prunign rule
  #
  # Argument:
  #  tree: a classification tree or a regression tree fit by rpart
  #     c: the pruning parameter, needs to be a scalar
  #
  # Returns:
  # a pruned tree
  tree <- tree
  c <- c
  mindex <- which.min(tree$cptable[,4])  # find the row of the minimum x-error
  cp.minse <- tree$cptable[mindex,4] + c*tree$cptable[mindex,5]  # the minimum x-error + c*SE
  cp.row <- min(which(tree$cptable[,4]<= cp.minse))  # find the smallest tree within the minimum x-error + c*SE
  cp.take <- tree$cptable[cp.row, 1]  # get the cp value for the smallest tree
  prune(tree, cp=cp.take)  # prune the tree
}

#' Fixed effects meta-regression tree
#'
#' Apply fixed effects meta-regression trees to meta-analytic data sets.
#' Note that the model is assuming fixed effects within and between subgroups.
#' In other words, it is assumed that no residual heterogeneity is present.
#' @param formula: A formula to specify the response variable (usually the effect size) and the interested moderators. It should be written like y ~ x1 + x2.
#' @param vi: The name of the column of the sampling variance of each study in the data set.
#' @param data: A meta-analytic data set. Needs to be a data frame.
#' @param c: The pruning parameter to prune the initial regression tree by "c-standard-error" rule. Needs to be a non-negative scalar.
#' @return If no moderator effect is detected, the function will perform a standard meta-analysis and return a list including the following objects:
#' @return no.: The total number of the studies
#' @return Q: The Q-statistics for the heterogeneity test
#' @return df: The degree of freedoms of the heterogeneity test
#' @return pval.Q: The p-value of the heterogeneity test
#' @return g: The overall effect size for all studies, based on Hedges'g
#' @return se: The standard error of the overall effect size
#' @return zval: The Z value of the z-test
#' @return pval: The p-value of the significance test for the overall effect size
#' @return ci.lb: The lower bound of the confidence interval for the overall effect size
#' @return ci.ub: The upper bound of the confidence interval for the overall effect size
#' @return formula: The formula that was specified in the model
#' @return If  moderator effect(s) is(are) detected, the function will perform a subgroup meta-analysis and return a list including the following objects:
#' @return tree: The tree that represents the moderator effects and interaction effects between the moderators. An rpart object.
#' @return labs: A data frame showing how the studies were split by the moderators into subgroups
#' @return no.: The number of the studies in each subgroup
#' @return Qb: The between-subgroups Q-statistics
#' @return df: The degree of freedoms of the test for between-subgroups heterogeneity
#' @return pval.Qb: The p-value of the test for between-subgroups heterogeneity
#' @return Qw: The within-subgroup Q-statistic in each subgroup
#' @return g: The overall effect size in each subgroup, based on Hedges'g
#' @return se: The standard error of each overall effect size
#' @return zval: The Z value of the z-test for each overall effect size
#' @return pval: The p-value of the significance test for the overall effect size of each group
#' @return ci.lb: The lower bounds of the confidence intervals
#' @return ci.ub: The upper bounds of the confidence intervals
#' @return formula: The formula that was specified in the model
#' @examples data(SimData)
#' test <- FEmrt(efk~m1+m2+m3+m4+m5, vark, data=SimData, c=0)
#' test
#' plot(test)
#' @export
FEmrt <- function(formula, vi, data, c = 1) {
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
  mf <- match.call()
  mf.vi <- mf[[match("vi", names(mf))]]
  vi <- eval(mf.vi, data, enclos = sys.frame(sys.parent()))
  assign("vi", vi, envir = .GlobalEnv)
  formula <- as.formula(formula)
  data <- as.data.frame(data)
  tree <- rpart(formula, weights = 1/vi, data = data)
  prunedtree <- treepruner(tree, c*sqrt(mean(1/vi)))

  if (length(unique(prunedtree$where)) < 2) {
    warning("no moderator effect was detected")
    y <- model.response(model.frame(formula, data))
    no. <- length(y)
    g <- sum(y/vi)/sum(1/vi)
    Q <- sum((y-g)^2/vi)
    df <- no. - 1
    pval.Q <- pchisq(Q, df, lower.tail = FALSE)
    se <- 1/sqrt(sum(1/vi))
    zval <- g/se
    pval <- pnorm(abs(zval),lower.tail=FALSE)*2
    ci.lb <- g - qnorm(0.975)*se
    ci.ub <- g + qnorm(0.975)*se

    res <- list(no. = no. ,  Q = Q,
                df = df, pval.Q = pval.Q, g = g, se = se, zval = zval,
                pval = pval, ci.lb = ci.lb, ci.ub = ci.ub, formula = formula)

  } else {
    treeframe <- prunedtree$frame
    no. <- treeframe[treeframe$var == "<leaf>", 2]
    Qw <- treeframe[treeframe$var == "<leaf>", 4]
    g <- treeframe[treeframe$var == "<leaf>", 5]
    Qb <- treeframe[1,4] - sum(Qw)
    df <- length(unique(prunedtree$where))-1
    pval.Qb <- pchisq(Qb, df, lower.tail = FALSE)
    se <- tapply(vi, prunedtree$where, function(x) sqrt(1/sum(1/x)))
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
                ci.ub = ci.ub, formula = formula)

  }
  class(res) <- "FEmrt"
  remove(vi, envir = .GlobalEnv)
  res
}

#' Pinrt function for FEmrt
#'
#' Print out the results of an FEmrt object
#'
#' @usage ## S3 method for class "FEmrt"
#' print(x, digits = 3)
#' @details If no moderator effect is detected,
#' the print function will show the standard meta-analysis results.
#' Otherwise, the print function will show the subgroup meta-analysis results,
#' with the significance test resutls for moderator effects, the splitting points of the moderators,
#' and the estimated effect sizes in all subgroups.
#' @export
print.FEmrt <- function(x, digits = 3){
  if (!is.element("FEmrt", class(x))) {
    stop("Argument 'x' must be an object of class \"FEmrt\".")
  } else {
    if (length(x$no.) == 1) {
      cat("\n")
      cat("Fixed effects meta-regression tree (K = ", sum(x$no.), " studies); ",
          sep = "")
      cat("\n")
      cat("formula = ",deparse(x$formula))
      cat("\n")
      cat("No moderator effect was detected" )
      cat("\n\n")
      cat("Test for Heterogeneity:")
      cat("\n")
      Q <- formatC(x$Q, digits=digits, format="f")
      pval <- format.pval(x$pval, eps = 10^(-digits-1))
      cat("Q = ", x$Q," (df = ", x$df, "), ", "p-value ", pval, ";", sep = "")
      cat("\n\n")
      cat("Meta-analysis Results:")
      cat("\n")
      sig <- symnum(x$pval, corr = FALSE, na = FALSE, cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                 symbols = c("***", "**", "*", ".", " "))
      res.table <- c(x$no.,
                     formatC(c(x$g, x$se, x$zval, x$pval, x$ci.lb, x$ci.ub), digits, format = "f"),
                     sig)
      names(res.table) <- c("no.", "g", "se", "zval", "pval", "ci.lb", "ci.ub", " ")
      print(res.table, quote = FALSE, right = TRUE)
      cat("---\nSignif. codes: ", attr(sig, "legend"), "\n\n")

    } else {
      cat("\n")
      cat("Fixed effects meta-regression tree (K = ", sum(x$no.), " studies); ",
          sep = "")
      cat("\n")
      cat("formula = ",deparse(x$formula))
      cat("\n")
      cat("A tree with ", length(x$no.), " terminal nodes was detected", sep="" )
      cat("\n")
      cat("Moderators were detected as: ", paste(colnames(x$labs), collapse = ", "), sep="")
      cat("\n\n")
      cat("Test for Between-Subgroups Heterogeneity:")
      cat("\n")
      Qb <- formatC(x$Qb, digits=digits, format="f")
      pval.Qb <- format.pval(x$pval.Qb, eps = 10^(-digits-1))
      cat("Qb = ", Qb," (df = ", x$df, "), ", "p-value ", pval.Qb, ";", sep = "")
      cat("\n\n")
      cat("Subgroup Meta-analysis Results:")
      cat("\n")
      sig <- symnum(x$pval, corr = FALSE, na = FALSE, cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                    symbols = c("***", "**", "*", ".", " "))
      res.table <- cbind(x$labs, x$no.,
                     formatC(cbind(x$Qw, x$g, x$se, x$zval, x$pval, x$ci.lb, x$ci.ub), digits, format = "f"),
                     sig)
      colnames(res.table) <- c(colnames(x$labs), "K", "Qw", "g", "se", "zval",  "pval",
                               "ci.lb", "ci.ub", " ")
      print(res.table, quote = FALSE, right = TRUE)
      cat("---\nSignif. codes: ", attr(sig, "legend"), "\n\n")
      }

  }
}

#' Visualisation of Fixed Effects Meta-Regression Tree
#'
#' Plot function for an FEmrt object. The plot function uses the plot method from the package rpart.plot
#'  of Stephen Milborrow. The plots shows the result of fixed effects meta-regression tree with subgroups
#'  split by moderators and the corresponding estimated overall effect sizes.
#' @param x: A FEmrt object.
#' @param type:
#' @usage ## S3 method for class "FEmrt"
#' plot(x, type=4)
#' @export
plot.FEmrt <- function(x, type=4){
  prp(x$tree)
}
