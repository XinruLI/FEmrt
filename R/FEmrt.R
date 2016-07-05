
require(rpart)
require(rpart.plot)
#' @export
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

    res <- list(tree =  prunedtree, no. = no. ,  Q = Q,
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

#' @export
plot.FEmrt <- function(x, type=4){
  prp(x$tree)
}
