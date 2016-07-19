
#' Summary function for metacart
#'
#' Summary the results of an metacart object
#'
#' @usage ## S3 method for class "metacart"
#' print(x, digits = 3)
#' @details If no moderator effect is detected,
#' the summary function will show the standard meta-analysis results.
#' Otherwise, the summary function will show the subgroup meta-analysis results,
#' with the significance test resutls for moderator effects, the splitting points of the moderators,
#' and the estimated subgroup overall effect sizes.
#' @export
summary.metacart <- function(x, digits = 3){
  if (!is.element("metacart", class(x))) {
    stop("Argument 'x' must be an object of class \"metacart\".")
  } else {
    if (length(x$no.) == 1) {
      cat("\n")
      cat("Fixed Effects Meta-", x$type , " Tree (K = ", sum(x$no.), " studies); ",
          sep = "")
      cat("\n")
      print(x$call)
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
      cat("Fixed Effects Meta-", x$type , " Tree (K = ", sum(x$no.), " studies); ",
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
#' @usage ## S3 method for class "metacart"
#' plot(x)
#' @export
plot.metacart <- function(x){
  if (length(x$no.) < 2) {stop("no tree was detected")}
  else {prp(x$tree,type=4)}
}

#' Print function for metacart
#'
#' Print the results of an metacart object
#'
#' @usage ## S3 method for class "metacart"
#' print(x)
#' @details
#' The function returns the objects concerning the analysis results.
#'
#' @export
print.metacart <- function(x){
  if (length(x$no.) < 2) {
    print(x)
    cat("\n")
    cat("Fixed Effects Meta-", x$type , " Tree (K = ", sum(x$no.), " studies); ",
        sep = "")
    cat("\n")
    print(x$call)
    cat("\n")
    cat("No moderator effect was detected" )


  } else {
    res <- list(tree = x$tree, labs = x$labs, no. = x$no.,
                Qb = x$Qb, df=x$df, pval.Qb=x$pval.Qb, Qw=x$Qw,
                g=x$g, se=x$se, zval=x$zval, pval=x$pval, ci.lb=x$ci.lb,
                ci.ub= x$ci.ub, call=x$call, type=x$type)
    print(res)
    cat("\n")
    cat("Fixed Effects Meta-", x$type , " Tree (K = ", sum(x$no.), " studies); ",
        sep = "")
    cat("\n")
    cat("formula = ",deparse(x$formula))
    cat("\n")
    cat("A tree with ", length(x$no.), " terminal nodes was detected", sep="" )
    cat("\n")
  }
}

#' @export
print.summary.metacart <- function(x, digits = 3){
  if (!is.element("summary.metacart", class(x))) {
    stop("Argument 'x' must be an object of class \"summary.metacart\".")
  } else {
    if (length(x$no.) == 1) {
      cat("\n")
      cat("Fixed Effects Meta-", x$type , " Tree (K = ", sum(x$no.), " studies); ",
          sep = "")
      cat("\n")
      print(x$call)
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

      print(list( call=x$call, type=x$type, tree = x$tree,
                 Q = x$Q, df=x$df, pval.Q=x$pval.Q,
                 g=x$g, se=x$se, zval=x$zval, pval=x$pval, ci.lb=x$ci.lb,
                 ci.ub= x$ci.ub))


    } else {
      cat("\n")
      cat("Fixed Effects Meta-", x$type , " Tree (K = ", sum(x$no.), " studies); ",
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
      print(list(call=x$call, type=x$type, tree = x$tree, labs = x$labs, no. = x$no.,
                 Qb = x$Qb, df=x$df, pval.Qb=x$pval.Qb, Qw=x$Qw,
                 g=x$g, se=x$se, zval=x$zval, pval=x$pval, ci.lb=x$ci.lb,
                 ci.ub= x$ci.ub,  node=x$node))


  }
}}
