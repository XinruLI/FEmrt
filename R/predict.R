#' Predictions from a fitted metacart object
#'
#' Returns a data frame of predicted effect sizes and moderators from a fitted metacart object
#'
#' @usage ## S3 method for class "metacart"
#' predict(x, newdata)
#' @param x: fitted model object of class "metacart".
#' @param newdata: data frame containing the values at which predictions are required.
#' @return  A data frame containing the predicted effect size and the predictors.
#' @export

predict.metacart <- function(x, newdata){
  if (!inherits(x, "metacart"))
    warning("calling predict.metacart(<fake-metacart-object>) ...")
  mf <- as.formula(x$call$formula)
  tt <- terms(mf)
  ms <- model.frame(delete.response(tt), newdata)
  tree <- x$tree
  if (x$type == "Regression"){
  pred.efk <- predict(tree, newdata, type = "vector")
  } else {
    inx <- match(predict(tree, newdata, type = "prob")[ ,1], predict(tree, type="prob")[, 1])
    node.pred <- tree$where[inx]
    pred.efk <- x$g[as.character(node.pred)]
  }
  data.frame(efk = pred.efk, ms)
}
