
# Subgroup Discovery
# Copyright (C) 2020  Jurian Baas
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.



#-------------------------------------------------------------------------------------------#
################################## PRIM INTERFACE FUNCTIONS #################################
#-------------------------------------------------------------------------------------------#

#' @title Prepare data for PRIM
#' @description Only numerical and factor data is allowed, this function helps prepare the data by converting
#' any columns that are not numerical or factor to factor
#' @param X The data frame to prepare
#' @return The same data frame, with only numerical and factor data
#' @author Jurian Baas
#' @examples
#' \donttest{
#'   data(credit)
#'   credit <- prim.data.prepare(credit)
#' }
#' @export
prim.data.prepare <- function(X) {
  violating.cols <- !(sapply(X, is.numeric) | sapply(X, is.factor))
  X[violating.cols] <- lapply(X[violating.cols], as.factor)
  return(X)
}


#' @title Prim Box Index
#' @description For a given box (defined as a set of peels) id, return the index of those rows that satisfy all conditions of the box
#' @param object An S3 object of class prim.peel or prim.predict
#' @param newdata A data frame on which to apply the conditions
#' @param box.index Optionally, the index of the peel which defines the box. If not provided, the simplest box with the highest quality is used
#' @return A vector of indices
#' @author Jurian Baas
#' @importFrom stats predict model.frame model.response complete.cases terms formula
#' @examples
#' \donttest{
#'   data(pima)
#'   pima.sample <- sample(nrow(pima), 0.75*nrow(pima))
#'   pima <- prim.data.prepare(pima)
#'   pima.model <- prim(class ~ ., data = pima[pima.sample,],
#'     peeling.quantile = 0.4, min.support = 0.4)
#'   pima.predict <- predict(pima.model, pima[-pima.sample,])
#'   pima.idx <- prim.box.index(pima.predict, pima)
#' }
#' @export
prim.box.index <- function(object, newdata, box.index) {

  if(!(class(object) != "prim.peel" | class(object) != "prim.validate")) {
    stop("Supplied argument is not of class prim.peel or prim.validate, aborting...")
  }

  if(!missing(box.index)) {
    if(box.index > length(object$peels)) {
      stop("Error: box.index is larger than number of peels")
    }
    if(box.index <= 0) {
      stop("Error: box.index must be positive non-zero integer")
    }
  }

  if(class(object) == "prim.peel") {
    X <- stats::model.frame(formula(stats::terms(object$formula, data = newdata, simplify = TRUE)), newdata)
  }else {
    X <- stats::model.frame(formula(stats::terms(object$peel.result$formula, data = newdata, simplify = TRUE)), newdata)
  }
  X <- X[,-1, drop = FALSE]

  # Remove rows with NA values
  test.complete <- stats::complete.cases(X)
  if(!all(test.complete)) {
    warning("Incomplete cases found in data, removing...")
    X <- X[test.complete,]
  }

  # Keep track of which columns are numerical and factors
  col.types <- sapply(X, function(col){
    if(is.numeric(col)) return(0)
    if(is.factor(col)) return(1)
    stop("Data contains invalid data types, only numeric and factors are allowed. Use prim.data.prepare() to make sure the data is properly set up.")
  })

  if(!identical(col.types, object$col.types)) stop("Columns in X differ from those in peeling result")
  if(length(object$peels) == 0) stop("No peeling steps in provided result")

  M <- X

  # Turn factors into numerical
  M[sapply(M, is.factor)] <- lapply(M[sapply(M, is.factor)], function(col) {as.numeric(col)-1})

  if(missing(box.index)) {
    # Exploit the fact that since we use the simplest box with the best quality,
    # so we can make use of the simplified rules (there are fewer of those)
    return(indexCpp(
      object$peels.simplified,
      as.matrix(M),
      length(object$peels.simplified) - 1)
    )

  } else {

    return(indexCpp(
      object$peels,
      as.matrix(M),
      box.index - 1)
    )
  }
}

#' @title Find the optimal sub-box using the PRIM peeling strategy
#' @description By iteratively removing a small portion of the data... Note that categorical data has to be in factor form.
#' @param formula Formula with a response and terms
#' @param data Data frame to find rules in
#' @param peeling.quantile Quantile to peel off for numerical variables
#' @param min.support Minimal size of a box to be valid
#' @return An S3 object of class prim.peel
#' @author Jurian Baas
#' @importFrom stats model.frame model.response complete.cases terms formula
#' @examples
#' \donttest{
#'   data(pima)
#'   pima <- prim.data.prepare(pima)
#'   pima.model <- prim(class ~. , pima, 0.3, 0.4)
#'   plot(pima.model)
#'   summary(pima.model)
#' }
#' @export
prim <- function (
  formula,
  data,
  peeling.quantile = 0.03,
  min.support = 0.05) {

  if(peeling.quantile <= 0) stop("Peeling quantile must be positive")
  if(peeling.quantile >= 1) stop("Peeling quantile must be a fraction smaller than 1")
  if(min.support <= 0) stop("Minimum support must be positive")
  if(min.support >= 1) stop("Minimum support must be a fraction smaller than 1")

  X <- stats::model.frame(formula(stats::terms(formula, data = data, simplify = TRUE)), data)
  y <- stats::model.response(X)
  X <- X[,-1, drop = FALSE]

  if(is.null(y)) stop("Data has no response variable, aborting...")

  # Remove rows with NA values
  test.complete <- stats::complete.cases(X)
  if(!all(test.complete)) {
    warning("Incomplete cases found in data, removing...")
    X <- X[test.complete,]
    y <- y[test.complete]
  }

  # Keep track of which columns are numerical and factors
  col.types <- sapply(X, function(col){
    if(is.numeric(col)) return(0)
    if(is.factor(col)) return(1)
    stop("Data contains invalid data types, only numeric and factors are allowed. Use prim.data.prepare() to make sure the data is properly set up.")
  })

  peel.result <- list()
  class(peel.result) <- "prim.peel"
  peel.result$call <- match.call()
  peel.result$formula <- formula
  peel.result$N <- nrow(X)
  peel.result$peeling.quantile <- peeling.quantile
  peel.result$min.support <- min.support
  peel.result$global.quality <- mean(y)
  peel.result$col.types <- col.types
  peel.result$col.names <- colnames(X)

  M <- X;

  # Turn factors into numerical
  M[sapply(M, is.factor)] <- lapply(M[sapply(M, is.factor)], function(col) {as.numeric(col)-1})

  peel.result$peels <- peelCpp(
    M = as.matrix(M),
    y = y,
    colTypes = col.types,
    alpha = peeling.quantile,
    minSup = min.support
  )

  peel.result$rules <- sapply(peel.result$peels, function(peel){
    paste(
      peel.result$col.names[peel$column + 1],
      switch(peel$type + 1, "\u2265", "\u2264", "!="),
      switch(
        peel$type + 1,
        peel$value,
        peel$value,
        levels(X[,peel$column + 1])[peel$value + 1]
      )
    )
  })

  # Many peels are redundant, we remove them here
  peel.result$peels.simplified <- simplifyCpp (
    peel.result$peels,
    peel.result$col.types,
    which.max(sapply(peel.result$peels, function(peel){peel$quality})) - 1)

  # We generate a string version of the simplified rules here since we need access to the factor levels
  # of the dataset X above
  peel.result$rules.simplified <- sapply(peel.result$peels.simplified , function(peel) {
    paste(
      peel.result$col.names[peel$column + 1],
      switch(peel$type + 1, "\u2265", "\u2264", "!="),
      switch(
        peel$type + 1,
        peel$value,
        peel$value,
        levels(X[,peel$column + 1])[peel$value + 1]
      )
    )
  })

  return(peel.result)
}


#' @title Validate peels
#' @description Validate the results taken from the PRIM peeling process
#' @details This function takes the result of the prim peeling process and applies it to new data. Usually the optimal box in the peeling process is not the best on unobserved data.
#' @param object An S3 object of class prim.peel
#' @param newdata A data frame in which to look for variables with which to predict
#' @param ... further arguments passed to or from other methods
#' @return An S3 object of type prim.predict
#' @author Jurian Baas
#' @importFrom stats predict model.frame model.response complete.cases terms formula
#' @examples
#' \donttest{
#'   data(ames)
#'   ames.sample <- sample(nrow(ames), nrow(ames) * 0.75)
#'   ames <- prim.data.prepare(ames)
#'   ames.model <- prim(SalePrice ~ . - PID - Order, ames[ames.sample,])
#'   ames.predict <- predict(ames.model, ames[-ames.sample,])
#'   plot(ames.model)
#'   plot(ames.predict)
#' }
#' @export
predict.prim.peel <- function(object, newdata, ...) {

  if(class(object) != "prim.peel") {
    stop("Supplied argument is not of class prim.peel, aborting...")
  }

  X <- stats::model.frame(formula(stats::terms(object$formula, data = newdata, simplify = TRUE)), newdata)
  y <- stats::model.response(X)
  X <- X[,-1, drop = FALSE]

  # Remove rows with NA values
  test.complete <- stats::complete.cases(X)
  if(!all(test.complete)) {
    warning("Incomplete cases found in data, removing...")
    X <- X[test.complete,]
    y <- y[test.complete]
  }

  # Keep track of which columns are numerical and factors
  col.types <- sapply(X, function(col){
    if(is.numeric(col)) return(0)
    if(is.factor(col)) return(1)
    stop("Data contains invalid data types, only numeric and factors are allowed. Use prim.data.prepare() to make sure the data is properly set up.")
  })

  if(!identical(col.types, object$col.types)) stop("Columns in X differ from those in peeling result")
  if(length(object$peels) == 0) stop("No peeling steps in provided result")

  predict.result <- list()
  class(predict.result) <- "prim.predict"
  predict.result$call <- match.call()
  predict.result$N <- nrow(X)
  predict.result$global.quality <- mean(y)
  predict.result$col.types <- col.types
  predict.result$col.names <- colnames(X)
  predict.result$peel.result <- object

  M <- X

  # Turn factors into numerical
  M[sapply(M, is.factor)] <- lapply(M[sapply(M, is.factor)], function(col) {as.numeric(col)-1})

  predict.result$peels <- predictCpp(object, as.matrix(M), y)

  # Many peels are redundant, we remove them here
  predict.result$peels.simplified <- simplifyCpp (
    predict.result$peels,
    predict.result$col.types,
    which.max(sapply(predict.result$peels, function(peel){peel$quality})) - 1)

  # We generate a string version of the simplified rules here since we need access to the factor levels
  # of the dataset X above
  predict.result$rules.simplified <- sapply(predict.result$peels.simplified , function(peel) {
    paste(
      predict.result$col.names[peel$column + 1],
      switch(peel$type + 1, "\u2265", "\u2264", "!="),
      switch(
        peel$type + 1,
        peel$value,
        peel$value,
        levels(X[,peel$column + 1])[peel$value + 1]
      )
    )
  })

  return(predict.result)
}

#-------------------------------------------------------------------------------------------#
################################### S3 PLOTTING FUNCTIONS ###################################
#-------------------------------------------------------------------------------------------#

#' @title Plot PRIM peel result
#' @description Plot an S3 object of class prim.peel
#' @param x An S3 object of class prim.peel
#' @param ... Optional arguments to pass on
#' @author Jurian Baas
#' @importFrom graphics par plot points text
#' @return Nothing, this function is called for its side-effects
#' @export
plot.prim.peel <- function(x, ...) {

  box.qualities <- sapply(x$peels, function(peel){peel$quality})
  box.supports <- sapply(x$peels, function(peel){peel$support})
  best.box.idx <- which.max(sapply(x$peels, function(peel){peel$quality}))

  graphics::par(bty = "l")
  graphics::plot (
    c(1, box.supports),
    c(x$global.quality, box.qualities),
    xlim = c(1, 0),
    type = "n",
    xlab = "Support", ylab = "Box quality",
    main = "PRIM peel result")
  graphics::abline(v = x$min.support, col = "red", lwd = 2)
  graphics::lines (
    c(1, box.supports),
    c(x$global.quality, box.qualities))
  graphics::points (
    c(1, box.supports[1:best.box.idx]),
    c(x$global.quality, box.qualities[1:best.box.idx]),
    col= "royalblue4", pch = 19, cex = 1, lty = "solid", lwd = 2)
  graphics::points (
    box.supports[-(1:best.box.idx)],
    box.qualities[-(1:best.box.idx)],
    col= "royalblue2", pch = 4, cex = 1, lty = "solid", lwd = 2)

  print.names <- sapply(1:best.box.idx, function(i) {
    if(i < length(x$peels)) {
      # Print labels up to the best box
      if(x$peels[[i]]$column == x$peels[[i + 1]]$column) {
        # Don't print a label if the next box overwrites this one
        return("")
      } else {
        return(x$rules[i])
      }
    } else {
      # Always print label of the best box
      return(x$rules[i])
    }
  })

  graphics::text (
    box.supports[1:best.box.idx],
    box.qualities[1:best.box.idx] + 0.001,
    labels = print.names,
    cex = 0.7, pos = 2, col = "orangered4", font = 2)
}

#' @title Plot PRIM test result
#' @description Plot an S3 object of class prim.predict
#' @param x An S3 object of class prim.predict
#' @param ... Optional arguments to pass on
#' @author Jurian Baas
#' @importFrom graphics par plot points text
#' @return Nothing, this function is called for its side-effects
#' @export
plot.prim.predict <- function(x, ...) {

  box.qualities <- sapply(x$peels, function(peel){peel$quality})
  box.supports <- sapply(x$peels, function(peel){peel$support})
  best.box.idx <- which.max(sapply(x$peels, function(peel){peel$quality}))

  graphics::par(bty = "l")
  graphics::plot (
    c(1, box.supports),
    c(x$global.quality, box.qualities),
    xlim = c(1, 0),
    type = "n",
    xlab = "Support", ylab = "Box quality",
    main = "PRIM predict result")
  graphics::abline(v = x$peel.result$min.support, col = "red", lwd = 2)
  graphics::lines (
    c(1, box.supports),
    c(x$global.quality, box.qualities))
  graphics::points (
    c(1, box.supports[1:best.box.idx]),
    c(x$global.quality, box.qualities[1:best.box.idx]),
    col= "royalblue4", pch = 19, cex = 1, lty = "solid", lwd = 2)
  graphics::points (
    box.supports[-(1:best.box.idx)],
    box.qualities[-(1:best.box.idx)],
    col= "royalblue2", pch = 4, cex = 1, lty = "solid", lwd = 2)

  print.names <- sapply(1:best.box.idx, function(i) {
    if(i < length(x$peels)) {
      # Print labels up to the best box
      if(x$peels[[i]]$column == x$peels[[i + 1]]$column) {
        # Don't print a label if the next box overwrites this one
        return("")
      } else {
        return(x$peel.result$rules[i])
      }
    } else {
      # Always print label of the best box
      return(x$peel.result$rules[i])
    }
  })

  graphics::text (
    box.supports[1:best.box.idx],
    box.qualities[1:best.box.idx] + 0.001,
    labels = print.names,
    cex = 0.7, pos = 2, col = "orangered4", font = 2)
}


#-------------------------------------------------------------------------------------------#
################################### S3 SUMMARY FUNCTIONS ####################################
#-------------------------------------------------------------------------------------------#

#' @title Summarize a PRIM peeling result object
#' @description Summarize a PRIM peeling result object
#' @param object An S3 object of class prim.peel
#' @param ... Optional arguments to pass on
#' @param round Optional setting to disable rounding
#' @param digits Optional setting to control number of digits to round
#' @author Jurian Baas
#' @return Nothing, this function is called for its side-effects
#' @export
summary.prim.peel <- function(object, ..., round = TRUE, digits = 2) {

  if(!round)  digits = 7

  cat("  ======================================", "\n")
  cat("  ========== PRIM PEEL RESULT ==========", "\n")
  cat("  ======================================", "\n")
  cat("\n")
  best.box.idx <- which.max(sapply(object$peels, function(peel){peel$quality}))

  cat("  Set size: ", object$N, "\n")
  cat("  Set quality: ", round(object$global.quality, digits), "\n")
  cat("\n")
  cat("  ============== BEST BOX ==============", "\n")
  cat("  Box quality: ", round(object$peels[[best.box.idx]]$quality, digits), "\n")
  cat("  Box support: ", round(object$peels[[best.box.idx]]$support, digits), "\n")
  cat("\n")
  cat("  ================ RULES ===============", "\n")
  cat(" ", paste0(object$rules.simplified, collapse = "\n  "))

}

#' @title Summarize a PRIM predict result object
#' @description Summarize a PRIM predict result object
#' @param object An S3 object of class prim.predict
#' @param ... Optional arguments to pass on
#' @param round Optional setting to disable rounding
#' @param digits Optional setting to control number of digits to round
#' @author Jurian Baas
#' @return Nothing, this function is called for its side-effects
#' @export
summary.prim.predict <- function(object, ..., round = TRUE, digits = 2) {

  if(!round)  digits = 7

  cat("  ======================================", "\n")
  cat("  ======== PRIM PREDICT RESULT ========", "\n")
  cat("  ======================================", "\n")
  cat("\n")
  best.box.idx <- which.max(sapply(object$peels, function(peel){peel$quality}))

  cat("  Set size: ", object$N, "\n")
  cat("  Set quality: ", round(object$global.quality, digits), "\n")
  cat("\n")
  cat("  ============== BEST BOX ==============", "\n")
  cat("  Box quality: ", round(object$peels[[best.box.idx]]$quality, digits), "\n")
  cat("  Box support: ", round(object$peels[[best.box.idx]]$support, digits), "\n")
  cat("\n")
  cat("  ================ RULES ===============", "\n")
  cat(" ", paste0(object$rules.simplified, collapse = "\n  "))

}


