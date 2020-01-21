#-------------------------------------------------------------------------------------------#
################################## PRIM INTERFACE FUNCTIONS #################################
#-------------------------------------------------------------------------------------------#

#' @title Prepare data for PRIM
#' @description Only numerical and factor data is allowed, this function helps prepare the data by converting logical and character columns to factor
#' @param X The data frame to prepare
#' @return The same data frame, with only numerical and factor data
#' @author Jurian Baas
#' @export
prim.data.prepare <- function(X) {
  violating.cols <- !(sapply(X, is.numeric) | sapply(X, is.factor))
  X[violating.cols] <- lapply(X[violating.cols], as.factor)
  return(X)
}

#' @title Find the optimal sub-box using the PRIM peeling strategy
#' @description By iteratively removing a small portion of the data... Note that categorical data has to be in factor form.
#' @param formula Formula with a response and terms
#' @param data Data frame to find rules in
#' @param peeling.quantile Quantile to peel off for numerical variables
#' @param min.support Minimal size of a box to be valid
#' @return An S3 object of class prim.peel
#' @author Jurian Baas
#' @importFrom stats model.frame model.response complete.cases terms
#' @export
prim.peel <- function (
  formula,
  data,
  peeling.quantile = 0.03,
  min.support = 0.05) {

  #pima.peel <- prim.peel(formula = class ~., data = pima[train.idx,], peeling.quantile = 0.05, min.support = 0.3)
#pima.validate <- prim.validate(pima.peel, pima[-train.idx, 1:8], pima[-train.idx, 9])

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
  peel.result$X <- X;
  peel.result$y <- y;

  M <- X;

  # Turn factors into numerical
  M[sapply(M, is.factor)] <- lapply(M[sapply(M, is.factor)], function(col) {as.numeric(col)-1})

  peel.result$peels <- peel(
    M = as.matrix(M),
    y = y,
    colTypes = col.types,
    alpha = peeling.quantile,
    minSup = min.support
  )

  peel.result$rules <- sapply(peel.result$peels, function(peel) {
    paste(
      peel.result$col.names[peel$column + 1],
      switch(peel$type + 1, ">", "<", "!="),
      switch(
        peel$type + 1,
        peel$value,
        peel$value,
        levels(X[,i])[peel$value + 1]
      )
    )
  })

  return(peel.result)
}


#' @title Validate peels
#' @description Validate the results taken from the PRIM peeling process
#' @details This function takes the result of the prim peeling process and applies it to new data. Usually the optimal box in the peeling process is not the best on unobserved data.
#' @param peel.result An S3 object of class prim.peel
#' @param data Data frame to find rules in
#' @return An S3 object of type prim.validate
#' @author Jurian Baas
#' @export
prim.validate <- function(peel.result, data) {

  if(class(peel.result) != "prim.peel") {
    stop("Supplied argument is not of class prim.peel, aborting...")
  }

  X <- stats::model.frame(formula(stats::terms(peel.result$formula, data = data, simplify = TRUE)), data)
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

  if(!identical(col.types, peel.result$col.types)) stop("Columns in X differ from those in peeling result")
  if(length(peel.result$peels) == 0) stop("No peeling steps in provided result")

  validation.result <- list()
  class(validation.result) <- "prim.validate"
  validation.result$call <- match.call()
  validation.result$N <- nrow(X)
  validation.result$global.quality <- mean(y)
  validation.result$col.types <- col.types
  validation.result$col.names <- colnames(X)
  validation.result$peel.result <- peel.result
  validation.result$X <- X
  validation.result$y <- y

  M <- X

  # Turn factors into numerical
  M[sapply(M, is.factor)] <- lapply(M[sapply(M, is.factor)], function(col) {as.numeric(col)-1})

  validation.result$peels <- validate(peel.result$peels, as.matrix(M), y)

  return(validation.result)
}

#-------------------------------------------------------------------------------------------#
################################### PRIM PRIVATE FUNCTIONS ##################################
#-------------------------------------------------------------------------------------------#

#' @title Find the optimal box depending on the strategy
#' @description Finds the box with the highest quality or the box closest to the maximum quality minus 2 times the standard error
#' @param prim.validate An object of type "prim.validate"
#' @author Jurian Baas
#' @importFrom utils tail
#' @return The index of the optimal box
prim.box.optimal <- function(prim.validate) {

  if(class(prim.validate) != "prim.validate")
    stop("Argument is not of class prim.validate")

  if(prim.validate$optimal.box == "best") {
    best.box.idx <- which.max(prim.validate$box.qualities)
  } else {
    best.box.idx <- which.max(prim.validate$box.qualities)
    # Find the ordering of boxes leading up and including the best box
    o <- order(prim.validate$box.qualities[1:best.box.idx])
    # The new optimal is 2 standard errors below the maximum
    new.optimal <- utils::tail(prim.validate$box.qualities[o], n = 1) - (2 * prim.validate$metrics$se)
    # Find the box which is closest to this new optimal
    cutoff.point <- o[which.min(abs(prim.validate$box.qualities[o] - new.optimal))]
    # There could be another point close by with a better quality
    # So we pick the optimal in the subset defined by the new best box
    best.box.idx <- which.max(prim.validate$box.qualities[1:cutoff.point])
  }
  if(length(best.box.idx) == 0) stop("Could not find a best box, try adding more features")
  return(best.box.idx)
}

#' @title Calculate statistical metrics
#' @description This function calculates the mean, standard deviation, standard error of the mean, 95% confidence intervals
#' @param prim.validate An object of type "prim validate"
#' @return A list with elements described above
#' @author Jurian Baas
#' @importFrom stats sd
prim.validate.metrics <- function(prim.validate) {

  if(class(prim.validate) != "prim.validate")
    stop("Supplied argument not of class prim.validate")

  # Name it x for clearer code
  x <- prim.validate$box.qualities
  n <- length(x)

  metrics <- list()

  # Calculate sample standard deviation
  ssd <-  sqrt(sum((x - mean(x))^2) / (n - 1))

  # Calculate standard error of the mean
  metrics$se <- ssd / sqrt(n)
  # Calculate standard deviation of sample mean
  metrics$sd <- stats::sd(x) / sqrt(n)


  return(metrics)
}

#-------------------------------------------------------------------------------------------#
################################### S3 PLOTTING FUNCTIONS ###################################
#-------------------------------------------------------------------------------------------#

#' @title Plot PRIM peel result
#' @description Plot an S3 object of class prim.peel
#' @param x An S3 object of class prim.peel
#' @param ... Optional arguments to pass on
#' @return Nothing, this function is called for its side-effects
#' @author Jurian Baas
#' @export
#' @importFrom graphics par plot points text
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
    c(x$global.quality, box.qualities)
  )
  graphics::points (
    c(1, box.supports[1:best.box.idx]),
    c(x$global.quality, box.qualities[1:best.box.idx]),
    col= "royalblue4", pch = 19, cex = 1, lty = "solid", lwd = 2)
  graphics::points (
    box.supports[-(1:best.box.idx)],
    box.qualities[-(1:best.box.idx)],
    col= "royalblue2", pch = 4, cex = 1, lty = "solid", lwd = 2)

  print.names <- character(length = best.box.idx)

  for(i in 1:(best.box.idx)) {
    if(i < length(x$peels) & x$peels[[i]]$column == x$peels[[i + 1]]$column) {
      print.names[i] <- ""
    } else {
      print.names[i] <- x$rules[i]
    }
  }
  print.names[best.box.idx] <- x$rules[best.box.idx]

  graphics::text (
    box.supports[1:best.box.idx],
    box.qualities[1:best.box.idx] + 0.001,
    labels = print.names,
    cex = 0.7, pos = 2, col = "orangered4", font = 2)
}

#' @title Plot PRIM test result
#' @description Plot an S3 object of class prim.validate
#' @param x An S3 object of class prim.validate
#' @param ... Optional arguments to pass on
#' @return Nothing, this function is called for its side-effects
#' @author Jurian Baas
#' @export
#' @importFrom graphics par plot points text
plot.prim.validate <- function(x, ...) {

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
    main = "PRIM validation result")
  graphics::abline(v = x$peel.result$min.support, col = "red", lwd = 2)
  graphics::lines (
    c(1, box.supports),
    c(x$global.quality, box.qualities)
  )
  graphics::points (
    c(1, box.supports[1:best.box.idx]),
    c(x$global.quality, box.qualities[1:best.box.idx]),
    col= "royalblue4", pch = 19, cex = 1, lty = "solid", lwd = 2)
  graphics::points (
    box.supports[-(1:best.box.idx)],
    box.qualities[-(1:best.box.idx)],
    col= "royalblue2", pch = 4, cex = 1, lty = "solid", lwd = 2)

  print.names <- character(length = best.box.idx)

  for(i in 1:(best.box.idx)) {
    if(i < length(x$peels) & x$peels[[i]]$column == x$peels[[i + 1]]$column) {
      print.names[i] <- ""
    } else {
      print.names[i] <- x$peel.result$rules[i]
    }
  }
  print.names[best.box.idx] <- x$peel.result$rules[best.box.idx]

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
  cat(" ", paste0(object$rules[1:best.box.idx], collapse = "\n  "))

}

#' @title Summarize a PRIM test result object
#' @description Summarize a PRIM test result object
#' @param object An S3 object of class prim.validate
#' @param ... Optional arguments to pass on
#' @param round Optional setting to disable rounding
#' @param digits Optional setting to control number of digits to round
#' @author Jurian Baas
#' @return Nothing, this function is called for its side-effects
#' @export
summary.prim.validate <- function(object, ..., round = TRUE, digits = 2) {

  if(!round)  digits = 7

  cat("  ======================================", "\n")
  cat("  ======== PRIM VALIDATE RESULT ========", "\n")
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
  cat(" ", paste0(object$peel.result$rules[1:best.box.idx], collapse = "\n  "))

}


