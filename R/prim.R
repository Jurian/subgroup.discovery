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
#' @param X Optionally instead of using a formula: Data frame to find rules in
#' @param y Optionally instead of using a formula: Response vector, usually of type numeric
#' @param peeling.quantile Quantile to peel off for numerical variables
#' @param min.support Minimal size of a box to be valid
#' @param parallel Peel columns in parallel, faster for many columns (defaults to true)
#' @return An S3 object of class prim.peel
#' @author Jurian Baas
#' @importFrom stats model.frame model.response complete.cases terms
#' @export
prim.peel <- function (
  formula,
  data,
  X, y,
  peeling.quantile = 0.03,
  min.support = 0.05,
  parallel = T) {

  using.formula <- missing(X) & missing(y)

  # Use the formula interface
  if(using.formula) {
    if(!is.data.frame(data)) stop("Data argument is not a data frame, aborting...")

    X <- stats::model.frame(formula(stats::terms(formula, data = data, simplify = TRUE)), data)
    y <- stats::model.response(X)
    X <- X[,-1, drop = FALSE]

    if(is.null(y)) stop("Data has no response variable, aborting...")

  } else {

    if(!is.data.frame(X)) stop("Parameter X has to be a data frame")
    if(!is.vector(y)) stop("Parameter y has to be a vector")
    if(nrow(X) != length(y)) stop("Parameters X and y are not of same size")

    # Remove rows with NA values
    test.complete <- stats::complete.cases(X)
    if(!all(test.complete)) {
      warning("Incomplete cases found in data, removing...")
      X <- X[test.complete,]
      y <- y[test.complete]
    }
  }

  if(peeling.quantile <= 0) stop("Peeling quantile must be positive")
  if(peeling.quantile >= 1) stop("Peeling quantile must be a fraction smaller than 1")
  if(min.support <= 0) stop("Minimum support must be positive")
  if(min.support >= 1) stop("Minimum support must be a fraction smaller than 1")

  # Keep track of which columns are numerical and factors
  col.types <- sapply(X, function(col){
    if(is.numeric(col)) return(0)
    if(is.factor(col)) return(1)
    stop("Data contains invalid data types, only numeric and factors are allowed. Use prim.data.prepare() to make sure the data is properly set up.")
  })

  # Turn factors into numerical
  X[sapply(X, is.factor)] <- lapply(X[sapply(X, is.factor)], function(col) {as.numeric(col)-1})

  result <- list()
  class(result) <- "prim.peel"

  result$N <- nrow(X)
  result$peeling.quantile <- peeling.quantile
  result$min.support <- min.support
  result$global.quality <- mean(y)
  result$col.types <- col.types

  result$peels <- peel(
    M = as.matrix(X),
    y = y,
    colTypes = col.types,
    alpha = peeling.quantile,
    minSup = min.support,
    parallel = parallel
  )

  if(using.formula) result$formula <- formula

  return(result)
}


#' @title Validate peels
#' @description Validate the results taken from the PRIM peeling process
#' @details This function takes the result of the prim peeling process and applies it to new data. Usually the optimal box in the peeling process is not the best on unobserved data.
#' @param peel.result An S3 object of class prim.peel
#' @param X A data frame with at least those columns that were used in creating the prim.peel S3 object
#' @param y Response vector, usually of type numeric
#' @return An S3 object of type prim.validate
#' @author Jurian Baas
#' @export
prim.validate <- function(peel.result, X, y) {

  if(class(peel.result) != "prim.peel") {
    stop("Supplied argument is not of class prim.peel, aborting...")
  }

  # Keep track of which columns are numerical and factors
  col.types <- sapply(X, function(col){
    if(is.numeric(col)) return(0)
    if(is.factor(col)) return(1)
    stop("Data contains invalid data types, only numeric and factors are allowed. Use prim.data.prepare() to make sure the data is properly set up.")
  })

  if(!identical(col.types, peel.result$col.types)) stop("Columns in X differ from those in peeling result")
  if(length(prim.result$peels) == 0) stop("No peeling steps in provided result")
  if(nrow(X) != length(y)) stop(paste("Parameters X and y are not of same size:", nrow(X), length(y)))



  result$call <- match.call()

  return(result)
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

#' @title Create box matching index
#' @description Generate a logical vector in which elements are true iff applying the rule to a record evaluates to true.
#' @param prim.object An S3 object of class prim.peel or prim.validate result
#' @param X A data frame with at least those columns that were used in creating the prim S3 object
#' @return A logical index of matching records
#' @author Jurian Baas
prim.rule.match <- function(prim.object, X) {

  supported.classes <- c("prim.peel", "prim.validate")

  if(!class(prim.object) %in% supported.classes) {
    stop("Supplied argument is not of class prim.peel or prim.validate, aborting...")
  }

  return(with(X, eval(parse(text = paste0(prim.object$superrule, collapse = " & ")))))
}

#' @title Condense multiple (redundant) rules
#' @details This function condenses the many (redundant) rules of an S3 object of class prim.peel or prim.validate to a single rule.
#' @description This function condenses the many (redundant) rules of an S3 object of class prim.peel or prim.validate to a single rule.
#' @param prim.object An S3 object of class prim.validate
#' @return The condensed rule as a single string
#' @author Jurian Baas
#' @importFrom stats aggregate
prim.rule.condense <- function(prim.object) {

  if(class(prim.object) != "prim.validate") stop("Supplied argument is not of class prim.validate")

  # Search for the rules leading up to the best box
  rule.idx <- 1:prim.object$best.box.idx

  # Combine the lists into a single data frame
  rules <- data.frame (
    name = factor(prim.object$rule.names[rule.idx]),
    operator = factor(prim.object$rule.operators[rule.idx]),
    value = as.character(prim.object$rule.values[rule.idx]),
    #quality = prim.object$box.qualities[rule.idx],
    type = prim.object$rule.types[rule.idx]
  )

  numerical.rules.idx <- rules$type == "numeric"

  if(sum(numerical.rules.idx) > 0) {

    numerical.rules <- rules[numerical.rules.idx,1:3]

    numerical.rules <- stats::aggregate (
      list(value = numerical.rules$value),
      list(name = numerical.rules$name, operator = numerical.rules$operator),
      function(x) rev(x)[1])

    other.rules <- rules[!numerical.rules.idx,1:3]
    condensed.rules <- rbind(numerical.rules, other.rules)

  } else {

    condensed.rules <- rules[, 1:3]

  }

  return(sort(apply(condensed.rules, 1, paste, collapse = " ")))
}

#' @title Intersection of multiple rules
#' @description This function applies the rules given by the parameters to a dataset and calculates the intersection, i.e. those observations where all rules evaluate to TRUE are returned as TRUE.
#' @param X Data frame to apply rules to
#' @param prim.objects A list of objects of class "prim.peel" and/or "prim.validate"
#' @param operation One of {"union", "intersect"}
#' @author Jurian Baas
#' @return Logical vector, true iff all rules evaluate to TRUE for a certain observation
prim.rule.operations <- function(X, prim.objects, operation = c("union", "intersect")) {

  supported.classes <- c("prim.peel", "prim.validate")
  used.classes <- sapply(prim.objects, class)

  if(!all(used.classes %in% supported.classes)) {
    stop("Supplied prim.object is not of class prim.peel or prim.validate")
  }

  # Match the dataset for each rule
  matches <- sapply(X = prim.objects, FUN = prim.rule.match, X)

  if(operation == "union") {
    # Apply logical OR row wise
    matches <- apply(matches, 1, any)
  }else {
    # Apply logical AND row wise
    matches <- apply(matches, 1, all)
  }

  return(matches)
}

#' @title Compare PRIM diversify results
#' @description Compares all attempts of a PRIM diversify operation with each other
#' @details Comparison is done by the following formula \eqn{\frac{|A \cap B|}{|A \cup B|}}
#' @param X Data frame to do intersect and union operations
#' @param p.div An S3 object of type "prim.diversify"
#' @author Jurian Baas
#' @return A matrix with the comparisons laid out by position
#' @importFrom utils combn
prim.diversify.compare <- function(X, p.div) {

  if(class(p.div) != "prim.diversify")
    stop("Argument is not of class prim.diversify")

  frontier <- rev(p.div$frontier)
  nr.of.attempts <- length(frontier)
  idx <- combn(1:nr.of.attempts, 2)
  scores <- apply(idx, 2, function(i) {
    sum(prim.rule.operations(X, p.div$attempts[frontier[i]], "intersect")) / sum(prim.rule.operations(X, p.div$attempts[frontier[i]], "union"))
  })
  m <- matrix(ncol = nr.of.attempts, nrow = nr.of.attempts)
  m[t(idx)] <- scores
  m[t(idx[ nrow(idx):1, ])] <- scores
  colnames(m) <- frontier
  rownames(m) <- frontier
  return(m)
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
################################### S3 PREDICT FUNCTIONS ####################################
#-------------------------------------------------------------------------------------------#

#' @title Predict method for PRIM Cover Fits
#' @description Predicted values based on the PRIM cover object
#' @param object Object of class prim.cover
#' @param ...  Further arguments passed to or from other methods
#' @param newdata A data frame in which to look for variables with which to predict.
#' @return Depends on the quality function used. In the case of base::mean, the mean of the target variable of the first matching box.
#' @author Jurian Baas
#' @export
predict.prim.cover <- function(object, newdata, ... ) {

  i <- 1

  y.hat <- as.numeric(rep(NA, nrow(newdata)))

  repeat {

    if(all(!is.na(y.hat))) break
    if(i > length(object$covers)) break

    cover <- object$covers[[i]]

    idx <- prim.rule.match(cover, newdata)

    y.hat[idx & is.na(y.hat)] <- cover$cover.box.quality

    i <- i + 1
  }

  y.hat[is.na(y.hat)] <- object$leftover$cover.global.quality

  return(y.hat)
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

  best.box.idx <- which.max(x$box.qualities)

  graphics::par(bty = "l")
  graphics::plot (
    c(1, x$supports),
    c(x$global.quality, x$box.qualities),
    xlim = c(1, 0),
    type = "n",
    xlab = "Support", ylab = "Box quality",
    main = "PRIM peel result",
    ...)
  graphics::abline(v = x$min.support, col = "red", lwd = 2)
  graphics::lines (
    c(1, x$supports),
    c(x$global.quality, x$box.qualities)
  )
  graphics::points (
    c(1, x$supports[1:best.box.idx]),
    c(x$global.quality, x$box.qualities[1:best.box.idx]),
    col= "royalblue4", pch = 19, cex = 1, lty = "solid", lwd = 2)
  graphics::points (
    x$supports[-(1:best.box.idx)],
    x$box.qualities[-(1:best.box.idx)],
    col= "royalblue2", pch = 4, cex = 1, lty = "solid", lwd = 2)

  print.names <- character(length = best.box.idx)

  for(i in 1:(best.box.idx)) {
    if(i < length(x$rule.names) & x$rule.names[i] == x$rule.names[i + 1]) {
      print.names[i] <- ""
    } else {
      print.names[i] <- paste(x$rule.names[i], x$rule.operators[i], x$rule.values[i])
    }
  }
  print.names[best.box.idx] <- paste(x$rule.names[best.box.idx], x$rule.operators[best.box.idx], x$rule.values[best.box.idx])

  graphics::text (
    c(1, x$supports[1:best.box.idx]),
    c(x$global.quality, x$box.qualities[1:best.box.idx]),
    labels = c("", print.names),
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

  best.box.idx <- x$best.box.idx

  graphics::par(bty = "l")

  graphics::plot (
    c(1, x$supports),
    c(x$global.quality, x$box.qualities),
    xlim = c(1, 0),
    type = "n",
    xlab = "Support", ylab = "Box quality",
    main = "PRIM validate result",
    ...)

  if(x$optimal.box == "2se")
    graphics::abline(h = max(x$box.qualities) - 2 * x$metrics$se, lty = 2)
  graphics::abline(v = x$min.support, col = "red", lwd = 2)
  graphics::lines (
    c(1, x$supports),
    c(x$global.quality, x$box.qualities)
  )
  graphics::points (
    c(1, x$supports[1:best.box.idx]),
    c(x$global.quality, x$box.qualities[1:best.box.idx]),
    col= "royalblue4", pch = 19, cex = 1, lty = "solid", lwd = 2)
  graphics::points (
    x$supports[-(1:best.box.idx)],
    x$box.qualities[-(1:best.box.idx)],
    col= "royalblue2", pch = 4, cex = 1, lty = "solid", lwd = 2)

  print.names <- character(length = best.box.idx)

  if(best.box.idx > 1)
    for(i in 1:(best.box.idx-1)) {
      if(i < length(x$rule.names) & x$rule.names[i] == x$rule.names[i + 1]) {
        print.names[i] <- ""
      } else {
        print.names[i] <- paste(x$rule.names[i], x$rule.operators[i], x$rule.values[i])
      }
    }
  print.names[best.box.idx] <- paste(x$rule.names[best.box.idx], x$rule.operators[best.box.idx], x$rule.values[best.box.idx])

  graphics::text (
    c(1, x$supports[1:best.box.idx]),
    c(x$global.quality, x$box.qualities[1:best.box.idx]),
    labels = c("", print.names),
    cex = 0.7, pos = 2, col = "orangered4", font = 2)
}

#' @title Plot PRIM cover result
#' @description Plot an S3 object of class prim.cover
#' @param x An S3 object of class prim.cover
#' @param ... Optional arguments to pass on
#' @return Nothing, this function is called for its side-effects
#' @author Jurian Baas
#' @export
#' @importFrom graphics par plot lines points text
plot.prim.cover <- function(x, ...) {

  dat <- data.frame(t(sapply(x$covers, function(c) {
    c( supports = c$cover.box.N / c$cover.N,
       box.qualities = c$cover.box.quality)
  })))
  dat <- rbind(dat, c(
    supports = x$leftover$cover.N / x$covers[[length(x$covers)]]$cover.N,
    box.qualities = x$leftover$cover.global.quality
  ))

  graphics::par(bty = "l")
  graphics::plot (
    dat$supports,
    dat$box.qualities,
    type = "n",
    xlim = c(0,1),
    xlab = "Relative support", ylab = "Relative box quality",
    main = "PRIM cover result",
    ...)
  graphics::abline(v = x$min.support, col = "red", lwd = 2)
  graphics::lines (
    x = dat$supports,
    y = dat$box.qualities,
    lwd = 2,
    col = "ivory4"
  )
  graphics::points (
    dat$supports,
    dat$box.qualities,
    col= "royalblue4", pch = 19, cex = 1, lty = "solid", lwd = 2)
  graphics::text (
    dat$supports,
    dat$box.qualities,
    labels = c(paste("Cover", 1 : length(x$covers)), "Leftover"),
    cex = 0.7, pos = 2, col = "orangered4", font = 2)
}

#' @title Plot PRIM diversify result
#' @description Plot an S3 object of class prim.diversify
#' @param x An S3 object of class prim.diversify
#' @param ... Optional arguments to pass on
#' @return Nothing, this function is called for its side-effects
#' @author Jurian Baas
#' @export
#' @importFrom graphics par plot points text
plot.prim.diversify <- function(x, ...) {

  dat <- data.frame(t(sapply(x$attempts, function(a) {
    c( supports = a$final.box.N / x$N,
       box.qualities = a$final.box.quality)
  })))

  frontier.index <- x$frontier
  frontier <- dat[frontier.index,]
  dominated <- dat[-frontier.index,]

  graphics::par(bty = "l")
  graphics::plot (
    dat$supports,
    dat$box.qualities,
    type = "n",
    xlab = "Support", ylab = "Box quality",
    main = "PRIM diversify result",
    ...)
  graphics::lines (
    x = frontier$supports,
    y = frontier$box.qualities,
    lwd = 2,
    col = "ivory4"
  )
  graphics::points (
    dominated$supports,
    dominated$box.qualities,
    col= "royalblue2", pch = 4, cex = 1, lty = "solid", lwd = 2)
  graphics::points (
    frontier$supports,
    frontier$box.qualities,
    col= "royalblue4", pch = 19, cex = 1.25, lty = "solid", lwd = 2)
  graphics::text (
    frontier$supports,
    frontier$box.qualities,
    labels = frontier.index,
    cex = 0.85, adj = -c(0.5,0.5), col = "orangered4", font = 2)


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
  best.box.idx <- which.max(object$box.qualities)

  cat("  Original Set size: ", object$N, "\n")
  cat("  Original Set quality: ", round(object$global.quality, digits), "\n")
  cat("\n")
  cat("  ============== BEST BOX ==============", "\n")
  cat("  Box quality: ", round(object$box.qualities[best.box.idx], digits), "(", round(object$box.qualities[best.box.idx] / object$global.quality, digits), ") \n")
  cat("  Box support: ", round(object$supports[best.box.idx], digits), " (", object$supports[best.box.idx] * object$N, ") \n")
  cat("\n")
  cat("  ================ RULES ===============", "\n")
  cat(" ", paste0(object$superrule, collapse = "\n  "))

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
  best.box.idx <- object$best.box.idx

  cat("  Original set size: ", object$N, "\n")
  cat("  Original set quality: ", round(object$global.quality, digits), "\n")
  cat("\n")
  cat("  ============== BEST BOX ==============", "\n")
  cat("  Box quality: ", round(object$box.qualities[best.box.idx], digits), "(", round(object$box.qualities[best.box.idx] / object$global.quality, digits), ") \n")
  cat("  Box support: ", round(object$supports[best.box.idx], digits), " (", object$supports[best.box.idx] * object$N, ") \n")
  cat("\n")
  cat("  ================ RULES ===============", "\n")
  cat(" ", paste0(object$superrule, collapse = "\n  "))

}

#' @title Summarize a PRIM cover result object
#' @description Summarize a PRIM cover result object
#' @param object An S3 object of class prim.cover
#' @param ... Optional arguments to pass on
#' @param round Optional setting to disable rounding
#' @param digits Optional setting to control number of digits to round
#' @author Jurian Baas
#' @return Nothing, this function is called for its side-effects
#' @export
summary.prim.cover <- function(object, ..., round = TRUE, digits = 2) {

  if(!round)  digits = 7

  cat("  ======================================", "\n")
  cat("  ========== PRIM COVER RESULT =========", "\n")
  cat("  ======================================", "\n")
  cat("  |\n")
  cat("  |  Peeling quantile:\t", object$peeling.quantile, "\n")
  cat("  |  Min support:\t", object$min.support, "\n")
  cat("  |  Train/test split:\t", object$train.fraction, "\n")
  cat("  |  Quality function:\t", object$quality.function.name, "\n")
  cat("\n")

  for(i in 1:length(object$covers)) {
    x <- object$covers[[i]]
    cat("\n")
    cat("  ======================================", "\n")
    cat("  ============== COVER", i,"===============", "\n")
    cat("  |  Cover size:\t", x$cover.N, "\n")
    cat("  |  Cover quality:\t", round(x$cover.global.quality, digits), "\n")
    cat("  |\n")
    cat("  |  Box quality:\t", round(x$cover.box.quality, digits), "\n")
    cat("  |  Box support:\t", round(x$cover.box.N / x$cover.N, digits), "\n")
    cat("  |  Box size:\t\t", x$cover.box.N, "\n")
    cat("\n")
    cat("  ================ RULES ===============", "\n")
    cat("  | ", paste0(x$superrule, collapse = "\n  |  "))
    cat("\n","\n","\n")

  }

  x <- object$leftover

  cat("\n")
  cat("  ======================================", "\n")
  cat("  ============== LEFTOVER ==============", "\n")
  cat("  |  Cover set size: ", x$cover.N, "\n")
  cat("  |  Cover set quality: ", round(x$cover.global.quality, digits), "\n")
  cat("\n","\n","\n")

}

#' @title Summarize a PRIM diversify object
#' @description Summarize a PRIM diversify result object
#' @param object An S3 object of class prim.diversify
#' @param ... Optional arguments to pass on
#' @param round Optional setting to disable rounding
#' @param digits Optional setting to control number of digits to round
#' @author Jurian Baas
#' @return Nothing, this function is called for its side-effects
#' @export
summary.prim.diversify <- function(object, ..., round = TRUE, digits = 2) {

  if(!round) digits = 7

  cat("  ======================================", "\n")
  cat("  ======== PRIM DIVERSIFY RESULT =======", "\n")
  cat("  ======================================", "\n")
  cat("  |\n")
  cat("  |  Peeling quantile:\t", object$peeling.quantile, "\n")
  cat("  |  Min support:\t", object$min.support, "\n")
  cat("  |  Train/test split:\t", object$train.fraction, "\n")
  cat("  |  Quality function:\t", object$quality.function.name, "\n")
  cat("  |\n")
  cat("  |  Set size:\t\t", object$N, "\n")
  cat("  |  Set quality:\t", round(object$global.quality, digits), "\n")
  cat("\n")
  cat("  Scores:", "\n  | ")
  cat(paste0(
    base::gsub(
      "NA",
      "    ",
      apply(formatC(round(object$scoreMatrix, 2), format = "f", digits = 2), 1, paste, collapse = "\t"))
    , collapse = "\n  |  ")
    )
  cat("\n\n\n")
  cat("  Dominating attempts:","\n")
  frontier <- base::sort(object$frontier)
  for(i in base::seq_along(frontier)) {
    x <- object$attempts[[frontier[i]]]
    cat("\n")
    cat("  ======================================", "\n")
    cat("  ============= ATTEMPT", frontier[i],"=============", "\n")
    cat("  |  Score:\t\t", round(object$scores[i], digits), "\n")
    cat("  |  Box quality:\t", round(x$final.box.quality, digits), "\n")
    cat("  |  Box support:\t", round(x$final.box.N / object$N, digits), "\n")
    cat("  |  Box size:\t\t", x$final.box.N, "\n")
    cat("\n")
    cat("  ================ RULES ===============", "\n")
    cat("  | ", paste0(x$superrule, collapse = "\n  |  "))
    cat("\n","\n","\n")

  }
}

#-------------------------------------------------------------------------------------------#
################################### MISC HELPER FUNCTIONS ###################################
#-------------------------------------------------------------------------------------------#

#' @title Calculate a frontier of dominating points
#' @description During the diversify process, we are really only interested in the attempts which dominate all others in performance.
#' @param p.div An object of type "prim.diversify"
#' @author William Huber & Jurian Baas
#' @return A vector of indexes for the dominating points
#' @seealso \url{https://stats.stackexchange.com/a/65157}
quasi.convex.hull <- function(p.div) {

  if(class(p.div) != "prim.diversify")
    stop("Parameter not of class prim.diversify")

  X <- data.frame(t(sapply(p.div$attempts, function(a) {
    c( supports = a$final.box.N / p.div$N,
       box.qualities = a$final.box.quality)
  })))

  i <- order(X[, 1], X[, 2], decreasing = TRUE)
  y <- X[i, 2]
  frontier <- which(cummax(y) <= y)
  #
  # Eliminate interior points along edges of the hull
  #
  y.0 <- y[frontier]
  frontier <- frontier[c(TRUE, y.0[-1] != y.0[-length(y.0)])]
  return(i[frontier])
}

