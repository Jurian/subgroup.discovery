#-------------------------------------------------------------------------------------------#
################################## PRIM INTERFACE FUNCTIONS #################################
#-------------------------------------------------------------------------------------------#


#' @title PRIM covering strategy
#' @description In bump hunting it is customary to follow a so-called covering strategy. This means that the same box construction (rule induction) algorithm is applied sequentially to subsets of the data.
#' @param formula Formula with a response and terms
#' @param data Data frame to find rules in
#' @param X Optionally instead of using a formula: Data frame to find rules in
#' @param y Optionally instead of using a formula: Response vector, usually of type numeric
#' @param peeling.quantile Quantile to peel off for numerical variables
#' @param min.support Minimal size of a box to be valid
#' @param train.fraction Train-test split fraction used in validation
#' @param max.boxes Optional maximum number of boxes
#' @param quality.function Optional setting for function to use for determining subset quality, defaults to mean
#' @param plot Optional setting to plot intermediate results, defaults to false
#' @return An S3 object of class prim.cover
#' @author Jurian Baas
#' @importFrom stats model.frame model.response
#' @importFrom graphics plot par
#' @export
prim.cover <- function(formula, data, X = NULL, y = NULL, peeling.quantile, min.support, train.fraction = 0.66, max.boxes = NA, quality.function = base::mean, plot = FALSE) {

  using.formula <- is.null(X) & is.null(y)

  # Use the formula interface
  if(using.formula) {
    if(!is.data.frame(data)) stop("Data argument is not a data frame, aborting...")

    X <- stats::model.frame(formula = formula, data = data)
    y <- stats::model.response(X)
    X <- X[,-1]

    if(is.null(y)) stop("Data has no response variable, aborting...")

  } else {

    if(!is.data.frame(X)) stop("Parameter X has to be a data frame")
    if(!is.vector(y)) stop("Parameter y has to be a vector")
    if(nrow(X) != length(y)) stop("Parameters X and y are not of same size")
  }

  if(peeling.quantile <= 0) stop("Peeling quantile must be positive")
  if(peeling.quantile >= 1) stop("Peeling quantile must be a fraction smaller than 1")
  if(min.support <= 0) stop("Minimum support must be positive")
  if(min.support >= 1) stop("Minimum support must be a fraction smaller than 1")
  if(train.fraction <= 0) stop("Training fraction must be positive")
  if(train.fraction >= 1) stop("Training fraction must be a fraction smaller than 1")
  if(!is.na(max.boxes) & max.boxes <= 0) stop("Maximum boxes must be a positive integer")

  N <- nrow(X) * min.support
  box.nr <- 1
  covers <- list()

  result <- list()
  result$peeling.quantile <- peeling.quantile
  result$min.support <- min.support
  result$train.fraction <- train.fraction
  result$quality.function <- quality.function
  if(!is.na(max.boxes)) result$max.boxes <- max.boxes
  class(result) <- "prim.cover"

  y.quality <- quality.function(y)

  repeat {

    y.sub.quality <- quality.function(y)

    # In case:
    # The user set a max nr of boxes and we reached that limit
    # Current cover has become too small
    # Current cover quality fell below overall quality
    if((!is.na(max.boxes) & box.nr > max.boxes) | (nrow(X) < N) | y.sub.quality < y.quality ) {
      p.leftover <- list()
      class(p.leftover) <- "prim.cover.leftover"

      p.leftover$cov.support <- nrow(X)
      p.leftover$cov.overall.quality <- quality.function(y)

      result$leftover <- p.leftover
      break
    }

    train <- sample(1:nrow(X), nrow(X) * train.fraction)

    p.train <- prim.peel(X = X[train,], y = y[train], peeling.quantile = peeling.quantile, min.support = min.support, quality.function = quality.function)
    p.validate <- prim.validate(p.train, X[-train,], y[-train])

    if(plot) {
      graphics::par(mfrow = c(1,2))
      graphics::plot(p.train)
      graphics::plot(p.validate)
    }

    idx <- prim.match.rule(p.validate, X)

    p.validate$cov.support <- nrow(X)
    p.validate$cov.box.support <- sum(idx)
    p.validate$cov.overall.quality <- quality.function(y)
    p.validate$cov.box.quality <- quality.function(y[idx])

    X <- X[!idx,]
    y <- y[!idx]
    box.nr <- box.nr + 1

    covers <- c(covers, list(p.validate))
  }

  result$covers <- covers

  if(using.formula) result$formula <- formula

  return(result)
}

#' @title PRIM diversify strategy
#' @description Provide a (hopefully) diverse number of box definitions
#' @details Because the final box depends on the data used, we re-run the PRIM peeling algorithm multiple times, each with a different random train/test split.
#' @param formula Formula with a response and terms
#' @param data Data frame to find rules in
#' @param X Optionally instead of using a formula: Data frame to find rules in
#' @param y Optionally instead of using a formula: Response vector, usually of type numeric
#' @param n Numer of attempts to run the PRIM algorithm
#' @param peeling.quantile Quantile to peel off for numerical variables
#' @param min.support Minimal size of a box to be valid
#' @param train.fraction Optional train-test split fraction used in validation, defaults to 0.66
#' @param quality.function Optional setting for function to use for determining subset quality, defaults to mean
#' @param plot Optional setting to plot intermediate results, defaults to false
#' @return An S3 object of type prim.diversify
#' @author Jurian Baas
#' @importFrom stats model.frame model.response
#' @importFrom graphics plot par
#' @export
prim.diversify <- function(formula, data, X = NULL, y = NULL, n, peeling.quantile, min.support, train.fraction = 0.66, quality.function = base::mean, plot = FALSE) {

  using.formula <- is.null(X) & is.null(y)

  # Use the formula interface
  if(using.formula) {
    if(!is.data.frame(data)) stop("Data argument is not a data frame, aborting...")

    X <- stats::model.frame(formula = formula, data = data)
    y <- stats::model.response(X)
    X <- X[,-1]

    if(is.null(y)) stop("Data has no response variable, aborting...")

  } else {

    if(!is.data.frame(X)) stop("Paremeter X has to be a data frame")
    if(!is.vector(y)) stop("Parameter y has to be a vector")
    if(nrow(X) != length(y)) stop("Parameters X and y are not of same size")
  }

  if(n < 2) stop("Must supply n >= 2")
  if(peeling.quantile <= 0) stop("Peeling quantile must be positive")
  if(peeling.quantile >= 1) stop("Peeling quantile must be a fraction smaller than 1")
  if(min.support <= 0) stop("Minimum support must be positive")
  if(min.support >= 1) stop("Minimum support must be a fraction smaller than 1")
  if(train.fraction <= 0) stop("Training fraction must be positive")
  if(train.fraction >= 1) stop("Training fraction must be a fraction smaller than 1")

  result <- list()
  attempts <- list()

  result$peeling.quantile = peeling.quantile
  result$min.support = min.support
  result$train.fraction = train.fraction
  result$quality.function <- quality.function

  result$support <- nrow(X)
  result$overall.quality <- quality.function(y)

  for(i in 1:n) {

    train <- sample(1:nrow(X), nrow(X) * train.fraction)

    p.train <- prim.peel(X = X[train,], y = y[train], peeling.quantile = peeling.quantile, min.support = min.support, quality.function = quality.function)
    p.validate <- prim.validate(p.train, X[-train,], y[-train])

    if(plot) {
      graphics::par(mfrow = c(1,2))
      graphics::plot(p.train)
      graphics::plot(p.validate)
    }

    idx <- prim.match.rule(p.validate, X)

    p.validate$att.box.support <- sum(idx)
    p.validate$att.box.quality <- quality.function(y[idx])

    attempts <- c(attempts, list(p.validate))
  }

  result$attempts <- attempts
  class(result) <- "prim.diversify"

  if(using.formula) result$formula <- formula

  return(result)
}







#-------------------------------------------------------------------------------------------#
################################### PRIM PRIVATE FUNCTIONS ##################################
#-------------------------------------------------------------------------------------------#

#' @title Bump hunting using the Patient Rule Induction Method
#' @description Peeling function for bump hunting using the Patient Rule Induction Method (PRIM).
#' @param X Data frame to find rules in
#' @param y Response vector, usually of type numeric
#' @param peeling.quantile Quantile to peel off for numerical variables
#' @param min.support Minimal size of a box to be valid, as a fraction
#' @param quality.function Which function to use to determine the quality of a box, defaults to mean
#' @return An S3 object of class prim.peel
#' @author Jurian Baas
#' @importFrom stats model.frame model.response
prim.peel <- function(X, y, peeling.quantile, min.support, quality.function = base::mean) {


  if(!is.data.frame(X)) stop("Paremeter X has to be a data frame")
  if(!is.vector(y)) stop("Parameter y has to be a vector")
  if(nrow(X) != length(y)) stop("Parameters X and y are not of same size")
  if(peeling.quantile <= 0) stop("Peeling quantile must be positive")
  if(peeling.quantile >= 1) stop("Peeling quantile must be a fraction smaller than 1")
  if(min.support <= 0) stop("Minimum support must be positive")
  if(min.support >= 1) stop("Minimum support must be a fraction smaller than 1")

  result <- list()
  result$box.qualities <- quality.function(y)
  result$supports <- 1

  result$rule.names <- character()
  result$rule.operators <- character()
  result$rule.types <- character()
  result$rule.values <- list() # A list because we store multiple types of values (i.e. numerical, logical and factors)
  result$quality.function <- quality.function
  result$peel.support <- nrow(X)
  result$peel.overall.quality <- quality.function(y)

  repeat {

    if(nrow(X) / result$peel.support <= min.support) break

    # Find box candidates
    candidates <- prim.candidates.find(X, y, peeling.quantile, min.support, result$peel.support, quality.function)

    if(length(candidates) == 0) break

    cf <- prim.candidates.best(candidates)

    if(cf$type == "factor") {
      cf$value <- paste0("'", cf$value, "'")
    }

    X <- X[-cf$idx,]
    y <- y[-cf$idx]

    result$box.qualities <- c(result$box.qualities, quality.function(y))
    result$supports <- c(result$supports, length(y) / result$peel.support)

    result$rule.values <- c(result$rule.values, cf$value)
    result$rule.names <- c(result$rule.names, cf$colname)
    result$rule.operators <- c(result$rule.operators, cf$operator)
    result$rule.types <- c(result$rule.types, cf$type)
  }

  result$rule.names <- factor(result$rule.names, levels = colnames(X))
  result$rule.operators <- factor(result$rule.operators, levels = c(">=", "<=", "==", "!="))
  result$rule.types <- factor(result$rule.types, levels = c("numeric", "logical", "factor"))

  class(result) <- "prim.peel"

  result$superrule <- prim.condense.rules(result)
  result$call <- match.call()

  return(result)
}

#' @title Bump hunting using the Patient Rule Induction Method
#' @description Validate the results taken from the PRIM peeling process
#' @details This function takes the result of the prim peeling process and applies it to new data. Usually the optimal box in the peeling process is not the best on unobserved data.
#' @param peel.result An S3 object of class prim.peel
#' @param X A data frame with at least those columns that were used in creating the prim.peel S3 object
#' @param y Response vector, usually of type numeric
#' @return An S3 object of type prim.validate
#' @author Jurian Baas
prim.validate <- function(peel.result, X, y) {

  if(class(peel.result) != "prim.peel") {
    stop("Supplied argument is not of class prim.peel, aborting...")
  }

  if(!is.data.frame(X)) stop("Paremeter X has to be a data frame")
  if(!is.vector(y)) stop("Parameter y has to be a vector")
  if(nrow(X) != length(y)) stop("Parameters X and y are not of same size")

  quality.function <- peel.result$quality.function

  i <- 1

  result <- list()
  class(result) <- "prim.validate"
  result$box.qualities <- quality.function(y)
  result$supports <- 1
  result$redundant <- integer()
  result$rule.names <- character()
  result$rule.operators <- character()
  result$rule.types <- character()
  result$rule.values <- list() # A list because we store multiple types of values (i.e. numerical, logical and factors)
  result$quality.function <- quality.function
  result$validate.support <- nrow(X)
  result$validate.overall.quality <- quality.function(y)
  result$peel.result <- peel.result

  repeat {

    # Stop if there are no more boxes
    if(i > length(peel.result$rule.names)) break

    rule.info <- list(
      name = as.character(peel.result$rule.names[i]),
      operator = as.character(peel.result$rule.operators[i]),
      value = peel.result$rule.values[i],
      type = as.character(peel.result$rule.types[i])
    )

    rule <- paste( rule.info$name, rule.info$operator, rule.info$value)
    idx <- !with(X, eval(parse(text = rule)))

    # Check if this rule has any observations in the test data
    if(sum(idx) > 0) {

      X <- X[!idx,]
      y <- y[!idx]

      result$box.qualities <- c(result$box.qualities, quality.function(y))
      result$supports <- c(result$supports, length(y) / result$validate.support)

      result$rule.values <- c(result$rule.values, rule.info$value)
      result$rule.names <- c(result$rule.names, rule.info$name)
      result$rule.operators <- c(result$rule.operators, rule.info$operator)
      result$rule.types <- c(result$rule.types, rule.info$type)

    } else {

      result$redundant <- c(result$redundant, i)

    }

    # This box could have all the remaining observations, no need to iterate further
    if(nrow(X) == 0) break

    i <- i + 1
  }

  result$rule.names <- factor(result$rule.names, levels = colnames(X))
  result$rule.operators <- factor(result$rule.operators, levels = c(">=", "<=", "==", "!="))
  result$rule.types <- factor(result$rule.types, levels = c("numeric", "logical", "factor"))

  result$superrule <- prim.condense.rules(result)
  result$call <- match.call()

  return(result)
}

#' @title PRIM find split candidates
#' @description Find all box candidates for a given (sub)set
#' @details This function goes through all columns of the dataset and tries to findbox candidates based on the quantile peeling.quantile and minimum support min.support. Note that the indexes returned are those that have to be removed in order to create the box!
#' @param X Data frame with observations (may be a subset of original data)
#' @param y Dependent variable, usually a numeric vector
#' @param peeling.quantile Quantile to peel off
#' @param min.support Minimal size of a box
#' @param support Support of subset
#' @param quality.function Function to use to determine box quality
#' @return A list of potential boxes
#' @author Jurian Baas
#' @importFrom stats quantile
prim.candidates.find <- function(X, y, peeling.quantile, min.support, support, quality.function) {

  i <- 1
  candidates <- list()
  cnames <- colnames(X)

  repeat {

    if(i > ncol(X)) break

    r <- list()

    col <- X[,i]

    # Do something different depending on data type
    if(is.numeric(col)) {

      quantiles <- stats::quantile(col, c(peeling.quantile, 1 - peeling.quantile) , names = FALSE, na.rm = T)

      quantile.min <- quantiles[1]
      quantile.plus <- quantiles[2]

      idx.min <- col < quantile.min
      idx.plus <- col > quantile.plus

      if(sum(idx.min) > 0 & (length(col) - sum(idx.min)) / support >= min.support) {
        r <- c(r,
               list(
                 min = list(
                   value = quantile.min,
                   operator = ">=",
                   type = "numeric",
                   quality = quality.function(y[!idx.min]),
                   idx = which(idx.min),
                   size = sum(idx.min) / length(col),
                   colname = cnames[i]
                 )
               )
        )
      }

      if(sum(idx.plus) > 0 & (length(col) - sum(idx.plus)) / support >= min.support) {
        r <- c(r,
               list(
                 plus = list(
                   value = quantile.plus,
                   operator = "<=",
                   type = "numeric",
                   quality = quality.function(y[!idx.plus]),
                   idx = which(idx.plus),
                   size = sum(idx.plus) / length(col),
                   colname = cnames[i]
                 )
               )
        )
      }

    } else if (is.logical(col)  ) {

      # Don't check if column consists of only T or F
      if(any(col) & any(!col)) {

        idx.min <- col
        idx.plus <- !col

        # Make sure the resulting subsets have enough support
        if((sum(!col) / support) >= min.support) {
          r <- c(r,
                 list(
                   min = list(
                     value = FALSE,
                     operator = "==",
                     type = "logical",
                     quality = quality.function(y[!idx.min]),
                     idx = which(idx.min),
                     size = sum(idx.min) / length(col),
                     colname = cnames[i]
                   )
                 )
          )
        }
        if((sum(col) / support) >= min.support) {
          r <- c(r,
                 list(
                   plus = list(
                     value = TRUE,
                     operator = "==",
                     type = "logical",
                     quality = quality.function(y[!idx.plus]),
                     idx = which(idx.plus),
                     size = sum(idx.plus) / length(col),
                     colname = cnames[i]
                   )
                 )
          )
        }
      }

    } else if (is.factor(col)) {

      # Reset levels just to be sure, so we don't use levels that no
      # longer exist due to them being removed in some previous iteration
      col <- factor(col)
      lvls <- levels(col)
      j <- 1

      repeat {

        if(j > length(lvls)) break

        lvl <- lvls[j]

        idx <- col == lvl

        if((length(col) - sum(idx)) / support >= min.support) {
          r[[lvl]] <- list(
            value = lvl,
            operator = "!=",
            type = "factor",
            quality = quality.function(y[!idx]),
            idx = which(idx),
            size = sum(idx) / length(col),
            colname = cnames[i]
          )
        }

        j <- j + 1
      }

    }


    if(length(r) > 0) {
      candidates[[cnames[i]]] <- r
    }

    i <- i + 1

  }

  return(candidates)
}

#' @title Find optimal box candidate
#' @description This function goes through the box candidate qualities and finds the optimal candidate
#' @param candidates List of candidate generated by prim.candidates.find()
#' @return A list with the optimal candidate information
#' @author Jurian Baas
prim.candidates.best <- function(candidates) {

  qualities <- sapply(candidates, function(col) {
    sapply(col, function(x) x$quality)
  })
  supports <- sapply(candidates, function(col) {
    sapply(col, function(x) x$size)
  })

  candidate.best <- NULL
  quality.max <- NA
  support.min <- NA

  # Find the best candidate by linear search
  for(i in 1:length(qualities)) {

    col <- qualities[[i]]
    sup <- supports[[i]]

    for(j in 1:length(col)) {

      if(is.na(quality.max) | col[j] > quality.max) {

        quality.max <- col[j]
        support.min <- sup[j]
        candidate.best <- candidates[[i]][[j]]

      } else if (!is.na(quality.max) & col[j] == quality.max & sup[j] < support.min) { # Break ties by their support

        quality.max <- col[j]
        support.min <- sup[j]
        candidate.best <- candidates[[i]][[j]]

      }

    }
  }
  return(candidate.best)
}

#' @title Create box matching index
#' @description Generate a logical vector in which elements are true iff applying the rule to a record evaluates to true.
#' @param prim.object An S3 object of class prim.peel or prim.validate result
#' @param X A data frame with at least those columns that were used in creating the prim S3 object
#' @return A logical index of matching records
#' @author Jurian Baas
prim.match.rule <- function(prim.object, X) {

  if(class(prim.object) != "prim.peel" & class(prim.object) != "prim.validate") {
    stop("Supplied argument is not of class prim.peel or prim.validate, aborting...")
  }
  return(with(X, eval(parse(text = paste0(prim.object$superrule, collapse = " & ")))))
}

#' @title Condense multiple (redundant) rules
#' @details This function condenses the many (redundant) rules of an S3 object of class prim.peel or prim.validate to a single rule.
#' @description This function condenses the many (redundant) rules of an S3 object of class prim.peel or prim.validate to a single rule.
#' @param prim.object An S3 object of class prim.peel or prim.validate
#' @return The condensed rule as a single string
#' @author Jurian Baas
prim.condense.rules <- function(prim.object) {

  if(class(prim.object) != "prim.peel" & class(prim.object) != "prim.validate"){
    stop("Supplied argument is not of class prim.peel or prim.validate, aborting...")
  }

  # Search for the rules leading up to the best box
  rule.idx <- 1:(which.max(prim.object$box.qualities)-1)

  # Combine the lists into a single data frame
  t <- data.frame (
    name = factor(prim.object$rule.names[rule.idx]),
    operator = factor(prim.object$rule.operators[rule.idx]),
    value = as.character(prim.object$rule.values[rule.idx]),
    score = prim.object$box.qualities[rule.idx],
    type = factor(prim.object$rule.types[rule.idx])
  )

  # Group by name first
  t <- by(t, t$name, function(x) {

    # Group by operator in the name subgroup
    by(x, factor( x$operator), function(y) {
      # Find the rule with the largest associated score
      y <- y[which.max( y$score),]
      # Return this rule as a proper string
      return(paste(y$name, y$operator, y$value))
    })
  })

  # Collapse to character vector
  return(unname(unlist(as.vector(t))))
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

    idx <- prim.match.rule(cover, newdata)

    y.hat[idx & is.na(y.hat)] <- cover$cov.box.quality

    i <- i + 1
  }

  y.hat[is.na(y.hat)] <- object$leftover$cov.overall.quality

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
  graphics::par(bty = "l")
  graphics::plot (
    x$supports,
    x$box.qualities,
    xlim = c(1, 0),
    type = "n",
    xlab = "Support", ylab = "Box quality",
    main = "PRIM peel result",
    ...)
  graphics::lines (
    x$supports,
    x$box.qualities
  )
  graphics::points (
    x$supports,
    x$box.qualities,
    col= "royalblue4", pch = 19, cex = 1, lty = "solid", lwd = 2)
  graphics::text (
    x$supports,
    x$box.qualities,
    labels = c("", paste(x$rule.names, x$rule.operators, x$rule.values) ),
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
  graphics::par(bty = "l")
  graphics::plot (
    x$supports,
    x$box.qualities,
    xlim = c(1, 0),
    type = "n",
    xlab = "Support", ylab = "Box quality",
    main = "PRIM validate result",
    ...)
  graphics::lines (
    x$supports,
    x$box.qualities
  )
  graphics::points (
    x$supports,
    x$box.qualities,
    col= "royalblue4", pch = 19, cex = 1, lty = "solid", lwd = 2)
  graphics::text (
    x$supports,
    x$box.qualities,
    labels = c("", paste(x$rule.names, x$rule.operators, x$rule.values) ),
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
    c( supports = c$cov.box.support / c$cov.support,
       box.qualities = c$cov.box.quality)
  })))

  graphics::par(bty = "l")
  graphics::plot (
    dat$supports,
    dat$box.qualities,
    type = "n",
    xlim = c(0,1),
    xlab = "Relative support", ylab = "Relative box quality",
    main = "PRIM cover result",
    ...)
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
    labels = paste("Cover", 1 : length(x$covers)),
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
    c( supports = a$att.box.support / x$support,
       box.qualities = a$att.box.quality)
  })))
  graphics::par(bty = "l")
  graphics::plot (
    dat$supports,
    dat$box.qualities,
    type = "n",
    xlab = "Support", ylab = "Box quality",
    main = "PRIM diversify result",
    ...)
  graphics::points (
    dat$supports,
    dat$box.qualities,
    col= "royalblue4", pch = 19, cex = 1, lty = "solid", lwd = 2)
  graphics::text (
    dat$supports,
    dat$box.qualities,
    labels = paste("Attempt", 1 : length(x$attempts)),
    cex = 0.7, pos = 1, col = "orangered4", font = 2)
}







#-------------------------------------------------------------------------------------------#
################################### S3 SUMMARY FUNCTIONS ####################################
#-------------------------------------------------------------------------------------------#

#' @title Summarize a PRIM peeling result object
#' @description Summarize a PRIM peeling result object
#' @param object An S3 object of class prim.peel
#' @param ... Optional arguments to pass on
#' @param round Optional setting to disable rounding
#' @param digits Optional setting to control nr of digits to round
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

  cat("  Set size: ", object$peel.support, "\n")
  cat("  Set quality: ", round(object$peel.overall.quality, digits), "\n")
  cat("\n")
  cat("  ============== BEST BOX ==============", "\n")
  cat("  Box quality: ", round(object$box.qualities[best.box.idx], digits), "(", round(object$box.qualities[best.box.idx] / object$peel.overall.quality, digits), ") \n")
  cat("  Box support: ", round(object$supports[best.box.idx], digits), " (", object$supports[best.box.idx] * object$peel.support, ") \n")
  cat("\n")
  cat("  ================ RULES ===============", "\n")
  cat(" ", paste0(object$superrule, collapse = "\n  "))

}

#' @title Summarize a PRIM test result object
#' @description Summarize a PRIM test result object
#' @param object An S3 object of class prim.validate
#' @param ... Optional arguments to pass on
#' @param round Optional setting to disable rounding
#' @param digits Optional setting to control nr of digits to round
#' @author Jurian Baas
#' @return Nothing, this function is called for its side-effects
#' @export
summary.prim.validate <- function(object, ..., round = TRUE, digits = 2) {

  if(!round)  digits = 7

  cat("  ======================================", "\n")
  cat("  ========== PRIM TEST RESULT ==========", "\n")
  cat("  ======================================", "\n")
  cat("\n")
  best.box.idx <- which.max(object$box.qualities)

  cat("  Original set size: ", object$validate.support, "\n")
  cat("  Original set quality: ", round(object$validate.overall.quality, digits), "\n")
  cat("\n")
  cat("  ============== BEST BOX ==============", "\n")
  cat("  Box quality: ", round(object$box.qualities[best.box.idx], digits), "(", round(object$box.qualities[best.box.idx] / object$validate.overall.quality, digits), ") \n")
  cat("  Box support: ", round(object$supports[best.box.idx], digits), " (", object$supports[best.box.idx] * object$validate.support, ") \n")
  cat("\n")
  cat("  ================ RULES ===============", "\n")
  cat(" ", paste0(object$superrule, collapse = "\n  "))

}

#' @title Summarize a PRIM cover result object
#' @description Summarize a PRIM cover result object
#' @param object An S3 object of class prim.cover
#' @param ... Optional arguments to pass on
#' @param round Optional setting to disable rounding
#' @param digits Optional setting to control nr of digits to round
#' @author Jurian Baas
#' @return Nothing, this function is called for its side-effects
#' @export
summary.prim.cover <- function(object, ..., round = TRUE, digits = 2) {

  if(!round)  digits = 7

  cat("  ======================================", "\n")
  cat("  ========== PRIM COVER RESULT =========", "\n")
  cat("  ======================================", "\n")
  cat("  |\n")
  cat("  |  Peeling quantile:", object$peeling.quantile, "\n")
  cat("  |  Min support:", object$min.support, "\n")
  cat("  |  Train/test split:", object$train.fraction, "\n")
  cat("\n")

  for(i in 1:length(object$covers)) {
    x <- object$covers[[i]]
    cat("\n")
    cat("  ======================================", "\n")
    cat("  ============== COVER", i,"===============", "\n")
    cat("  |  Cover set size: ", x$cov.support, "\n")
    cat("  |  Cover set quality: ", round(x$cov.overall.quality, digits), "\n")
    cat("  |\n")
    cat("  |  Box relative quality: ", round(x$cov.box.quality, digits), "(", round(x$cov.box.quality / x$cov.overall.quality, digits), ") \n")
    cat("  |  Box relative support: ", round(x$cov.box.support / x$cov.support, digits) , " (", x$cov.box.support, ") \n")
    cat("\n")
    cat("  ================ RULES ===============", "\n")
    cat("  | ", paste0(x$superrule, collapse = "\n  |  "))
    cat("\n","\n","\n")

  }

  x <- object$leftover

  cat("\n")
  cat("  ======================================", "\n")
  cat("  ============== LEFTOVER ==============", "\n")
  cat("  |  Cover set size: ", x$cov.support, "\n")
  cat("  |  Cover set quality: ", round(x$cov.overall.quality, digits), "\n")
  cat("\n","\n","\n")

}

#' @title Summarize a PRIM diversify object
#' @description Summarize a PRIM diversify result object
#' @param object An S3 object of class prim.diversify
#' @param ... Optional arguments to pass on
#' @param round Optional setting to disable rounding
#' @param digits Optional setting to control nr of digits to round
#' @author Jurian Baas
#' @return Nothing, this function is called for its side-effects
#' @export
summary.prim.diversify <- function(object, ..., round = TRUE, digits = 2) {

  if(!round) digits = 7

  cat("  ======================================", "\n")
  cat("  ======== PRIM DIVERSIFY RESULT =======", "\n")
  cat("  ======================================", "\n")
  cat("  |\n")
  cat("  |  Peeling quantile:", object$peeling.quantile, "\n")
  cat("  |  Min support:", object$min.support, "\n")
  cat("  |  Train/test split:", object$train.fraction, "\n")
  cat("  |\n")
  cat("  |  Set support: ", object$support, "\n")
  cat("  |  Set quality: ", round(object$overall.quality, digits), "\n")
  cat("\n")

  for(i in 1:length(object$attempts)) {
    x <- object$attempts[[i]]
    cat("\n")
    cat("  ======================================", "\n")
    cat("  ============= ATTEMPT", i,"==============", "\n")
    cat("  |  Box quality: ", round(x$att.box.quality, digits), "(", round(x$att.box.quality / object$overall.quality, digits), ") \n")
    cat("  |  Box support: ", round(x$att.box.support / object$support, digits) , " (", x$att.box.support, ") \n")
    cat("\n")
    cat("  ================ RULES ===============", "\n")
    cat("  | ", paste0(x$superrule, collapse = "\n  |  "))
    cat("\n","\n","\n")

  }
}
