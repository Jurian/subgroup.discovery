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
#' @param max.peel Maximal size of a peel, as a fraction. Defaults to 0.1
#' @param train.fraction Train-test split fraction used in validation, defaults to 0.66
#' @param max.boxes Maximum number of boxes, NA or leave out for no limit
#' @param quality.function Function to use for determining set quality, defaults to mean
#' @param plot Plot intermediate results, defaults to false
#' @param optimal.box During validation, choose the box with the highest quality or a simpler box, two standard errors from the optimum
#' @return An S3 object of class prim.cover
#' @author Jurian Baas
#' @importFrom stats model.frame model.response complete.cases terms
#' @importFrom graphics plot par
#' @examples
#'   data(pima)
#'   p.cov <- prim.cover(
#'       class ~ .,
#'       data = pima,
#'       peeling.quantile = 0.05,
#'       min.support = 0.1,
#'       plot = TRUE,
#'       optimal.box = "2se"
#'   )
#'
#'   \donttest{
#'   summary(p.cov)
#'   plot(p.cov)
#'   }
#' @export
prim.cover <- function (
  formula,
  data,
  X, y,
  peeling.quantile = 0.03,
  min.support = 0.05,
  max.peel = 0.1,
  train.fraction = 0.66,
  max.boxes = NA,
  quality.function = base::mean,
  plot = FALSE,
  optimal.box = c("best", "2se")) {

  using.formula <- missing(X) & missing(y)

  # Use the formula interface
  if(using.formula) {
    if(!is.data.frame(data)) stop("Data argument is not a data frame, aborting...")

    X <- stats::model.frame(formula(stats::terms(formula, data = data, simplify = TRUE)), data)
    y <- stats::model.response(X)
    X <- X[,-1]

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
  if(train.fraction <= 0) stop("Training fraction must be positive")
  if(train.fraction >= 1) stop("Training fraction must be a fraction smaller than 1")
  if(!is.na(max.boxes) & max.boxes <= 0) stop("Maximum boxes must be a positive integer")
  if(length(optimal.box) > 1) optimal.box <- optimal.box[1]
  if(!optimal.box %in% c("best", "2se")) stop("Invalid optimal box paramter")

  result <- list()
  class(result) <- "prim.cover"

  result$N <- nrow(X)
  result$peeling.quantile <- peeling.quantile
  result$min.support <- min.support
  result$train.fraction <- train.fraction
  result$quality.function <- quality.function
  result$global.quality <- quality.function(y)

  if(!is.na(max.boxes)) result$max.boxes <- max.boxes


  covers <- list()

  y.sub.quality <- result$global.quality
  min.N <- result$N * min.support
  box.nr <- 1

  repeat {

    cur.N <- nrow(X)

    # In case:
    # The user set a max number of boxes and we reached that limit
    # Current cover has become too small
    if((!is.na(max.boxes) & box.nr > max.boxes) | (cur.N < min.N)) break

    box.nr <- box.nr + 1

    train <- sample(1:cur.N, cur.N * train.fraction)

    p.peel <- prim.peel(
      X = X[train,],
      y = y[train],
      N = result$N,
      peeling.quantile = peeling.quantile,
      min.support = min.support,
      max.peel = max.peel,
      quality.function = quality.function
    )
    p.peel$global.quality <- quality.function(y[train])
    p.validate <- prim.validate(p.peel, X[-train,], y[-train], optimal.box)

    idx <- prim.rule.match(p.validate, X)

    if(sum(idx) == 0) {
      warning("Validating yielded no box, try increasing max.peel and/or decreasing min.support")
      break
    }

    p.validate$cover.N <- nrow(X)
    p.validate$cover.box.N <- sum(idx)
    p.validate$cover.global.quality <- quality.function(y)
    p.validate$cover.box.quality <- quality.function(y[idx])

    y.sub.quality <- p.validate$cover.box.quality

    # This box does not meet the minimal quality
    if(y.sub.quality < result$global.quality )  break

    if(plot) {
      graphics::par(mfrow = c(1,2))
      graphics::plot(p.peel)
      graphics::plot(p.validate)
      graphics::par(mfrow = c(1,1))
    }

    X <- X[!idx,]
    y <- y[!idx]

    covers <- c(covers, list(p.validate))

  }

  result$covers <- covers

  # Add the leftover points
  p.leftover <- list()
  class(p.leftover) <- "prim.cover.leftover"
  p.leftover$cover.N <- nrow(X)
  p.leftover$cover.global.quality <- quality.function(y)
  result$leftover <- p.leftover

  if(using.formula) result$formula <- formula

  return(result)
}

#' @title PRIM diversify strategy
#' @description Provide a (hopefully) diverse number of box definitions
#' @details Because the final box depends on the data used, we re-run the PRIM peeling algorithm multiple times, each with a different random train/test split. Each run is independent from the others, so this algorithm is run in parallel by default.
#' @param formula Formula with a response and terms
#' @param data Data frame to find rules in
#' @param X Optionally instead of using a formula: Data frame to find rules in
#' @param y Optionally instead of using a formula: Response vector, usually of type numeric
#' @param n Numer of attempts to run the PRIM algorithm
#' @param peeling.quantile Quantile to peel off for numerical variables
#' @param min.support Minimal size of a box to be valid
#' @param max.peel Maximal size of a peel, as a fraction. Defaults to 0.1
#' @param train.fraction Train-test split fraction used in validation, defaults to 0.66
#' @param quality.function Function to use for determining subset quality, defaults to mean
#' @param plot Plot intermediate results, defaults to false. Note that intermediate plotting is unavailable when running in parallel
#' @param parallel Compute each run in parallel, defaults to TRUE. This will use all but one core. Note that intermediate plotting is unavailable when running in parallel
#' @param optimal.box During validation, choose the box with the highest quality or a simpler box, two standard errors from the optimum
#' @return An S3 object of type prim.diversify
#' @author Jurian Baas
#' @importFrom stats model.frame model.response complete.cases terms
#' @importFrom graphics plot par
#' @importFrom parallel detectCores makeCluster parLapply stopCluster
#' @examples
#'   data(ames)
#'   p.div <- prim.diversify(
#'       SalePrice ~ . - PID - Order,
#'       data = ames,
#'       n = 5,
#'       plot = TRUE,
#'       parallel = FALSE,
#'       optimal.box = "best"
#'   )
#'
#'   \donttest{
#'   summary(p.div)
#'   plot(p.div)
#'   }
#' @export
prim.diversify <- function (
  formula,
  data,
  X, y, n,
  peeling.quantile = 0.03,
  min.support = 0.05,
  max.peel = 0.1,
  train.fraction = 0.66,
  quality.function = base::mean,
  plot = FALSE,
  parallel = TRUE,
  optimal.box = c("best", "2se")) {

  using.formula <- missing(X) & missing(y)

  # Use the formula interface
  if(using.formula) {
    if(!is.data.frame(data)) stop("Data argument is not a data frame, aborting...")

    X <- stats::model.frame(formula(stats::terms(formula, data = data, simplify = TRUE)), data)
    y <- stats::model.response(X)
    X <- X[,-1]

    if(is.null(y)) stop("Data has no response variable, aborting...")

  } else {

    if(!is.data.frame(X)) stop("Paremeter X has to be a data frame")
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

  if(n < 2) stop("Must supply n >= 2")
  if(peeling.quantile <= 0) stop("Peeling quantile must be positive")
  if(peeling.quantile >= 1) stop("Peeling quantile must be a fraction smaller than 1")
  if(min.support <= 0) stop("Minimum support must be positive")
  if(min.support >= 1) stop("Minimum support must be a fraction smaller than 1")
  if(train.fraction <= 0) stop("Training fraction must be positive")
  if(train.fraction >= 1) stop("Training fraction must be a fraction smaller than 1")
  if(length(optimal.box) > 1) optimal.box <- optimal.box[1]
  if(!optimal.box %in% c("best", "2se")) stop("Invalid optimal box paramter")

  if(plot & parallel) message("Note: Plotting intermediate results is unavailable while running in parallel")

  result <- list()
  class(result) <- "prim.diversify"

  attempts <- list()

  result$peeling.quantile = peeling.quantile
  result$min.support = min.support
  result$train.fraction = train.fraction
  result$quality.function <- quality.function

  result$N <- nrow(X)
  result$global.quality <- quality.function(y)

  # Define a seperate function so we can use a parallel approach if wanted
  # but also fall back on a default method
  peel.and.validate <- function(i) {
    train <- sample(1:nrow(X), nrow(X) * train.fraction)

    p.peel <- prim.peel(
      X = X[train,],
      y = y[train],
      N = result$N,
      peeling.quantile = peeling.quantile,
      min.support = min.support,
      max.peel = max.peel,
      quality.function = quality.function)
    p.peel$global.quality <- result$global.quality

    p.validate <- prim.validate(p.peel, X[-train,], y[-train], optimal.box)

    if(plot) {
      graphics::par(mfrow = c(1,2))
      graphics::plot(p.peel)
      graphics::plot(p.validate)
      graphics::par(mfrow = c(1,1))
    }

    idx <- prim.rule.match(p.validate, X)

    p.validate$final.box.N <- sum(idx)
    p.validate$final.box.quality <- quality.function(y[idx])

    return(p.validate)
  }

  if(parallel) {
    # Calculate the number of cores
    nr_cores <- parallel::detectCores() - 1

    # Initiate cluster
    cl <- parallel::makeCluster(nr_cores)

    result$attempts <- tryCatch({
      parallel::parLapply(cl, 1:n, peel.and.validate)
    }, finally = {
      # Make sure we ALWAYS stop the cluster neatly
      parallel::stopCluster(cl)
    })

  } else {
    # Fall back on the default method
    result$attempts <- lapply(1:n, peel.and.validate)
  }

  result$frontier <- quasi.convex.hull(result)
  result$scoreMatrix <- prim.diversify.compare(X, result)
  result$scores <- rowMeans(result$scoreMatrix, na.rm = TRUE)

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
#' @param N Size of entire data set
#' @param peeling.quantile Quantile to peel off for numerical variables
#' @param min.support Minimal size of a box to be valid, as a fraction
#' @param max.peel Maximal size of a peel, as fraction
#' @param quality.function Which function to use to determine the quality of a box
#' @return An S3 object of class prim.peel
#' @author Jurian Baas
#' @importFrom stats model.frame model.response
prim.peel <- function(X, y, N, peeling.quantile, min.support, max.peel, quality.function) {

  result <- list()
  result$box.qualities <- numeric()
  result$supports <- numeric()

  result$rule.names <- character()
  result$rule.operators <- character()
  result$rule.types <- character()
  result$rule.values <- list() # A list because we store multiple types of values (i.e. numerical, logical and factors)
  result$quality.function <- quality.function
  result$N <- N
  peel.N <- nrow(X)

  repeat {

    if(nrow(X) / N <= min.support) break

    # Find box candidates
    candidates <- prim.candidates.find(X, y, peeling.quantile, min.support, max.peel, quality.function)

    if(length(candidates) == 0) break

    cf  <- candidates[[which.max(sapply(candidates, function(c){c$quality}))]]

    if(cf$type == "factor") {
      cf$value <- paste0("'", cf$value, "'")
    }

    X <- X[-cf$idx,]
    y <- y[-cf$idx]

    result$box.qualities <- c(result$box.qualities, quality.function(y))
    result$supports <- c(result$supports, length(y) / peel.N)

    result$rule.values <- c(result$rule.values, cf$value)
    result$rule.names <- c(result$rule.names, cf$colname)
    result$rule.operators <- c(result$rule.operators, cf$operator)
    result$rule.types <- c(result$rule.types, cf$type)
  }

  result$rule.names <- factor(result$rule.names, levels = colnames(X))
  result$rule.operators <- factor(result$rule.operators, levels = c(">=", "<=", "==", "!="))
  result$rule.types <- factor(result$rule.types, levels = c("numeric", "logical", "factor"))

  class(result) <- "prim.peel"

  #result$superrule <- prim.rule.condense(result)
  result$call <- match.call()

  return(result)
}

#' @title Bump hunting using the Patient Rule Induction Method
#' @description Validate the results taken from the PRIM peeling process
#' @details This function takes the result of the prim peeling process and applies it to new data. Usually the optimal box in the peeling process is not the best on unobserved data.
#' @param peel.result An S3 object of class prim.peel
#' @param X A data frame with at least those columns that were used in creating the prim.peel S3 object
#' @param y Response vector, usually of type numeric
#' @param optimal.box Choose the box with the highest quality or a simpler box, two standard errors from the optimum
#' @return An S3 object of type prim.validate
#' @author Jurian Baas
prim.validate <- function(peel.result, X, y, optimal.box) {

  if(class(peel.result) != "prim.peel") {
    stop("Supplied argument is not of class prim.peel, aborting...")
  }

  if(nrow(X) != length(y)) stop(paste("Parameters X and y are not of same size:", nrow(X), length(y)))

  quality.function <- peel.result$quality.function

  result <- list()
  class(result) <- "prim.validate"
  result$box.qualities <- numeric()
  result$supports <- numeric()
  result$redundant <- integer()
  result$rule.names <- character()
  result$rule.operators <- character()
  result$rule.types <- character()
  result$rule.values <- list() # A list because we store multiple types of values (i.e. numerical, logical and factors)
  result$quality.function <- peel.result$quality.function
  result$N <- peel.result$N
  result$global.quality <- peel.result$global.quality
  result$optimal.box <- optimal.box
  validate.N <- nrow(X)

  i <- 1

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
      result$supports <- c(result$supports, length(y) / validate.N)

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

  result$metrics <- prim.validate.metrics(result)
  result$best.box.idx <- prim.box.optimal(result)
  result$superrule <- prim.rule.condense(result)
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
#' @param max.peel Maximal size of a peel
#' @param quality.function Function to use to determine box quality
#' @return A list of potential boxes
#' @author Jurian Baas
#' @importFrom stats quantile
prim.candidates.find <- function(X, y, peeling.quantile, min.support, max.peel, quality.function) {

  i <- 1
  candidates <- list()
  cnames <- colnames(X)

  N <- nrow(X)

  repeat {

    if(i > ncol(X)) break

    r <- list()

    col <- X[,i]

    # Do something different depending on data type
    if(is.numeric(col)) {

      quantiles <- stats::quantile(col, c(peeling.quantile, 1 - peeling.quantile), names = FALSE)

      quantile.left <- quantiles[1]
      quantile.right <- quantiles[2]

      idx.left <- col < quantile.left
      idx.right <- col > quantile.right

      left.support <- sum(idx.left)
      right.support <- sum(idx.right)

      left.support.complement <- N - left.support
      right.support.complement <- N - right.support

      if(left.support > 0 & left.support / N <= max.peel & left.support.complement / N >= min.support) {
        r <- c(r, list( list (
          value = quantile.left,
          operator = ">=",
          type = "numeric",
          quality = quality.function(y[!idx.left]),
          idx = which(idx.left),
          colname = cnames[i]
        )))
      }

      if(right.support > 0 & right.support / N <= max.peel  & right.support.complement / N >= min.support) {
        r <- c(r, list( list (
          value = quantile.right,
          operator = "<=",
          type = "numeric",
          quality = quality.function(y[!idx.right]),
          idx = which(idx.right),
          colname = cnames[i]
        )))
      }

    } else if (is.logical(col)  ) {

      # Don't check if column consists of only T or F
      if(any(col) & any(!col)) {

        idx.true <- col
        idx.false <- !col

        support.true <- sum(col)
        support.false <- sum(!col)

        if(support.true / N <= max.peel & support.false / N >= min.support) {
          r <- c(r, list( list (
            value = FALSE,
            operator = "==",
            type = "logical",
            quality = quality.function(y[!idx.true]),
            idx = which(idx.true),
            colname = cnames[i]
          )))
        }
        if(support.false / N <= max.peel & support.true / N >= min.support) {
          r <- c(r, list( list (
            value = TRUE,
            operator = "==",
            type = "logical",
            quality = quality.function(y[!idx.false]),
            idx = which(idx.false),
            colname = cnames[i]
          )))
        }
      }

    } else if (is.factor(col)) {

      lvls <- levels(col)
      j <- 1

      repeat {

        if(j > length(lvls)) break

        lvl <- lvls[j]

        idx <- col == lvl

        support.lvl <- sum(idx)

        if(support.lvl > 0) {

          support.lvl.complement <- N - support.lvl

          if(support.lvl / N <= max.peel & support.lvl.complement / N >= min.support) {
            r <- c(r, list( list (
              value = lvl,
              operator = "!=",
              type = "factor",
              quality = quality.function(y[!idx]),
              idx = which(idx),
              colname = cnames[i]
            )))
          }

        }


        j <- j + 1
      }

    }


    if(length(r) > 0) {
      candidates <- c(candidates, r)
    }

    i <- i + 1

  }

  return(candidates)
}


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

  frontier <- sort(p.div$frontier)
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

  for(i in 1:(best.box.idx-1)) {
    if(x$rule.names[i] == x$rule.names[i + 1]) {
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
      if(x$rule.names[i] == x$rule.names[i + 1]) {
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
  cat("  |  Peeling quantile:", object$peeling.quantile, "\n")
  cat("  |  Min support:", object$min.support, "\n")
  cat("  |  Train/test split:", object$train.fraction, "\n")
  cat("\n")

  for(i in 1:length(object$covers)) {
    x <- object$covers[[i]]
    cat("\n")
    cat("  ======================================", "\n")
    cat("  ============== COVER", i,"===============", "\n")
    cat("  |  Cover set size: ", x$cover.N, "\n")
    cat("  |  Cover set quality: ", round(x$cover.global.quality, digits), "\n")
    cat("  |\n")
    cat("  |  Box relative quality: ", round(x$cover.box.quality, digits), "(", round(x$cover.box.quality / x$cover.global.quality, digits), ") \n")
    cat("  |  Box relative support: ", round(x$cover.box.N / x$cover.N, digits) , " (", x$cover.box.N, ") \n")
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
  cat("  |  Peeling quantile:", object$peeling.quantile, "\n")
  cat("  |  Min support:", object$min.support, "\n")
  cat("  |  Train/test split:", object$train.fraction, "\n")
  cat("  |\n")
  cat("  |  Set size: ", object$N, "\n")
  cat("  |  Set quality: ", round(object$global.quality, digits), "\n")
  cat("\n")
  cat("  Scores:", "\n  | ")
  cat(paste0(
    gsub(
      "NA",
      "    ",
      apply(formatC(round(object$scoreMatrix, 2), format = "f", digits = 2), 1, paste, collapse = "\t"))
    , collapse = "\n  |  ")
    )
  cat("\n\n\n")
  cat("  Dominating attempts:","\n")
  frontier <- sort(object$frontier)
  for(i in seq_along(frontier)) {
    x <- object$attempts[[frontier[i]]]
    cat("\n")
    cat("  ======================================", "\n")
    cat("  ============= ATTEMPT", frontier[i],"==============", "\n")
    cat("  |  Score:", round(object$scores[i], digits), "\n")
    cat("  |  Box quality: ", round(x$final.box.quality, digits), "(", round(x$final.box.quality / object$global.quality, digits), ") \n")
    cat("  |  Box support: ", round(x$final.box.N / object$N, digits) , " (", x$final.box.N, ") \n")
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

