
data.column.create <- function(data) {

  data.column <- list()
  class(data.column) <- "data.column"

  data.column$ordering    <- base::order(data)
  data.column$data        <- data[data.column$ordering]
  data.column$orderindex  <- base::order(data.column$ordering)

  return(data.column)
}

data.column.update <- function(data.column, index) {

  local.index <- data.column$orderindex[index]

  # First idea, replace with NA
  data.column$data[local.index] <- NA
  data.column$ordering[local.index] <- NA
  data.column$orderindex[index] <- NA

  # Second idea, remove values and reconstruct orderindex
  #data.column$data <- data.column$data[-local.index]
  #data.column$ordering <- data.column$ordering[-local.index]
  # TODO: replace with something more efficient
  #data.column$orderindex <- base::match(1:length(data.column$orderindex), data.column$ordering)

  return(data.column)
}

as.matrix.data.column <- function(data.column) {
  base::cbind (
    data = data.column$data,
    ordering = data.column$ordering,
    orderindex = data.column$orderindex
  )
}

print.data.column <- function(data.column) {
  base::print(as.matrix.data.column(data.column))
}

