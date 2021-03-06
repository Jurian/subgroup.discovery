% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prim.R
\name{prim.data.prepare}
\alias{prim.data.prepare}
\title{Prepare data for PRIM}
\usage{
prim.data.prepare(X)
}
\arguments{
\item{X}{The data frame to prepare}
}
\value{
The same data frame, with only numerical and factor data
}
\description{
Only numerical and factor data is allowed, this function helps prepare the data by converting
any columns that are not numerical or factor to factor
}
\examples{
\donttest{
  data(credit)
  credit <- prim.data.prepare(credit)
}
}
\author{
Jurian Baas
}
