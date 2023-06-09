% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dataSplitAlgorithm.R
\name{DUPLEX}
\alias{DUPLEX}
\title{DUPLEX algorithm}
\usage{
DUPLEX(data, control)
}
\arguments{
\item{data}{Raw data set entered by the user in matrix or data.frame format}

\item{control}{The list of parameters given by the user}
}
\value{
A list of three subsets
}
\description{
The DUPLEX data splitting method was developed by (Snee,1977) based on one of the earliest data splitting algorithms called CADEX or Kennard-Stone sampling (Kennard and Stone, 1969).
DUPLEX draws samples based on Euclidean distances. When applying DUPLEX, the two points which are farthest apart in terms of the Euclidean distance are assigned to the first data set.
The next pair of points that are farthest apart in the remaining list are assigned to the second data set. This process is repeated until both data sets are filled (Snee, 1977).
The original DUPLEX algorithm was used to divide data into two sets. (May et al., 2010) modified the original DUPLEX algorithm to sample data into three data sets in turn based on the proportions of the data sets specified by the user.
Thus, DUPLEX can be used to generate the training, testing, and validation data sets for ANN model development.
}
\references{
Snee, R. D. (1977), Validation of regression models: Methods and examples, Technometrics, 19(4), 415–428.
}
