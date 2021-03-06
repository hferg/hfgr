% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/repData.R
\name{repData}
\alias{repData}
\title{repData}
\usage{
repData(dat, n, prefix)
}
\arguments{
\item{dat}{The data file to be replicated.}

\item{n}{The number of times to write the file to the disk, or a vector of suffixes that the data is replicated using.}

\item{prefix}{The prefix for each of the copies of the data file.}
}
\description{
This function will write a vector of data points to the working directory n times, with a given prefix and the number of the replicate as the suffix. Ideally, the names of the vector are the tip labels of the tree!
}

