% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/niceHist.R
\name{niceHist}
\alias{niceHist}
\title{niceHist}
\usage{
niceHist(parameter, fill = "dodgerblue", breaks = "scott", title = "",
  vline = "mean", textsize = 12, face = "bold")
}
\arguments{
\item{parameter}{The a vecotr of values that you want a histogram of.}

\item{fill}{The colour to fill the bars with}

\item{breaks}{The algorithgm to calculate breaks. Defaults to "scott", can also be "fandd"}

\item{title}{The title of the plot, defaults to off.}
}
\value{
A histogram with nice ggplot aesthetics, and good bin widths.
}
\description{
A function to plot a nice histogram - one that uses ggplot2
but has sensible break behaviour in R. Produces a single
histogram but is mostly called internally for things like
plotPosts.
}

