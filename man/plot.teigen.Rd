\name{plot.teigen}
\alias{plot.teigen}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
plot.teigen: Plotting Function for tEIGEN Objects
}
\description{
Outputs marginal contour or uncertainty plots to the graphics device for objects of class \code{\link{teigen}}.
}
\usage{
\method{plot}{teigen}(x, xmarg = 1, ymarg = 2, res = 200, levels = c(seq(0.01, 1, by = 0.01), 0.001),
    what = c("contour", "uncertainty"), main=NULL, xlab=NULL, legend=TRUE, \dots) 
}

\arguments{
  \item{x}{
An object of class \code{\link{teigen}}
}
  \item{xmarg}{
Scalar argument giving the number of the variable to be used on the x-axis
}
 \item{ymarg}{
Scalar argument giving the number of the variable to be used on the y-axis
}
 \item{res}{
Scalar argument giving the resolution for the calculation grid required for the contour plot. Default is 200, which results in a 200x200 grid.
}
 \item{levels}{
Numeric vector giving the levels at which contours should be drawn. Default is to draw a contour in 0.01 steps, plus a contour at 0.001. This may result in more/less contours than desired depending on the resulting density.
}
 \item{what}{
Character vector stating which plots should be sent to the graphics device. Choices are \code{"contour"} or \code{"uncertainty"}.
}
 \item{main}{
  Optional character string for title of plot. Useful default if left as \code{NULL}.
}
 \item{xlab}{
Optional character string for x-axis label.
}
 \item{legend}{
Logical for a default generation of a legend in the top left corner of the plot(s).
}

\item{\dots}{Options to be passed to \code{plot}.}
}


\details{
%%  ~~ If necessary, more details than the description above ~~
\code{"contour"} plots the marginal distribution of the mixture distribution. For univariate data, or if \code{ymarg} is \code{NULL}, a univariate marginal is provided that includes the kernel density estimate from \code{density()}, the mixture distribution, and colour-coded component densities.

\code{"uncertainty"} plots the uncertainty of each observation's classification - the larger the point, the more uncertainty associated with that observation. Uncertainty in this context refers to the probability that the observation arose from the mixture component specified by the colour in the plot rather than the other components.
}
\author{
%%  ~~who you are~~
Jeffrey L. Andrews
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{teigen}}
}
