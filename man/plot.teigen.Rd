\name{plot.teigen}
\alias{plot.teigen}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
plot.teigen: Plotting Function for tEIGEN Objects
}
\description{
For multivariate data, outputs marginal contour or uncertainty plots to the graphics device for objects of class \code{\link{teigen}}. For univariate data, plot a univariate density plot.
}
\usage{
\method{plot}{teigen}(x, xmarg = 1, ymarg = 2, res = 200, 
	what = c("contour", "uncertainty"),	alpha = 0.4, col = rainbow(x$G),
	pch = 21, cex = NULL, bg = NULL, lty = 1, uncmult = 0,
	levels = c(seq(0.01, 1, by = 0.025), 0.001), 
	main=NULL, xlab=NULL, draw.legend=TRUE, \dots) 
}

\arguments{
  \item{x}{
An object of class \code{\link{teigen}}
}
  \item{xmarg}{
Scalar argument giving the number of the variable to be used on the x-axis
}
 \item{ymarg}{
Scalar argument giving the number of the variable to be used on the y-axis. If \code{NULL}, the teigen object will be interpreted as univariate using x[,xmarg] as the data.
}
 \item{res}{
Scalar argument giving the resolution for the calculation grid required for the contour plot. Default is 200, which results in a 200x200 grid. Also determines how smooth the univariate density curves are (higher res, smoother curves). Ignored for uncertainty plots.
}
 \item{what}{
Only available if the model provided by x is multivariate. Character vector stating which plots should be sent to the graphics device. Choices are \code{"contour"} or \code{"uncertainty"}. Default is to plot both (see Details).
}
  \item{alpha}{
A factor modifying the opacity alpha for the plotted points. Typically in [0,1].
}
 \item{col}{
A specification for the default plotting color. See section 'Color Specification' in  \code{\link{par}}. Note that the number of colors provided must equal to the number of groups in the teigen object (extra colors ignored).
}
 \item{pch}{
Either an integer specifying a symbol or a single character to be used as the default in plotting points. See \code{\link{points}} for possible values and their interpretation. If pch is one of 21:25, see bg for coloring.
}
 \item{cex}{
A numerical value giving the amount by which plotting text and symbols should be magnified relative to the default. For uncertainty plots, cex changes the size of the smallest sized point. The relative sizes amongst the points remains the same. As a result, the sizes of all the points change.
}
 \item{bg}{
Background (fill) color for the open plot symbols if pch is one of 21:25. If pch is in 21:25 point color will be black (col = "black" will be used). If no bg is given to color the inside of the points, col will be used. 
}
 \item{lty}{
The line type for univariate plotting. See \code{\link{par}} for more information. Only updates the group curves, not the density or mixture curves.
}
 \item{uncmult}{
A multiplier for the points on the uncertainty plot. The larger the number, the more the size difference becomes magnified. Large points will get larger faster than smaller points.
}
 \item{levels}{
Numeric vector giving the levels at which contours should be drawn. Default is to draw a contour in 0.25 steps, plus a contour at 0.001. This may result in more/less contours than desired depending on the resulting density.
}
 \item{main}{
  Optional character string for title of plot. Useful default if left as \code{NULL}.
}
 \item{xlab}{
Optional character string for x-axis label.
}
 \item{draw.legend}{
Logical for a default generation of a legend to the right of the plot.
}

\item{\dots}{Options to be passed to \code{plot}.}
}


\details{
%%  ~~ If necessary, more details than the description above ~~
\code{"contour"} plots the marginal distribution of the mixture distribution. For univariate data, or if \code{ymarg} is \code{NULL}, a univariate marginal is provided that includes the kernel density estimate from \code{density()}, the mixture distribution, and colour-coded component densities.

\code{"uncertainty"} plots the uncertainty of each observation's classification - the larger the point, the more uncertainty associated with that observation. Uncertainty in this context refers to the probability that the observation arose from the mixture component specified by the colour in the plot rather than the other components.

The default behavior of the function is to specify both plot types. This opens up an interactive menu from which the user may switch back and forth between both graphs. On exiting the menu, the graph that was plotted last will remain in the open device. 
}
\author{
%%  ~~who you are~~
Jeffrey L. Andrews, Jaymeson R. Wickins, Nicholas M. Boers
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{teigen}}
}

\examples{
set.seed(2521)
tfaith <- teigen(faithful, models = "CCCC", Gs = 1:4, verbose = FALSE)

plot(tfaith, what = "uncertainty", cex = 1.5, uncmult = 1.5)
plot(tfaith, what = "contour")
plot(tfaith, ymarg = NULL, lwd = 2)
}
