\name{predict.teigen}
\alias{predict.teigen}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
predict.teigen: Predicting Function for tEIGEN Objects
}
\description{
Provides the fuzzy probability matrix and classification vector for inputted observations assuming the model provided by the \code{\link{teigen}} object.
}
\usage{
\method{predict}{teigen}(object, newdata=NULL, modelselect="BIC", \dots) 
}

\arguments{
  \item{object}{
An object of class \code{\link{teigen}}
}
  \item{newdata}{
Data frame or matrix of new observations on the same variables used in the fitting of the \code{\link{teigen}} object. For predicting one observation, a vector is permitted. If \code{NULL}, then the observations used in the fitting of the \code{\link{teigen}} object are inputted.
}
 \item{modelselect}{
A character string of either \code{"BIC"} (default) or \code{"ICL"} indicating the desired model-selection criteria to apply to the \code{\link{teigen}} object.
}
\item{\dots}{Arguments to be passed to other functions.}
}


\details{
%%  ~~ If necessary, more details than the description above ~~
Note that the scale argument from the \code{\link{teigen}} object is passed along to the \code{predict}function. See examples below for plotting.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
\item{fuzzy}{Matrix of fuzzy classification probabilities}
\item{classification}{Vector of maximum a posteriori classifications}

}
\author{
%%  ~~who you are~~
Jeffrey L. Andrews
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{teigen}}
}

\examples{
set.seed(2521)
ind <- sample(1:nrow(faithful), 20)
set.seed(256)
tfaith_unscaled <- teigen(faithful[-ind,], models = "UUUU", Gs = 2, verbose = FALSE, scale=FALSE)
pred_unscaled <- predict(tfaith_unscaled, faithful[ind,])
set.seed(256)
tfaith_scaled <- teigen(faithful[-ind,], models = "UUUU", Gs = 2, verbose = FALSE, scale=TRUE)
pred_scaled <- predict(tfaith_scaled, faithful[ind,])
identical(pred_unscaled$classification, pred_scaled$classification)

##Plotting UNSCALED
plot(tfaith_unscaled, what="contour")
points(faithful[ind,1], faithful[ind,2], pch=15)
plotcolours <- rainbow(tfaith_unscaled$G)
points(faithful[ind,1], faithful[ind,2], pch=20, col=plotcolours[pred_unscaled$classification])

##Plotting SCALED
plot(tfaith_scaled, what="contour")
points((faithful[ind,1]-tfaith_scaled$info$scalemeans[1])/tfaith_scaled$info$scalesd[1],
       (faithful[ind,2]-tfaith_scaled$info$scalemeans[2])/tfaith_scaled$info$scalesd[2],
       pch=15)
plotcolours <- rainbow(tfaith_scaled$G)
points((faithful[ind,1]-tfaith_scaled$info$scalemeans[1])/tfaith_scaled$info$scalesd[1],
       (faithful[ind,2]-tfaith_scaled$info$scalemeans[2])/tfaith_scaled$info$scalesd[2],
       pch=20, col=plotcolours[pred_scaled$classification])

}
