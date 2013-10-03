\name{teigen.parallel}
\alias{teigen.parallel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
teigen.parallel: Parallelized implementation of the teigen function
}
\description{
Fits multivariate t-distribution mixture models under a clustering paradigm (default) or classification paradigm.
}
\usage{
teigen.parallel(x, Gs = 1:9, numcores = NULL, models = "all", init = "kmeans",
	 scale = TRUE, dfstart = 50, clas = 0, known= NULL, training = NULL,
	 gauss = FALSE, dfupdate = TRUE, eps = c(0.001, 0.1), anneal=NULL,
	 maxit=c(20,1000) ) 
}

\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
A numeric matrix, data frame, or vector (for univariate data) .
}
\item{Gs}{
%%     ~~Describe \code{Gs} here~~
A number or vector indicating the number of groups to fit. Default is 1-9.
}

  \item{numcores}{
%%     ~~Describe \code{numcores} here~~
Scalar argument giving how many cores for the function to utilize. If \code{NULL} (default), then the function discerns the number of cores available and uses all of them.
}

  \item{models}{
%%     ~~Describe \code{models} here~~
A character vector giving the models to fit. See \code{\link{teigen}} details for a comprehensive list of choices.
}
  \item{init}{
%%     ~~Describe \code{init} here~~
A list of initializing classification. See \code{\link{teigen}} for details .
}
  \item{scale}{
%%     ~~Describe \code{scale} here~~
Logical indicating whether or not the function should scale the data. Default is \code{TRUE} and is the prescribed method --- tEIGEN models are not scale invariant.
}
  \item{dfstart}{
%%     ~~Describe \code{dfstart} here~~
The initialized value for the degrees of freedom. The default is 50.
}
  \item{clas}{
%%     ~~Describe \code{clas} here~~
Value between 0-100 indicating the percentage of data taken to be known. Note that a vector of known classifications is needed. See next argument for an alternative. Default is 0 and performs clustering, otherwise the algorithm chooses the training index randomly (and will return it via \code{index}).
}
  \item{training}{
%%     ~~Describe \code{training} here~~
Optional indexing vector for the observations whose classification is taken to be known. 
}
  \item{known}{
%%     ~~Describe \code{known} here~~
A vector of known classifications that can be numeric or character - optional for clustering, necessary for classification. Must be the same length as the number of rows in the data set. If using in a true classification sense, give samples with unknown classification the value \code{NA} within \code{known} (see training example below).
}
  \item{gauss}{
Logical indicating if the algorithm should use the gaussian distribution. If \code{models="mclust"} or \code{"gaussian"} then \code{gauss=TRUE} is forced.
}
  \item{dfupdate}{
%%     ~~Describe \code{dfupdate} here~~
Logical indicating whether or not the degrees of freedom should be estimated. If \code{FALSE}, the value from \code{dfstart} is used instead.
}
  \item{eps}{
%%     ~~Describe \code{eps} here~~
Vector (of size 2) giving tolerance values for the convergence criterion. First value is the tolerance level for iterated M-steps. Second value is tolerance for the EM algorithm: convergence is based on Aitken's acceleration, see cited papers for more information.
}
\item{anneal}{
%%     ~~Describe \code{anneal} here~~
Optional vector giving the deterministic annealing schedule.
}
\item{maxit}{
%%     ~~Describe \code{maxit} here~~
Vector (of size 2) giving maximum iteration number for the iterated M-steps and EM algorithm, respectively. A warning is displayed if either of these maximums are met.
}

}


\details{
%%  ~~ If necessary, more details than the description above ~~
This function is a parallelized wrapper of the \code{\link{teigen}} function. Please refer to details of that function. 
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
An object of class \code{\link{teigen}}.
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
###Note: numcores set to 2 in order to comply
###with CRAN submission policies (set to higher
###number or NULL to automatically use all available cores)

#Clustering old faithful data with tEIGEN
tfaith <- teigen.parallel(faithful, models="UUUU",
numcores=2, Gs=1:3, init="disc")
plot(tfaith)

#Classification with the iris data set via percentage of
#data taken to have known membership
tiris <- teigen.parallel(iris[,-5], numcores=2, models="CUUU",
init="uniform",clas=50,known=iris[,5])
tiris$tab

}
