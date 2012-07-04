\name{teigen}
\alias{teigen}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
teigen: Model-based clustering and classification with the multivariate t-distribution
}
\description{
Fits multivariate t-distribution mixture models (with eigen-decomposed covariance structure) to the given data underneath a clustering paradigm (default) or classification paradigm (by giving either training index or percentage of data taken to be known).
}
\usage{
teigen(x, Gs = 1:9, models = "all", init ="kmeans", scale = TRUE, dfstart = 50, 
      clas = 0, known= NULL, training = NULL, gauss = FALSE, dfupdate = TRUE, eps = 0.1)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
A numeric matrix, data frame, or vector (for univariate data) .
}
  \item{Gs}{
%%     ~~Describe \code{Gs} here~~
A number or vector indicating the number of groups to fit. Default is 1-9.
}
  \item{models}{
%%     ~~Describe \code{models} here~~
A character vector giving the models to fit. Alternatively, \code{"all"} runs all tEIGEN models (default), \code{"dfunconstrained"} runs all the unconstrained degrees of freedom models, \code{"dfconstrained"} runs all the constrained degrees of freedom models,  \code{"mclust"} approximates MCLUST models and restricts degrees of freedom estimation to approximately normal, \code{"gaussian"} is similar but includes two further mixture models than MCLUST, and finally \code{"univariate"} runs the univariate models.
}
  \item{init}{
%%     ~~Describe \code{init} here~~
A list of initializing classification of the form that \code{init[[G]]} contains the initializing vector for all G considered. Alternatively, the user can use a character string indicating initialization method. Currently the user can choose from \code{"kmeans"} (default), `hard' random - \code{"hard"}, `soft' random - \code{"soft"}, and \code{"uniform"} (classification only).
}
  \item{scale}{
%%     ~~Describe \code{scale} here~~
Logical indicating whether or not the function should scale the data. Default is \code{TRUE} and is the prescribed method.
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
Logical indicating if the algorithm should approximate the gaussian distribution. If \code{models="mclust"} then \code{gauss=TRUE} is forced.
}
  \item{dfupdate}{
%%     ~~Describe \code{dfupdate} here~~
Logical indicating whether or not the degrees of freedom should be estimated. If \code{FALSE}, the value from \code{dfstart} is used instead.
}
  \item{eps}{
%%     ~~Describe \code{eps} here~~
This value is the tolerance for the convergence criterion. Convergence is based on Aitken's acceleration; see cited papers for more information.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
Note that for G=1, several models are equivalent (for example, UUUU and CCCC). Thus, for G=1 only one model from each set of equivalent models will be run.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
\item{x}{Data used for clustering/classification.}
\item{index}{Indexing vector giving observations taken to be known (only available when \code{clas} is set greater than 0 or \code{training} is given). }
\item{classification}{Vector of group classifications as determined by the BIC.}
\item{allbic}{Matrix of BIC values according to model and G. A value of -Inf is returned when the model did not converge. }
\item{bestmodel}{Character string giving best model (BIC) details.}
\item{G}{Value corresponding to the number of components chosen by the BIC. }
\item{tab}{Classification table for BIC model (only available when \code{known} is given). When classification is used the ``known'' observations are left out of the table. }
\item{fuzzy}{The fuzzy clustering matrix for the model selected by the BIC. }
\item{logl}{The log-likelihood corresponding to the model with the best BIC. }
\item{parameters}{List containing the fitted parameters: \code{mean} - matrix of means where the rows correspond to the component and the columns are the variables; \code{sigma} - array of covariance matrices (multivariate) or variances (univariate); \code{lambda} - vector of scale parameters, or constants of proportionality; \code{d} - eigenvectors, or orientation matrices; \code{a} - diagonal matrix proportional to eigenvalues, or shape matrices; \code{df} - vector containing the degrees of freedom for each component. }
\item{iclresults}{List containing all the previous outputs, except \code{x} and \code{index}, pertaining to the model chosen by the best ICL (all under the same name except \code{allicl} and \code{icl} are the equivalent of \code{allbic} and \code{bic}, respectively).}
}
\references{
Andrews JL and McNicholas PD. ``Model-based clustering, classification, and discriminant analysis with the multivariate \emph{t}-distribution: The \emph{t}EIGEN family'' \emph{Statistics and Computing}. To appear: http://dx.doi.org/10.1007/s11222-011-9272-x

Andrews JL, McNicholas PD, and Subedi S (2011) ``Model-based classification via mixtures of multivariate t-distributions'' \emph{Computational Statistics and Data Analysis} 55, 520--529.
}
\author{
%%  ~~who you are~~
Jeffrey L. Andrews \code{<andrewsj@uoguelph.ca>} 
with contributions from Paul D. McNicholas
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See package manual \code{\link{tEIGEN}}
}
\examples{

data(faithful)

#Clustering old faithful data with tEIGEN
tfaith <- teigen(faithful, Gs=1:3, init="disc")
plot(faithful, col=c("blue","green","red")[tfaith$clas])

#Clustering old faithful data (MCLUST approximation)
mfaith <- teigen(faithful, Gs=1:3, init="disc", models="mclust")
plot(faithful, col=c("blue","green","red")[mfaith$clas])

#Classification with the iris data set via percentage of data taken to have known membership
tiris <- teigen(iris[,-5],init="uniform",clas=50,known=iris[,5])
tiris$tab

#Classification with the iris data set via training set
irisknown <- iris[,5]
#Introducing NAs is not required; this is to illustrate a `true' classification scenario
irisknown[134:150] <- NA
triris <- teigen(iris[,-5],init="uniform",known=irisknown,training=1:133)
}