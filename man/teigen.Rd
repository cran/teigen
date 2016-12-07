\name{teigen}
\alias{teigen}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
teigen: Function for Model-Based Clustering and Classification with the Multivariate t Distribution
}
\description{
Fits multivariate t-distribution mixture models (with eigen-decomposed covariance structure) to the given data within a clustering paradigm (default) or classification paradigm (by giving either training index or percentage of data taken to be known). Can be run in parallel.
}
\usage{
teigen(x, Gs = 1:9, models = "all", init = "kmeans", scale = TRUE, dfstart = 50, 
	known = NULL, training = NULL, gauss = FALSE, dfupdate = "approx", 
	eps = c(0.001, 0.1), verbose = TRUE, maxit = c(Inf,Inf), 
	convstyle = "aitkens", parallel.cores = FALSE, 
	ememargs = list(25, 5, "UUUU", "hard"))}
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
A character vector giving the models to fit. See details for a comprehensive list of choices.
}
  \item{init}{
%%     ~~Describe \code{init} here~~
A list of initializing classification of the form that \code{init[[G]]} contains the initializing vector for all G considered (see example below). Alternatively, the user can use a character string indicating initialization method. Currently the user can choose from \code{"kmeans"} (default), `hard' random - \code{"hard"}, `soft' random - \code{"soft"},  \code{``emem''} (see \code{ememargs} below for description), and \code{"uniform"} (classification only).
}
  \item{scale}{
%%     ~~Describe \code{scale} here~~
Logical indicating whether or not the function should scale the data. Default is \code{TRUE} and is the prescribed method --- tEIGEN models are not scale invariant.
}
  \item{dfstart}{
%%     ~~Describe \code{dfstart} here~~
The initialized value for the degrees of freedom. The default is 50.
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
Character string or logical indicating how the degrees of freedom should be estimated. The default is \code{"approx"} indicating a closed form approximation be used. Alternatively, \code{"numeric"} can be specified which makes use of \code{\link{uniroot}}. If \code{FALSE}, the value from \code{dfstart} is used and the degrees of freedom are not updated. If \code{TRUE}, \code{"numeric"} will be used for back-compatibility.
}
  \item{eps}{
%%     ~~Describe \code{eps} here~~
Vector (of size 2) giving tolerance values for the convergence criterion. First value is the tolerance level for iterated M-steps. Second value is tolerance for the EM algorithm: convergence is based on Aitken's acceleration, see cited papers for more information.
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
Logical indicating whether the running output should be displayed. This option is not available in parallel. What is displayed depends on the width of the R window. With a width of 80 or larger: time run, estimated time remaining, percent complete are all displayed.
}
\item{maxit}{
%%     ~~Describe \code{maxit} here~~
Vector (of size 2) giving maximum iteration number for the iterated M-steps and EM algorithm, respectively. A warning is displayed if either of these maximums are met, default for both is Inf (aka, no limit).
}
  \item{convstyle}{
%%     ~~Describe \code{convstyle} here~~
Character string specifying the method of determining convergence. Default is "aitkens" which uses a criteria based on Aitken's acceleration, but "lop" (lack of progress) may be used instead.
}
  \item{parallel.cores}{
%%     ~~Describe \code{parallel.cores} here~~
Logical indicating whether to run teigen in parallel or not. If \code{TRUE}, then the function determines the number of cores available and uses all of them. Alternatively, a positive integer may be provided indicating the number of cores the user wishes to use for running in parallel. 
}
  \item{ememargs}{
    A list of the controls for the emEM initialization with named elements: 
      \code{numstart} - numeric, number of starts (default 25)
      \code{iter} - numeric, number of EM iterations (default 5)
      \code{model} - character string, model name to be used (default ``UUUU'' from C,U,I nomenclature...see details below) 
      \code{init} - character string, initialization method for emEM (default \code{hard}, or \code{soft}, or \code{kmeans}).
    The emEM initialiazation will run multiple, randomized initialization attempts for a limited number of iterations, and then continue the model-fitting process.
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
Model specification (via the \code{models} argument) follows either the nomenclature discussed in Andrews and McNicholas (2012), or via the nomenclature popularized in other packages. In both cases, the nomenclature refers to the decomposition and constraints on the covariance matrix: \deqn{\Sigma_g = \lambda_g D_g A_g D_g'}

The nomenclature from Andrews and McNicholas (2012) gives four letters, each letter referring to (in order) \eqn{\lambda}, D, A, and the degrees of freedom. Possible letters are \code{"U"} for unconstrained, \code{"C"} for constrained (across groups), and \code{"I"} for when the parameter is replaced by the appropriately sized identity matrix (where applicable). As an example, the string \code{"UICC"} would refer to the model where \eqn{\Sigma_g = \lambda_g A} with degrees of freedom held equal across groups.

The alternative nomenclature describes (in order) the volume (\eqn{\lambda}), shape (A), orientation (D), and degrees of freedom in terms of \code{"V"}ariable, \code{"E"}qual, or the \code{"I"}dentity matrix. The example model discussed in the previous paragraph would then be called by \code{"VEIE"}.

Possible univariate models are \code{c("univUU", "univUC", "univCU", "univCC")} where the first capital letter describes \code{"U"}nconstrained or \code{"C"}onstrained variance and the second capital letter refers to the degrees of freedom. Once again, \code{"V"}ariable or \code{"E"}qual can replace U and C, but this time the orders match between the nomenclatures.

As many models as desired can be selected and ran via the vector supplied to \code{models}. More commonly, subsets can be called by the following character strings: 
	\code{"all"} runs all 28 tEIGEN models (default),
	\code{"dfunconstrained"} runs the 14 unconstrained degrees of freedom models, 
	\code{"dfconstrained"} runs the 14 constrained degrees of freedom models,  
	\code{"mclust"} runs the 10 MCLUST models using the multivariate Gaussian distribution rather than the multivariate t,
	\code{"gaussian"} is similar but includes four further mixture models than MCLUST,
	\code{"univariate"} runs the univariate models - will automatically be called if one of the previous shortcuts is used on univariate data. 
	
Note that adding \code{"alt"} to the beginning of those previously mentioned characters strings will run the same models, but return results with the V-E-I nomenclature.

Also note that for G=1, several models are equivalent (for example, UUUU and CCCC). Thus, for G=1 only one model from each set of equivalent models will be run.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
\item{x}{Data used for clustering/classification.}

\item{index}{Indexing vector giving observations taken to be known (only available when classification is performed). }

\item{classification}{Vector of group classifications as determined by the BIC.}

\item{bic}{BIC of the best fitted model.}

\item{modelname}{Name of the best model according to the BIC.}

\item{allbic}{Matrix of BIC values according to model and G. A value of -Inf is returned when the model did not converge. }

\item{bestmodel}{Character string giving best model (BIC) details.}

\item{G}{Value corresponding to the number of components chosen by the BIC. }

\item{tab}{Classification table for BIC-selected model (only available when \code{known} is given). When classification is used the ``known'' observations are left out of the table. }

\item{fuzzy}{The fuzzy clustering matrix for the model selected by the BIC. }

\item{logl}{The log-likelihood corresponding to the model with the best BIC. }

\item{iter}{The number of iterations until convergence for the model selected by the BIC.}

\item{parameters}{List containing the fitted parameters: \code{mean} - matrix of means where the rows correspond to the component and the columns are the variables; \code{sigma} - array of covariance matrices (multivariate) or variances (univariate); \code{lambda} - vector of scale parameters, or constants of proportionality; \code{d} - eigenvectors, or orientation matrices; \code{a} - diagonal matrix proportional to eigenvalues, or shape matrices; \code{df} - vector containing the degrees of freedom for each component; \code{weights} - matrix of the expected value of the characteristic weights; \code{pig} - a vector giving the mixing proportions. }

\item{iclresults}{List containing all the previous outputs, except \code{x} and \code{index}, pertaining to the model chosen by the best ICL (all under the same name except \code{allicl} and \code{icl} are the equivalent of \code{allbic} and \code{bic}, respectively).}

\item{info}{List containing a few of the original user inputs, for use by other dedicated functions of the \code{teigen} class.}
}
\references{
Andrews JL and McNicholas PD. ``Model-based clustering, classification, and discriminant analysis with the multivariate \emph{t}-distribution: The \emph{t}EIGEN family'' \emph{Statistics and Computing} 22(5), 1021--1029. 

Andrews JL, McNicholas PD, and Subedi S (2011) ``Model-based classification via mixtures of multivariate t-distributions'' \emph{Computational Statistics and Data Analysis} 55, 520--529.
}
\author{
%%  ~~who you are~~
Jeffrey L. Andrews, Jaymeson R. Wickins, Nicholas M. Boers, Paul D. McNicholas
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See package manual \code{\link{tEIGEN}}
}
\examples{
###Note that only one model is run for each example
###in order to reduce computation time

#Clustering old faithful data with hard random start
tfaith <- teigen(faithful, models="UUUU", Gs=1:3, init="hard")
plot(tfaith, what = "uncertainty")
summary(tfaith)

#Clustering old faithful with hierarchical starting values
initial_list <- list()
clustree <- hclust(dist(faithful))
for(i in 1:3){
	initial_list[[i]] <- cutree(clustree,i)
}
tfaith <- teigen(faithful, models="CUCU", Gs=1:3, init=initial_list)
print(tfaith)

#Classification with the iris data set
#Introducing NAs is not required; this is to illustrate a `true' classification scenario
irisknown <- iris[,5]
irisknown[134:150] <- NA
triris <- teigen(iris[,-5], models="CUUU", init="uniform", known=irisknown)

##Parallel examples:
###Note: parallel.cores set to 2 in order to comply
###with CRAN submission policies (set to higher
###number or TRUE to automatically use all available cores)

#Clustering old faithful data with tEIGEN
tfaith <- teigen(faithful, models="UUUU",
parallel.cores=2, Gs=1:3, init="hard")
plot(tfaith, what = "contour")

#Classification with the iris data set 
irisknown <- iris[,5]
irisknown[sample(1:nrow(iris),50)] <- NA
tiris <- teigen(iris[,-5], parallel.cores=2, models="CUUU",
init="uniform", known=irisknown)
tiris$tab
}


