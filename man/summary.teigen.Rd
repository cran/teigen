\name{summary.teigen}
\alias{summary.teigen}
\alias{print.summary.teigen}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
summary.teigen: Summary Function for tEIGEN Objects
}
\description{
Summary method for class "teigen". Gives summary information.
}
\usage{
\method{summary}{teigen}(object, \dots)
\method{print}{summary.teigen}(x, \dots)
}

\arguments{
  \item{object}{
An object of class \code{\link{teigen}}
}
  \item{x}{
An object of class "\code{summary.teigen}".
}
\item{\dots}{
Options to be passed to \code{summary}.
}
}

\value{
An object of class "\code{summary.teigen}" that has a specialized \code{print} method. The object is a list containing the BIC and ICL values, as well as loglik value, model number and group number for the BIC. These values are also stored for ICL if it disagrees with the BIC value. 
}
\author{
%%  ~~who you are~~
Jaymeson R. Wickins, Nicholas M. Boers, Jeffrey L. Andrews
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{teigen}}
}
