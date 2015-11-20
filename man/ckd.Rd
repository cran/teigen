\name{ckd}
\alias{ckd}
\docType{data}
\title{
Indian Chronic Kidney Disease Data
%%   ~~ data name/kind ... ~~
}
\description{
This is a cleaned up version of the Chronic Kidney Disease data set available from the UCI learning repository:

http://archive.ics.uci.edu/ml/datasets/Chronic_Kidney_Disease

Nominal variables have been removed and rows with missing values recorded on the remaining variables have also been removed.
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("ckd")}
\format{
  This data frame contains 203 rows (observations) and 13 columns (variables):
  1.) ckdclass: There are 2 classes, ckd or notckd
  2.) age: in years
  3.) blood.pressure: in mm/Hg
  4.) blood.glucose.random:  in mgs/dl 
  5.) blood.urea: in mgs/dl
  6.) serum.creatinine: in mgs/dl 
  7.) sodium: in mEq/L
  8.) potassium: in mEq/L
  9.) hemoglobin: in gms
  10.) packed.cell.volume
  11.) white.blood.cell.count: in cells/cmm
  12.) red.blood.cell.count: in cells/cmm
}
\source{
See http://archive.ics.uci.edu/ml/datasets/Chronic_Kidney_Disease for original source.
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
Lichman, M. (2013). UCI Machine Learning Repository [http://archive.ics.uci.edu/ml]. Irvine, CA: University of California, School of Information and Computer Science.
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(ckd)
hclustres <- cutree(hclust(dist(ckd[,-1])),3)
table(ckd[,1], hclustres)
}
\keyword{datasets}
