\name{check_filter}
\alias{check_filter}
\title{
Check and filter valid base of nucleotide sequences
}
\description{
 In general, whitespaces and non-base characters are removed and characters are converted to uppercase in given method.
}
\usage{
check_filter(ntseq, method)
}
\arguments{
  \item{ntseq}{
The primer/probe sequence as string or vector of characters
}
  \item{method}{
TM_Wallace: check and return "A","B","C","D","G","H","I","K","M","N","R","S","T","V","W" and "Y"

TM_GC: check and return "A","B","C","D","G","H","I","K","M","N","R","S","T","V","W", "X" and "Y"

TM_NN: check and return "A","C","G","I" and "T"
}
}
\value{
Return a sequence which fullfils the requirements of the given method.
}
\author{
Junhui Li
}
\examples{
ntseq <- c("ATCGBDHKMNRVYWSqq")
check_filter(ntseq,method='Tm_Wallace')
check_filter(ntseq,method='Tm_NN')
}