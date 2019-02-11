\name{complement}
\alias{complement}
\title{
reverse complement and complement base of nucleotide sequences
}
\description{
get  reverse complement and complement nucleotide sequences
}
\usage{
complement(ntseq, reverse = FALSE)
}
\arguments{
  \item{ntseq}{
The primer/probe sequence as string or vector of characters
}
  \item{reverse}{
Get reverse complemet sequence or not
}
}
\value{
nucleotide sequences
}
\author{
Junhui Li
}
\examples{
complement("ATCGYCGYsWwsaVv")
complement("ATCGYCGYsWwsaVv",reverse=TRUE)
}