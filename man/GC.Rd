\name{GC}
\alias{GC}
\title{
Calculate G and C content of nucleotide sequences
}
\description{
Calculate G and C content of nucleotide sequences.The number of G and C in sequence is divided by length of sequence(when totalnt is TRUE) or the number of all A,T,C,G and ambiguous base.
}
\usage{
GC(ntseq, ambiguous = FALSE, totalnt = FALSE)
}
\arguments{
  \item{ntseq}{
The primer/probe sequence as string or vector of characters
}
  \item{ambiguous}{
Ambiguous bases are taken into account to compute the G and C content when ambiguous is TRUE. 
}
  \item{totalnt}{
Sum of 'G' and 'C' bases divided by the length of the sequence when totalnt is TRUE
}
}
\value{
Content of G and C(range from 0 to 100)
}
\author{
Junhui Li
}
\examples{
GC(c("a","t","c","t","g","g","g","c","c","a","g","t","a"))#53.84615
GC("GCATSWSYK",ambiguous = TRUE)#55.55556
}

