\name{Tm_Wallace}
\alias{Tm_Wallace}
\title{
Calculate the melting temperature using the 'Wallace rule'
}
\description{
The Wallace rule is often used as rule of thumb for approximate melting temperature calculations for primers with 14 to 20 nt length.
}
\usage{
Tm_Wallace(ntseq, ambiguous = FALSE)
}
\arguments{
  \item{ntseq}{
The primer/probe sequence as string or vector of characters
}
  \item{ambiguous}{
Ambiguous bases are taken into account to compute the G and C content when ambiguous is TRUE. 
}
}
\value{
a numeric melting temperature
}
\references{
Thein S L , Lynch J R , Weatherall D J , et al. DIRECT DETECTION OF HAEMOGLOBIN E WITH SYNTHETIC OLIGONUCLEOTIDES[J]. The Lancet, 1986, 327(8472):93.
}
\author{
Junhui Li
}
\examples{
ntseq = c('acgtTGCAATGCCGTAWSDBSY') #for wallace rule
Tm_Wallace(ntseq,ambiguous = TRUE)
Tm_Wallace(ntseq,ambiguous = FALSE)
}
