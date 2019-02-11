\name{chem_correction}
\alias{chem_correction}
\title{
Correct melting temperature with chemical substances
}
\description{
Correct a given melting temperature with DMSO and formamide, and these corrections are rough approximations.
}
\usage{
chem_correction(Tm, DMSO = 0, fmd = 0, DMSOfactor = 0.75, 
fmdfactor = 0.65, fmdmethod = "concentration", ptGC = NULL)
}
\arguments{
  \item{Tm}{
Melting temperature
}
  \item{DMSO}{
Percent DMSO
}
  \item{fmd}{
Formamide concentration in percentage (fmdmethod="concentration") or molar (fmdmethod="molar").
}
  \item{DMSOfactor}{
Coeffecient of Tm decreases per percent DMSO. Default=0.75 von Ahsen N (2001) <PMID:11673362>. Other published values are 0.5, 0.6 and 0.675.
}
  \item{fmdfactor}{
Coeffecient of Tm decrease per percent formamide. Default=0.65. Several papers report factors between 0.6 and 0.72.
}
  \item{fmdmethod}{
"concentration" method for formamide concentration in percentage and "molar" for formamide concentration in molar.
}
  \item{ptGC}{
Content of GC
}
}
\details{
fmdmethod = "concentration"
Tm = Tm - factor*percentage_of_formamide

fmdmethod = "molar"
Tm = Tm + (0.453(f(GC)) - 2.88) x [formamide]
}
\references{
von Ahsen N, Wittwer CT, Schutz E , et al. Oligonucleotide melting temperatures under PCR conditions: deoxynucleotide Triphosphate and Dimethyl sulfoxide concentrations with comparison to alternative empirical formulas. Clin Chem 2001, 47:1956-C1961.
}
\author{
Junhui Li
}
\examples{
chem_correction(70, DMSO=3) #67.75
chem_correction(70, fmd=5) #66.75
chem_correction(70, fmdmethod="molar", fmd=1.25,ptGC=50) #66.68
}
