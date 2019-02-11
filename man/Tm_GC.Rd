\name{Tm_GC}
\alias{Tm_GC}
\title{
Calculate the melting temperature using empirical formulas based on GC content
}
\description{
Calculate the melting temperature using empirical formulas based on GC content with differenct rules
}
\usage{
Tm_GC(ntseq, ambiguous = FALSE, userset = NULL, variant = "Primer3Plus",
Na = 50, K = 0, Tris = 0, Mg = 0, dNTPs = 0, saltcorr = 0, mismatch = TRUE)
}
\arguments{
  \item{ntseq}{
The primer/probe sequence as string or vector of characters
}
  \item{ambiguous}{
Ambiguous bases are taken into account to compute the G and C content when ambiguous is TRUE. 
}
  \item{userset}{
A vector of fou coefficientr values. Usersets override valuesets.
}
  \item{variant}{
Empirical constants coefficient with 8 variant: Chester1993, QuikChange, Schildkraut1965, Wetmur1991_MELTING, Wetmur1991_RNA, Wetmur1991_RNA/DNA, Primer3Plus and vonAhsen2001
}
  \item{Na}{
Millimolar concentration of Na
}
  \item{K}{
Millimolar concentration of K
}
  \item{Tris}{
Millimolar concentration of Tris
}
  \item{Mg}{
Millimolar concentration of Mg
}
  \item{dNTPs}{
Millimolar concentration of dNTPs
}
  \item{saltcorr}{
Type of salt correction. Default=5.
}
  \item{mismatch}{
 If 'True' (default) every 'X' in the sequence is counted as mismatch
}
}
\details{
Empirical constants coefficient with 8 variant:
Chester1993: Tm = 69.3 + 0.41(Percentage_GC) - 650/N

QuikChange: Tm = 81.5 + 0.41(Percentage_GC) - 675/N - Percentage_mismatch

Schildkraut1965: Tm = 81.5 + 0.41(Percentage_GC) - 675/N + 16.6 x log[Na+]

Wetmur1991_MELTING: Tm = 81.5 + 0.41(Percentage_GC) - 500/N + 16.6 x log([Na+]/(1.0 + 0.7 x [Na+])) - Percentage_mismatch

Wetmur1991_RNA: Tm = 78 + 0.7(Percentage_GC) - 500/N + 16.6 x log([Na+]/(1.0 + 0.7 x [Na+])) - Percentage_mismatch

Wetmur1991_RNA/DNA: Tm = 67 + 0.8(Percentage_GC) - 500/N + 16.6 x log([Na+]/(1.0 + 0.7 x [Na+])) - Percentage_mismatch

Primer3Plus: Tm = 81.5 + 0.41(Percentage_GC) - 600/N + 16.6 x log[Na+]

vonAhsen2001: Tm = 77.1 + 0.41(Percentage_GC) - 528/N + 11.7 x log[Na+]
}
\references{
Marmur J , Doty P . Determination of the base composition of deoxyribonucleic acid from its thermal denaturation temperature.[J]. Journal of Molecular Biology, 1962, 5(1):109-118.

Schildkraut C . Dependence of the melting temperature of DNA on salt concentration[J]. Biopolymers, 2010, 3(2):195-208.

Wetmur J G . DNA Probes: Applications of the Principles of Nucleic Acid Hybridization[J]. CRC Critical Reviews in Biochemistry, 1991, 26(3-4):33.

Untergasser A , Cutcutache I , Koressaar T , et al. Primer3--new capabilities and interfaces[J]. Nucleic Acids Research, 2012, 40(15):e115-e115.

von Ahsen N, Wittwer CT, Schutz E , et al. Oligonucleotide melting temperatures under PCR conditions: deoxynucleotide Triphosphate and Dimethyl sulfoxide concentrations with comparison to alternative empirical formulas. Clin Chem 2001, 47:1956-1961.
}
\author{
Junhui Li
}
\examples{
ambiguous=TRUE
userset=NULL
variant="Primer3Plus"
Na=50
K=0
Tris=0
Mg=0
dNTPs=0
saltcorr=0
mismatch=TRUE
ntseq <- c("ATCGTGCGTAGCAGTACGATCAGTAG")
Tm_GC(ntseq,ambiguous,userset,variant,Na, K,Tris, Mg, dNTPs, saltcorr, mismatch)
}
