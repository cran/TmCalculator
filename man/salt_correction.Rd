\name{salt_correction}
\alias{salt_correction}
\title{
Calculate the corrected melting temperature with salt ions
}
\description{
It will correct melting temperature or entropy based different operations
}
\usage{
salt_correction(Na = 0, K = 0, Tris = 0, Mg = 0, dNTPs = 0,
method = 1, ntseq = NULL, ambiguous = FALSE)
}
\arguments{
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
  \item{method}{
Which method to be applied. Methods 1-4 correct Tm, method 5 corrects deltaS, methods 6 and 7 correct 1/Tm.
}
  \item{ntseq}{
The primer/probe sequence as string or vector of characters
}
  \item{ambiguous}{
Ambiguous bases are taken into account to compute the G and C content when ambiguous is TRUE. 
}
}
\details{
methods 1-4: Tm(new) = Tm(old) + corr

method 5: deltaS(new) = deltaS(old) + corr

methods 6+7: Tm(new) = 1/(1/Tm(old) + corr)

The methods are:

1 Schildkraut C (2010) <doi:10.1002/bip.360030207>

2 Wetmur J G (1991) <doi:10.3109/10409239109114069>

3 SantaLucia J (1996) <doi:10.1021/bi951907q>

4 SantaLucia J (1998) <doi:10.1073/pnas.95.4.1460>

5 SantaLucia J (1998) <doi:10.1073/pnas.95.4.1460>

6 Owczarzy R (2004) <doi:10.1021/bi034621r>

7 Owczarzy R (2008) <doi:10.1021/bi702363u>
}

\references{
Schildkraut C . Dependence of the melting temperature of DNA on salt concentration[J]. Biopolymers, 2010, 3(2):195-208.

Wetmur J G . DNA Probes: Applications of the Principles of Nucleic Acid Hybridization[J]. CRC Critical Reviews in Biochemistry, 1991, 26(3-4):3

Santalucia , J , Allawi H T , Seneviratne P A . Improved Nearest-Neighbor Parameters for Predicting DNA Duplex Stability, [J]. Biochemistry, 1996, 35(11):3555-3562.

SantaLucia, J. A unified view of polymer, dumbbell, and oligonucleotide DNA nearest-neighbor thermodynamics[J]. Proceedings of the National Academy of Sciences, 1998, 95(4):1460-1465.

Owczarzy R , You Y , Moreira B G , et al. Effects of Sodium Ions on DNA Duplex Oligomers: Improved Predictions ofMelting Temperatures[J]. Biochemistry, 2004, 43(12):3537-3554.

Owczarzy R , Moreira B G , You Y , et al. Predicting Stability of DNA Duplexes in Solutions Containing Magnesium and Monovalent Cations[J]. Biochemistry, 2008, 47(19):5336-5353.
}
\author{
Junhui Li
}
\examples{
ntseq <- c("acgtTGCAATGCCGTAWSDBSYXX")
salt_correction(Na=390, K=20, Tris=0, Mg=10, dNTPs=25, method=7, ntseq)
mySeq <- c("A", "C", "G", "T", "T", "G", "C", "A", "A", "T", "G", 
"C", "C", "G", "T", "A", "W", "S", "D", "B","S", "Y", "X", "X")
salt_correction(Na=390, K=20, Tris=0, Mg=10, dNTPs=25, method=7, mySeq)
}