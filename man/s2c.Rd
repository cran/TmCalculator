\name{s2c}
\alias{s2c}
\title{
convert a string into a vector of characters
}
\description{
Simply convert  a single string such as "HelloWorld" into a vector of characters such as c("H","e","l","l","o","W","o","r","l","d")
}
\usage{
s2c(strings)
}
\arguments{
  \item{strings}{
A single string such as "HelloWorld"
}
}
\value{
Retrun a vector of characters
}
\references{
\code{citation("TmCalculator")}
}
\author{Junhui Li}
\seealso{
\code{\link{c2s}}
}
\examples{
s2c(c("HelloWorld"))
}