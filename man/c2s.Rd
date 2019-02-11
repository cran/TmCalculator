\name{c2s}
\alias{c2s}
\title{
convert a vector of characters into a string
}
\description{
Simply convert a vector of characters such as c("H","e","l","l","o","W","o","r","l","d") into a single string such as "HelloWorld".
}
\usage{
c2s(characters)
}
\arguments{
  \item{characters}{
A vector of characters
}
}
\value{
Retrun a strings
}
\references{
\code{citation("TmCalculator")}
}
\author{Junhui Li}
\seealso{
\code{\link{s2c}}
}
\examples{
c2s(c("H","e","l","l","o","W","o","r","l","d"))
}