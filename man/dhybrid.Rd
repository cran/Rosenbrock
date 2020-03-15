\name{dhybrid}
\alias{dhybrid}
\title{Density Function "Hybrid Rosenbrock Function"}
\usage{
dhybrid(xprime,x,a,b,mu)
}
\description{
Density function for the hybrid rosenbrock function, with two list of vectors parameters b and input x/ xprime.
}
\arguments{
  \item{xprime}{Central input.}
  \item{x}{List of input vector blocks.}
  \item{mu}{Density shift.}
  \item{a}{Parameters for xprime.}
  \item{b}{List of parameters corresponding to input blocks.}
}
\value{
  Returns the density value of the "hybrid" Rosenbrock distribution at point \code{x}/\code{xprime} for parameters \code{mu,a} and \code{b}.
}
\examples{
a <- 1
xprime <- 1
mu <- - 1

x <- list(c(1,2),
          c(2,3))

b <- list(c(3,2),
          c(2,2))

dhybrid(xprime = xprime,x = x,a = a,b = b,mu = mu)
}
