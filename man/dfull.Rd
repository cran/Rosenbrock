\name{dfull}
\alias{dfull}
\title{Density Function "Full Rosenbrock Function"}
\usage{
dfull(x,a,b,mu)
}
\description{
Density function for the full rosenbrock function, with two vector parameters a and b.
}
\arguments{
  \item{x}{Input vector.}
  \item{mu}{Density shift.}
  \item{a}{First parameters.}
  \item{b}{Second parameters.}
}
\value{
  Returns the density value of the "full" Rosenbrock distribution at point \code{x} for parameters \code{mu,a} and \code{b}.
}
\examples{
x = rep(1,4)
a = b = rep(2,3)
mu = 1

dfull(x = x,mu = mu,a = a,b = b)
}
