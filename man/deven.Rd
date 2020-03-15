\name{deven}
\alias{deven}
\title{Density Function "Even Rosenbrock Function"}
\usage{
deven(x,a,b,mu)
}
\arguments{
  \item{x}{Input vector.}
  \item{mu}{Density shift.}
  \item{a}{Parameters for odd index input.}
  \item{b}{Parameters for even index input.}
}
\value{
  Returns the density value of the "even" Rosenbrock distribution at point \code{x} for parameters \code{mu,a} and \code{b}.
}
\description{
Density function for the full rosenbrock function, with two vector parameters a and b.
}
\examples{
x = rep(1,4)
a = b = rep(2,2)
mu = rep(1,2)

deven(x = x,a = a,b = b,mu = mu)
}
