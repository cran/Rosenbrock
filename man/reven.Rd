\name{reven}
\alias{reven}
\title{Random Sampler for the "Even Rosenbrock Function"}
\usage{
reven(n,a,b,mu)
}
\arguments{
  \item{n}{Sample size.}
  \item{mu}{Density shift.}
  \item{a}{Parameters for odd index input.}
  \item{b}{Parameters for even index input.}
}
\description{
Random sampler for the full Rosenbrock function, with two vector parameters a and b.
}
\value{
  Returns a sample of size \code{n} from the "full" Rosenbrock distribution.
}
\examples{
n = 100
a = b = rep(2,2)
mu = rep(1,2)

reven(n = n,mu = mu,a = a,b = b)
}
