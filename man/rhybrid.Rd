\name{rhybrid}
\alias{rhybrid}
\title{Random Sampler for the "Hybrid Rosenbrock Function"}
\usage{
rhybrid(n,a,b,mu)
}
\description{
Random Sampler for the hybrid rosenbrock function, with two list of vectors parameters b and input x/ xprime.
}
\arguments{
  \item{n}{Sample size.}
  \item{mu}{Density shift.}
  \item{a}{Parameters for xprime.}
  \item{b}{List of parameters corresponding to input blocks.}
}
\value{
  Returns a sample of size \code{n} from the "hybrid" Rosenbrock distribution.
}
\examples{
a <- 1
mu <- - 1

b <- list(c(3,2),
          c(2,2))

n <- 100

rhybrid(n = n,a = a,b = b,mu = mu)
}
