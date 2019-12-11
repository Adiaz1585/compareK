\name{compareidealK1D}
\alias{compareidealK1D}
\title{Sampling Bayesian Spatial Model}
\usage{
compareidealK1D(votes, group, demleader, repleader, nsamples, burn, thin, printevery, varalpprop)
}
\arguments{
\item{votes}{This arguement can take two types of objects:matrix or ord data.}

\item{group}{This arguement takes a vector with length the number of votes.}

\item{demleader}{This is the index position for the democratic party leader.}

\item{repleader}{This is the index position for the republican party leader.}

\item{nsamples}{The number of samples.}

\item{burn}{The amount of samples to burn.}

\item{thin}{The amount of sample skipped before each save.}

\item{printevery}{The rate at which an update is printed.}

\item{varalpprop}{Tuning parameter for the random walk Metropolis Hastings algorithm.}
}
\value{
Returns a list of matrices consisting of the sampled perameters.
}

\references{
reference to the paper.
}

\examples{
## Not run
## this example can be run if package pscl from CRAN is installed.
h111    <- readKH("https://voteview.com/static/data/out/votes/H111_votes.ord")
groups  <- sample(seq(0,16), dim(h111$votes)[2], replace=T)
demlead <- 193
replead <- 321
compare <- (votes=h111, group=groups, demleader=demlead, repleader=replead, nsamples = 30000L, burn = 1000L, thin = 1L, printevery = 250L, varalpprop = 0.06)
}