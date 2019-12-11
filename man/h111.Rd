\name{h111}
\alias{h111}
\docType{data}
\title{
111th House of Representatives Data
}
\description{
This is roll call data on the 111th(2009 - 2011) House of Representatives. It consists of the all motions voted on and vote casted for each individual. This object is of type "rollcall".
}
\usage{data("h111")}
\format{
The format is:

List of 8 variables

$votes     : num [1:447, 1:1647] 9 6 1 6 6 1 0 6 1 6 

$codes     :List of 4

  ..$ yea       : num [1:3] 1 2 3
  
  ..$ nay       : num [1:3] 4 5 6
  
  ..$ notInLegis: num 0
  
  ..$ missing   : num [1:3] 7 8 9
  
 $n         : int 447
 
 $m         : int 1647
 
 $legis.data:'data.frame':	447 obs. of  6 variables:
 
  ..$state     : Factor w/ 50 levels "AK","AL","AR",..: NA 2 2 2 2 2 2 2 2 1...
  
  ..$icpsrState: num [1:447] 99 41 41 41 41 41 41 41 41 81...
  
  ..$cd         : num [1:447] 0 1 2 3 4 5 5 6 7 1...
  
  ..$icpsrLegis: num [1:447] 99911 20300 20900 20301 29701...
  
  ..$party     : Factor w/ 2 levels "D","R": 1 2 1 2 2 1 2 2 1 2...
  
  ..$partyCode : num [1:447] 100 200 100 200 200 100 200 200 100 200...
  
 $vote.data : NULL
 
 $desc      : NULL
 
 $source    : chr "https://voteview.com/static/data/out/votes/H111_votes.ord"
 
 - attr(*, "class")= chr "rollcall"
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
https://voteview.com/static/data/out/votes/H111_votes.ord
}
\examples{
data(h111)
}
\keyword{datasets}