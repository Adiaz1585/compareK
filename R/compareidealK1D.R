#' MCMC for Bayesian Spatial Model using Probit link
#'
#' @param path      Path to csv file.
#' @param nsamples      Number of samples to run.
#' @param burn      Number of samples to burn before saving parameter samples.
#' @param thin      Rate to save parameter samples.
compareidealK1D <- function(votes, group, demleader, repleader, nsamples = 30000L, burn = 1000L, thin = 1L, printevery = 250L,  varalpprop = 0.06) {

    #Take the objet votes and create y and missing.  Two different approaches depending on wheters is a roll call object (ord file) or a matrix
    if(is(votes,"rollcall")){
        print("here")
        y       <- as.matrix(votes$votes)
        missing <- as.matrix(votes$votes)
        II      <- dim(y)[1]
        JJ      <- dim(y)[2]
        for(i in     1:II){
            for(j in 1:JJ){
                if(y[i,j]%in%votes$codes$yea){
                    y[i,j]       <- 1
                    missing[i,j] <- 0
                }else{
                    if(y[i,j]%in%votes$codes$nay){
                        y[i,j]       <- 0
                        missing[i,j] <- 0
                    }else{
                        y[i,j]       <- 0
                        missing[i,j] <- 1
                    }
                }
            }
        }
        votes <- y
    }else{
        if(is(votes,"matrix")){
            if(any(!(votes)%in%c(0,1,NA))){
                stop("If votes is a matrix, it can only contain values that are 1 (Yea), 0 (Nay), or NA (Missing/Absent)")
            }
            ###########
            y       <- votes
            missing <- votes
            II      <- dim(y)[1]
            JJ      <- dim(y)[2]
            for(i in 1:II){
                for(j in 1:JJ){
                    if(is.na(y[i,j])){
                        y[i,j]       <- 0
                        missing[i,j] <- 1
                    }else{
                        missing[i,j] <- 0
                    }
                }
            }
        }else{
            stop("The class of the object 'votes' is invalid ")
        }
    }
    group      <- as.matrix(group)
    ngroups    <- length(unique(group))
    nsamples   <- as.integer(nsamples)
    burn       <- as.integer(burn)
    thin       <- as.integer(thin)
    printevery <- as.integer(printevery)
    demleader  <- as.integer(demleader)
    repleader  <- as.integer(repleader)


    if(ngroups!=(max(group)+1))   stop("groups must be labeled continuously starting at 0")
    if(any(dim(missing)!=dim(y))) stop("The dimensions of votes and missing must match")
    if(dim(y)[2]!=length(group))  stop("The length of groups does not match the number of votes")
    if(nsamples <= burn)          stop("nsamples must be greater than burn")
    if(nsamples <= printevery)    stop("nsamples must be greater than printevery")

     .Call(`_compareKgroups`, ngroups, y, missing, demleader, repleader, group, nsamples, burn, thin, printevery,  varalpprop)
   
}

