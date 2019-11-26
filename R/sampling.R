#' MCMC for Bayesian Spatial Model using Probit link
#'
#' @param path      Path to csv file.
#' @param nsamples      Number of samples to run.
#' @param burn      Number of samples to burn before saving parameter samples.
#' @param thin      Rate to save parameter samples.
sampling <- function(ngroup, y, missing, group, nsamples = 30000L, burn = 1000L, thin = 1L, demleader, repleader) {
    
    missing = as.matrix(missing)
    y = as.matrix(y)
    group = as.matrix(group)
    
    if(nsamples <=10) stop("The number samples should be greater than 10")
    if(nsamples - burn < 10) stop("nsamples - burn must be 10 or greater")
    
    #Ord file manipulation will go here
    

     .Call(`_compareKgroups`, ngroup, nsamples, burn, y, missing, demleader, repleader, group, thin)
   
}

