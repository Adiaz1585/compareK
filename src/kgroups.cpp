#include <RcppArmadillo.h>
#include "samplers.hpp"
using namespace Rcpp;

//[[Rcpp::depends(RcppArmadillo)]]

//[[Rcpp::export]]
//////////////////////////////////////////////////////////
List compareKgroups(int ngroup, int nsamples, int burn, IntegerMatrix y, IntegerMatrix missing, int demleader, int repleader, IntegerVector group, int thin)
{
    
    
    int II = y.nrow();
    int JJ = y.ncol();
    int KK = ngroup;
    
    int printevery = 1;
    if(nsamples > 10) printevery = nsamples/10;
    
    int nsave = (nsamples - burn)/thin;
    
    Rcout << "I: " << II << " J: " << JJ << " \nIterations: " << nsamples << " Burn: " << burn << " Thin: " << thin << " \nPrint every: " << printevery << " \nDemLeader: " << demleader << " RepLeader: " << repleader<< std::endl;
    
    /* Initialization */
    NumericMatrix z(II, JJ);
    NumericVector mu(JJ);
    NumericVector alpha(JJ);
    NumericMatrix beta(II, KK);
    NumericVector latent(JJ);    //latent parms to sample CRP parm
    IntegerMatrix omega(II, ngroup);
    
    //int group[JJ];
    for(int i=0; i<II; i++){
        for(int k=0; k<KK; k++){
            beta(i, k)=0;
            omega(i, k)=0;  //all same group
        }
        latent[i] = 0.5;
        for(int j=0; j<JJ; j++){
            z(i, j)=0;
        }
    }
    for(int j=0; j<JJ; j++){
        mu[j]=0;
        alpha[j]=0;
        //group[j]=0;
        //if(j>=170) group[j] = 1;
        //if(j>=400) group[j] = 2;    //FOR TESTING K = 3
    }
    
    arma::cube betaOut = arma::zeros<arma::cube>(II, KK, nsave);
    arma::cube omegaOut = arma::zeros<arma::cube>(II, KK, nsave);
    NumericMatrix alphaOut(nsave, JJ);
    NumericMatrix muOut(nsave, JJ);
    
    
    /* Hyper */
    double alp = 1.0;        //crp parm - a priori prob of a new group = alp/(alp+1) = 0.1, ie prob of a different preference among types of bills
    double eta=0.0;
    double kappa2=1.0;
    double w2=1.0;
    double rho=0.0;
    double tau2=1.0;
    double pialpha=0.9;
    
    int index;
    
    
    /*mcmc starts here*/
    for(int m = 0; m<nsamples; m++){
        
        R_CheckUserInterrupt();
        /* Sample! */
        sample_alpha(II,JJ,KK,z,mu,alpha,beta,group,w2,pialpha);
        sample_mu(II,JJ,KK,z,mu,alpha,beta,group,eta,kappa2);
        sample_omega(II,JJ,omega,alp,z,mu,alpha,beta,group,rho,tau2,KK);
        sample_beta(II,JJ,KK,beta,alpha,mu,z,tau2,rho,group,omega);
        
        normalize2(II,JJ,KK,group,beta,alpha,mu,demleader,repleader);
        
        sample_z(II,JJ,KK,y,z,mu,alpha,beta,group);
        
        //sample_latent(II, KK, latent, alp);
        sample_alp(II, KK, latent, omega, &alp);
        
        sample_eta(JJ, &eta, kappa2, mu);
        sample_kappa2(JJ,&kappa2,eta,mu);
        sample_w2(JJ,&w2,alpha);
        sample_rho(II,KK,&rho,tau2,beta,omega);
        sample_tau2(II,KK,&tau2,beta,rho,omega);
        sample_pi(JJ,&pialpha,alpha);
        
        
        sample_missings(II,JJ,KK,y,missing,mu,alpha,beta,group);
        
        R_CheckUserInterrupt();
        
        if(m % printevery == 0){
            Rcout << "\n m: " << m << "\n";
            /* printf("Beta0_98: %f,Beta1_98: %f, Alpha98: %f, Mu98: %f, z1010: %lf \n", beta[98][0],beta[98][1],alpha[98],mu[98], z[10][10]);
             if(KK==3) printf("Beta2_98:  %f\n",beta[98][2]);
             */  Rcout << "beta0_0:  " << beta(0, 0) << " beta1_0:  " << beta(0, 1) << " beta2_0:  " << beta(0, 2) << "\n beta0_1:  " << beta(1, 0) << " beta1_1:  " << beta(1, 1) << " beta2_1: " << beta(1, 2)  << "\n ";
            /*  printf("Z0631: %f, Z0632: %f, Z10: %f, Z11: %f, Z12: %f \n", z[0][631],z[0][632],z[1][0],z[1][1],z[1][2]);
             printf("Eta: %f, Kappa2: %f, W2: %f, \n Rho: %f, Tau2: %f, PiAlpha: %f Alp: %f \n", eta, kappa2, w2, rho, tau2,pialpha, alp);
             */
        }
        if((m % printevery) % (printevery / 10) == 0)
        {
            printf("#+");
        }
        
        
        if(m>=burn && m%thin==0){
            index = (m-burn)/thin;
            for(int i=0; i<II; i++){
                for(int k=0; k<KK; k++){
                    betaOut(i, k, index) = beta(i, k);
                    omegaOut(i, k, index) = omega(i, k);
                    R_CheckUserInterrupt();
                }
            }
            
            alphaOut(index,_) = alpha;
            muOut(index,_) = mu;
            
            R_CheckUserInterrupt();
        }
    }
    
    return List::create(Named("beta") = betaOut, Named("alpha") = alphaOut, Named("mu") = muOut, Named("omega") = omegaOut);
    
}




