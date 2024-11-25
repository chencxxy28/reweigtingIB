##Varying coefficient and GLM


library(rootSolve)
library(MASS)
library(gam)
library(np)
library(dplyr)
library(glmnet)
library(ncvreg)
library(kableExtra)
library(ggplot2)
library(boot)
###code for information borrowing learner
#estimating function in the constrain
##this is the h_i estimating function. input has theta.
ee_ib<-function(theta,y,xm)
{ 
  n<-length(y)
  ee<-rep()
  xm.int<-cbind(1,xm)
  for(i in 1:n)
  {
    ee_i<-xm.int[i,]*((y-inv.logit(xm.int%*%theta))[i])
    ee<-cbind(ee,ee_i)
  }
  return(ee)
}
##first derivative of -log EL
R1der<-function(lambda,ZZ)
{
  apply(ZZ,2,function(xx)
  {as.matrix(xx,ncol=1)/as.vector((1+t(lambda)%*%as.matrix(xx,ncol=1)))})%*%rep(1,ncol(ZZ))
}

##second derivative of -log EL
R2der<-function(lambda,ZZ)
{
  r2der<-0
  for(i in 1:ncol(ZZ))
  {
    r2der_i<--as.matrix(ZZ[,i],ncol=1)%*%t(as.matrix(ZZ[,i],ncol=1))/as.vector(1+t(lambda)%*%as.matrix(ZZ[,i],ncol=1))^2
    r2der<-r2der+r2der_i
  }
  r2der
}
##-log EL
R0der<-function(lambda,ZZ)
{
  apply(ZZ,2, function (xx) {log(as.vector(1+t(lambda)%*%as.matrix(xx,ncol=1)))})%*%rep(1,ncol(ZZ))
}
#function to find lambda based on empirical likelihood!!!!!Maybe we need some change.!!!!!!!!!!!!!!!!!!!
lambda_find<-function(theta,y,xm)
{
  ZZ<-ee_ib(theta,y,xm)
  dim(ZZ)
  apply(ZZ,1,mean)
  
  gamma<-1
  c<-0
  lambda<-rep(0,nrow(ZZ))
  tol<-10e-8
  ssss=0
  
  repeat{
    ssss=ssss+1
    rl<-R1der(lambda,ZZ)
    rll<-R2der(lambda,ZZ) 
    Delta<--ginv(rll)%*%rl
    if(mean(abs(Delta))<tol|ssss>100)
    {break}else{
      repeat{
        mm<-0
        repeat{
          delta<-gamma*Delta
          index_1<-apply(ZZ,2,function (xx)
          {ifelse(1+t(lambda+delta)%*%as.matrix(xx,ncol=1)<=0,1,0)}
          )
          if (sum(index_1)>0)
          {gamma<-gamma/2
          mm<-mm+1}else{break}}
        index_2<-ifelse(R0der(lambda+delta,ZZ)-R0der(lambda,ZZ)<0,1,0)
        if (index_2==1)
        {gamma<-gamma/2}else{break}
      }
    }
    lambda<-lambda+delta
    c<-c+1
    gamma<-(c)^(-0.5)
  }
  lambda
}

Khn<-function(x,h){k=dnorm(x,mean=0,sd=h)
k}

Khn<-function(x,h){
  abb=(1-(x/h)^2)
  k=0.75*(abs(abb)+abb)/2
  k}


Khn_m<-function(x,y,h){
  abb=(1-(x/h)^2*(y/h)^2)
  k=0.75*(abs(abb)+abb)/2
  k}


ee_main<-function(alpha,h){
  n=length(alpha)
  fh=rep(0,n)
  for(i in 1:n)
  {fh[i]=mean(Khn(alpha-alpha[i],h))}
  
  fhh=rep(fh, 100)
  Sden<-matrix(fhh, nrow=n,byrow=FALSE)
  
  ad=kronecker(alpha, alpha, FUN="-")
  Snum<-matrix(Khn(ad,h)/n, nrow=n, byrow=TRUE)
  SS<-Snum/Sden
  
  IS=diag(1,nrow=n)-SS
  xs=(IS)%*%(xm)
  ys=IS%*%y
  betah=solve(t(xs)%*%xs)%*%t(xs)%*%ys
}



###calculate the estimate for beta by incorporating the information borrowing agent
ee_main_nt<-function(beta)
{
  n<-length(y)
  ee<-rep()
  for(i in 1:n)
  {
    ee_i<-as.vector(IS[i,]%*%xm)*(IS[i,]%*%(y-xm%*%beta))*Prop_scores[i]
    #ee_i<-c(IS[i,]%*%xm)*(IS[i,]%*%((y-xm%*%beta)*Prop_scores))
    ee<-cbind(ee,ee_i)
  }
  
  
  return(ee%*%rep(1,n))
}

ee_main<-function(beta)
{
  n<-length(y)
  ee<-rep()
  for(i in 1:n)
  {
    ee_i<-xm[i,]*((y-inv.logit(xm%*%beta))[i])*Prop_scores[i]
    ee<-cbind(ee,ee_i)
  }
  return(ee%*%rep(1,n))
}


SSmat_vcf<-function(u_vec, U_vec, X_vec,h){
  SS<-matrix(rep(0, length(u_vec)*length(U_vec)), nrow=length(u_vec))
  for(k in 1:length(u_vec))
  {
    Taumat<-matrix(c(rep(1,length(X_vec)),X_vec,(U_vec-u_vec[k]),(U_vec-u_vec[k])*X_vec), byrow=FALSE, ncol=4)
    TauW<-diag(Khn(U_vec-u_vec[k],h))
    tt=c(1,X_vec[k],0,0)%*%ginv(t(Taumat)%*%TauW%*%Taumat)%*%(t(Taumat)%*%TauW)
    SS[k,]=tt
  }
  SS
}



###############################################################################
start_time<-Sys.time()

pls_fp_t<-NULL
set.seed(1234567)
#n=200
n=600

#q=20
#q=1000
#beta<-c(1,1,1, 0.5,0.5)
#beta<-c(1,1,1)
beta<-c(1,1,1,1,1)
iteration<-100

# Mu=rep(0,q)
# Sigma=matrix(rep(0.2, q^2), nrow=q)
# diag(Sigma)=1


beta_ib_pool<-NULL
beta_pls_pool<-NULL

# ib2_fp_t<-NULL
#
rat=c(0.75,1.5,5)
bb=c(1)
for(r in rat)
{
  N<-n*r
  for(s in bb)
  {
    beta_ib<-NULL
    beta_ib_las<-NULL
    beta_initial<-NULL
    beta_initial_las<-NULL
    # beta2_pls_las<-NULL
    # beta2_ib_las<-NULL
    a0_mse_ib<-NULL
    a1z_mse_ib<-NULL
    a0_a1z_mse_ib<-NULL
    a0_mse_ini<-NULL
    a1z_mse_ini<-NULL
    a0_a1z_mse_ini<-NULL
    ib_fp_a<-NULL
    ini_fp_a<-NULL
    pls_fp_a<-NULL
    # ib2_fp_a<-NULL
    
    bound.pool<-rep()
    boundpls.pool<-rep()
    
    Varalls<-matrix(rep(0,25), nrow=5)
    Varallsp<-matrix(rep(0,25), nrow=5)
    
    for(i in 1:iteration)
    {
      ##external data No Penalty.
      tau_external=runif(N,0,1)
      xm1<-rnorm(N,tau_external,1)###############################################!!!!!!!!!!!!!!!!!!!!!!
      xm2<-rbinom(N,1,exp(tau_external)/(1+exp(tau_external)))
      xm3<-rnorm(N,0,1)
      #xm_external<-cbind(xm1,xm2,xm3)
      xm_external<-cbind(1, tau_external,xm1,xm2,xm3)
      ymean=xm_external%*%beta
      y_external=rbinom(N, 1, exp(ymean)/(1+exp(ymean)))
      
      ##internal data

      #Adjust to different external. xm2 is Binom(N, 1/(1+exp(tau))) and xm3 to rnorm(N, 1,2)
      tau=runif(n,0,1)
      xm1<-rnorm(n,tau,1)#############################################!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      xm2<-rbinom(n,1,exp(tau)/(1+exp(tau)))
      xm3<-runif(n, -1,1)
      #xm3<-rnorm(n, 0,1)
      
      #xm<-cbind(xm1,xm2,xm3)
      xm<-cbind(1,tau, xm1,xm2,xm3)
      ymean_i=xm%*%beta
      y=rbinom(n, 1, exp(ymean_i)/(1+exp(ymean_i)))
      #hh=s*n^(-0.2)
      
      ymall<-c(y_external, y)
      xmall<-rbind(xm_external, xm)[,-2]


      #ib2
      fit<-glm(ymall~xmall, family="binomial")
      theta<-fit$coefficients
      #theta<-theta_o[c(3,4,5)]
      #theta1<-theta_o[c(1,3,4,5)]
      #cst=N/(n+N)
      
      #Use sandwich for external V
      fit_e<-glm(y_external~xm_external[,3:5], family="binomial")   
      V_tilde<-vcov(fit_e)
      theta_e<-fit_e$coefficients
      
      xmm_e<-cbind(1, xm_external[,3:5])
      h<-0
      ep<-0
      for(k in 1:N)
      {
        hti=xmm_e[k,]*(y_external[k]-inv.logit(xmm_e[k,]%*%(theta_e)))
        sti=hti%*%t(hti)
        h=h+sti
        prob=inv.logit(xmm_e[k,]%*%theta_e)[1]
        ep=ep+xmm_e[k,]%*%t(xmm_e[k,])*(prob*(1-prob))
      }
      
      #SSt_e=ginv(h/N)
      #EPH_e=t(xmm_e)%*%xmm_e/N####!!!!!!!!!!!!!!!!!!!!!!!!!
      EPH_e=ep/N
      V_tilde=ginv(EPH_e)%*%(h/N)%*%ginv(EPH_e)/N
      
      #Use sandwich for internal V
      fit_i<-glm(y~xm[,3:5], family="binomial")   
      V_i<-vcov(fit_i)
      theta_i<-fit_i$coefficients
      
      xmm_i<-cbind(1,  xm[,3:5])
      h<-0
      ep<-0
      for(k in 1:n)
      {
        hti=xmm_i[k,]*(y[k]-inv.logit(xmm_i[k,]%*%(theta_i)))
        sti=hti%*%t(hti)
        h=h+sti
        prob=inv.logit(xmm_i[k,]%*%theta_i)[1]
        ep=ep+xmm_i[k,]%*%t(xmm_i[k,])*(prob*(1-prob))
      }
      #SSt_i=ginv(h/n)
      EPH_i=ep/n#!!!!!!!!!!!!!!!!!!!!!!
      V_i=ginv(EPH_i)%*%(h/n)%*%ginv(EPH_i)/n
      
      
      #V_t=(N/n)*ginv(V_tilde)+ginv(V_i)
      #V_t=(1/n)*ginv(V_tilde)+1/n*ginv(V_i)
      #V_t=V_tt[-3,-3]
      
      theta=ginv(n*ginv(V_i)+N*ginv(V_tilde))%*%t(n*fit_i$coefficients%*%ginv(V_i)+N*fit_e$coefficients%*%ginv(V_tilde))
      #theta=fit$coefficients
      
      
      ###information borrowing learner
      #calculate information borrowing agent
      #tuozi into xm matrix, linearly estimate it.
      #External model's XM
      #remove the intercept if generation process has one.
      xm_ee<-cbind(xm[,3:5])
      #xm_ee<-xm[,c(1,2,3)]
      lambda<-lambda_find(theta,y,xm_ee)
      ZZ<-ee_ib(theta,y,xm_ee)
      Prop_scores<-apply(ZZ,2,function(xx){1/(1+t(matrix(lambda,ncol=1))%*%xx)/n})
      
      ib_m<-glm(y~xm[,-1], family="binomial", weight=Prop_scores)
      beta_ib<-ib_m$coefficients[3:5]
      
      pls_m<-glm(y~xm[,-1], family="binomial")
      beta_pls<-pls_m$coefficients[3:5]
      
      beta_ib_pool<-c(beta_ib_pool, beta_ib)
      beta_pls_pool<-c(beta_pls_pool, beta_pls)
      
      print(i)
    }
   
    print(s)
  }
  print(r)
}


apply(matrix(beta_ib_pool, ncol=12, byrow=TRUE)[,-c(1,5,9)], 2, var)
apply(matrix(beta_pls_pool, ncol=12, byrow=TRUE)[,-c(1,5,9)], 2, var)

apply(matrix(beta_ib_pool, ncol=12, byrow=TRUE)[,-c(1,5,9)], 2, sd)
apply(matrix(beta_pls_pool, ncol=12, byrow=TRUE)[,-c(1,5,9)], 2, sd)


apply(matrix(beta_ib_pool-1, ncol=12, byrow=TRUE)[,-c(1,5,9)], 2, mean)
apply(matrix(beta_pls_pool-1, ncol=12, byrow=TRUE)[,-c(1,5,9)], 2, mean)


end_time<-Sys.time()

end_time-start_time
