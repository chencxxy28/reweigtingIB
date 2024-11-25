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
  #add a stop after 100 iteration.
  repeat{
    ssss=ssss+1
    rl<-R1der(lambda,ZZ)
    rll<-R2der(lambda,ZZ) 
    Delta<--ginv(rll)%*%rl
    if(mean(abs(Delta))<tol |ssss>100)
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

##########
ee_main_all<-function(theta,y,xm)
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

ee_main_ib<-function(beta_a)
{
  n<-length(y_obs_i)
  ee<-NULL
  for(i in 1:n)
  {
    beta<-beta_a[1:2]
    eta<-beta_a[3:7]
    Am<-cbind(1, Assignment_i)
    Xm<-cbind(1, xm_internal)
    
    #ee_i<-xm[i,]*((y-inv.logit(xm%*%beta))[i])*Prop_scores[i]*inv_weight_i[i]
    ee1_i<-Am[i,]*((y_obs_i-inv.logit(Am%*%beta))[i])*Prop_scores[i]*inv_weight_i[i]
    ee2_i<-Xm[i,]*((Assignment_i-inv.logit(Xm%*%eta))[i])*Prop_scores[i]
    ee_i<-c(ee1_i,ee2_i)
    ee<-cbind(ee,ee_i)
  }
  return(ee%*%rep(1,n))
}


ee_main_iptw<-function(beta_a)
{
  n<-length(y_obs_i)
  ee<-NULL
  for(i in 1:n)
  {
    beta<-beta_a[1:2]
    eta<-beta_a[3:7]
    Am<-cbind(1, Assignment_i)
    Xm<-cbind(1, xm_internal)
    
    #ee_i<-xm[i,]*((y-inv.logit(xm%*%beta))[i])*Prop_scores[i]*inv_weight_i[i]
    ee1_i<-Am[i,]*((y_obs_i-inv.logit(Am%*%beta))[i])*inv_weight_i[i]
    ee2_i<-Xm[i,]*((Assignment_i-inv.logit(Xm%*%eta))[i])
    ee_i<-c(ee1_i,ee2_i)
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
set.seed(123456)
#n=200
n=600

#q=20
#q=1000
#beta<-c(1,1,1, 0.5,0.5)
beta_1<-c(-0.5,0.5,-0.5,0.5,-0.5)
beta_0<-c(0.5,-0.5, 0.5,-0.5,0.5)
beta_A<-c(0.5,0.5,0.5,0.5,0.5)
iteration<-1000

# Mu=rep(0,q)
# Sigma=matrix(rep(0.2, q^2), nrow=q)
# diag(Sigma)=1


beta_ib_pool<-NULL
beta_causal_pool<-NULL
tor<-NULL
# ib2_fp_t<-NULL
#
rat=c(0.75, 1.5,5)
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
      ##external data
      tau_external=runif(N,0,1)
      xm1<-rnorm(N,tau_external,1)###############################################!!!!!!!!!!!!!!!!!!!!!!
      xm2<-rbinom(N,1,exp(tau_external)/(1+exp(tau_external)))
      xm3<-rnorm(N,0,1)
      xm_external<-cbind(tau_external,xm1,xm2,xm3)####check
      ymean_e1=cbind(1,xm_external)%*%beta_1
      ymean_e0=cbind(1,xm_external)%*%beta_0
      y_external_1=rbinom(N, 1, exp(ymean_e1)/(1+exp(ymean_e1)))
      y_external_0=rbinom(N, 1, exp(ymean_e0)/(1+exp(ymean_e0)))
      
      A_e=cbind(1,xm_external)%*%beta_A
      Assignment=rbinom(N, 1, inv.logit(A_e))
      y_obs_e=ifelse(Assignment==1, y_external_1, y_external_0)
      
      # odds_1=mean(y_external_1)/(1-mean(y_external_1))
      # odds_0=mean(y_external_0)/(1-mean(y_external_0))
      # odds_1/odds_0
      ##Estimation.
      am1<-glm(Assignment~xm_external, family = "binomial")
      raw_weight<-am1$fitted.values
      inv_weight<-1/ifelse(Assignment==1, raw_weight, 1-raw_weight)
      fm1<-glm(factor(y_obs_e)~Assignment, family = "binomial", weights = inv_weight)
      # 
      #       tor<-NULL
      #       eor<-NULL
      #       for(s in 1:150)
      #       {
      #         ##internal data
      tau_internal=runif(n,0,1)
      xm1<-rnorm(n,tau_internal,1)###############################################!!!!!!!!!!!!!!!!!!!!!!
      xm2<-rbinom(n,1,exp(tau_internal)/(1+exp(tau_internal)))
      xm3<-rnorm(n,0,1)
      #xm3<-runif(n,-1,1)
      xm_internal<-cbind(tau_internal,xm1,xm2,xm3)
      ymean_i1=cbind(1,xm_internal)%*%beta_1
      ymean_i0=cbind(1,xm_internal)%*%beta_0
      y_internal_1=rbinom(n, 1, exp(ymean_i1)/(1+exp(ymean_i1)))
      y_internal_0=rbinom(n, 1, exp(ymean_i0)/(1+exp(ymean_i0)))
      
      A_i=cbind(1,xm_internal)%*%beta_A
      Assignment_i=rbinom(n, 1, inv.logit(A_i))
      y_obs_i=ifelse(Assignment_i==1, y_internal_1, y_internal_0)
      
      am_int<-glm(Assignment_i~xm_internal, family="binomial")
      raw_weight_i<-am_int$fitted.values
      inv_weight_i<-1/ifelse(Assignment_i==1, raw_weight_i, 1-raw_weight_i)
      fm_log<-glm(factor(y_obs_i)~Assignment_i, family = "binomial", weights = inv_weight_i)
      
      #True odds
      #   odds_1=mean(y_internal_1)/(1-mean(y_internal_1))
      #   odds_0=mean(y_internal_0)/(1-mean(y_internal_0))
      #   true_odds=odds_1/odds_0
      #   tor=c(tor, true_odds)
      #   eor=c(eor, exp(fm_log$coefficients[2]))
      #   
      # }
      # 
      # mean(log(tor))
      # mean(log(eor))
      # 
      
      #Now calculate propensity score
      ##First pool means using meta
      
      #Use sandwich for external V
      fit_e<-glm(y_obs_e~Assignment+xm_external, family="binomial")
      #fit_e<-glm(y_obs_e~Assignment, family="binomial")
      V_tilde<-vcov(fit_e)
      theta_e<-fit_e$coefficients
      
      xmm_e<-cbind(1, Assignment,xm_external)
      #xmm_e<-cbind(1, Assignment)
      h<-0
      ep<-0
      for(k in 1:N)
      {
        hti=xmm_e[k,]*(y_obs_e[k]-inv.logit(xmm_e[k,]%*%(theta_e)))
        sti=hti%*%t(hti)
        h=h+sti
        prob=inv.logit(xmm_e[k,]%*%theta_e)[1]
        ep=ep+xmm_e[k,]%*%t(xmm_e[k,])*(prob*(1-prob))
      }
      
      #SSt_e=ginv(h/N)
      EPH_e=t(xmm_e)%*%xmm_e/N
      EPH_e=ep/N
      V_tilde=ginv(EPH_e)%*%(h/N)%*%ginv(EPH_e)/N
      
      #Use sandwich for internal V
      fit_i<-glm(y_obs_i~Assignment_i+xm_internal, family="binomial")
      #fit_i<-glm(y_obs_i~Assignment_i, family="binomial")
      V_i<-vcov(fit_i)
      theta_i<-fit_i$coefficients
      
      xmm_i<-cbind(1, Assignment_i , xm_internal)
      #xmm_i<-cbind(1, Assignment_i )
      h<-0
      ep<-0
      for(k in 1:n)
      {
        hti=xmm_i[k,]*(y_obs_i[k]-inv.logit(xmm_i[k,]%*%(theta_i)))
        sti=hti%*%t(hti)
        h=h+sti
        prob=inv.logit(xmm_i[k,]%*%theta_i)[1]
        ep=ep+xmm_i[k,]%*%t(xmm_i[k,])*(prob*(1-prob))
      }
      #SSt_i=ginv(h/n)
      #EPH_i=t(xmm_i)%*%xmm_i/n
      EPH_i=ep/n
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
      xm_ee<-cbind(Assignment_i, xm_internal)
      #xm_ee<-cbind(Assignment_i)
      lambda<-lambda_find(theta,y_obs_i,xm_ee)
      ZZ<-ee_ib(theta,y_obs_i,xm_ee)
      Prop_scores<-apply(ZZ,2,function(xx){1/(1+t(matrix(lambda,ncol=1))%*%xx)/n})
      
      ib_m<-glm(y_obs_i~Assignment_i, family="binomial", weight=Prop_scores*inv_weight_i)
      beta_ib_ini_a<-ib_m$coefficients
      beta_ib_ini_b<-am_int$coefficients
      beta_ib_ini=c(beta_ib_ini_a, beta_ib_ini_b)
      beta_ib<-multiroot(f = ee_main_ib, start = as.vector(beta_ib_ini))$root
      beta_iptw<-multiroot(f = ee_main_iptw, start = as.vector(beta_ib_ini))$root
      beta_ib
      beta_iptw
      
      
      
      
      
      # causal_m<-glm(y_obs_i~Assignment_i, family="binomial", weight=inv_weight_i)
      # beta_causal<-causal_m$coefficients
      
      beta_ib_pool<-rbind(beta_ib_pool, beta_ib)
      beta_causal_pool<-rbind(beta_causal_pool, beta_iptw)
      
      print(i)
      
      odds_1=mean(y_internal_1)/(1-mean(y_internal_1))
      odds_0=mean(y_internal_0)/(1-mean(y_internal_0))
      true_odds=odds_1/odds_0
      tor=c(tor, log(true_odds))
      
    }
    
    print(s)
  }
  print(r)
}


# var_ib<-apply(matrix(beta_ib_pool, ncol=6, byrow=TRUE)[,-c(1,3,5)], 2, var)
# var_causal<-apply(matrix(beta_causal_pool, ncol=6, byrow=TRUE)[,-c(1,3,5)], 2, var)
# 
# apply(matrix(beta_ib_pool, ncol=6, byrow=TRUE)[,-c(1,3,5)], 2, sd)
# apply(matrix(beta_causal_pool, ncol=6, byrow=TRUE)[,-c(1,3,5)], 2, sd)

true_odds<-apply(matrix(tor, nrow=3, byrow=TRUE), 1, mean)

# bias_ib<-apply(matrix(beta_ib_pool, ncol=6, byrow=TRUE)[,-c(1,3,5)], 2, mean)-true_odds
# bias_causal<-apply(matrix(beta_causal_pool, ncol=6, byrow=TRUE)[,-c(1,3,5)], 2, mean)-true_odds
# 
# (bias_causal^2+var_causal)/(bias_ib^2+var_ib)
# 
# bias_ib
# bias_causal


bias_ib<-apply(matrix(beta_ib_pool[,2], nrow=3, byrow=TRUE), 1, mean)-true_odds
bias_causal<-apply(matrix(beta_causal_pool[,2], nrow=3, byrow=TRUE), 1, mean)-true_odds

apply(matrix(beta_ib_pool[,2], nrow=3, byrow=TRUE), 1, sd)
apply(matrix(beta_causal_pool[,2], nrow=3, byrow=TRUE), 1, sd)

var_ib<-apply(matrix(beta_ib_pool[,2], nrow=3, byrow=TRUE), 1, var)
var_causal<-apply(matrix(beta_causal_pool[,2], nrow=3, byrow=TRUE), 1, var)

re=(var_causal+bias_causal^2)/(var_ib+bias_ib^2)

re


end_time<-Sys.time()

end_time-start_time


##TEST

var_ib<-apply(matrix(beta_ib_pool, ncol=4, byrow=TRUE)[,-c(1,3)], 2, var)
var_causal<-apply(matrix(beta_causal_pool, ncol=4, byrow=TRUE)[,-c(1,3)], 2, var)

apply(matrix(beta_ib_pool, ncol=4, byrow=TRUE)[,-c(1,3)], 2, sd)
apply(matrix(beta_causal_pool, ncol=4, byrow=TRUE)[,-c(1,3)], 2, sd)

true_odds<-apply(matrix(tor, nrow=2, byrow=TRUE), 1, mean)

bias_ib<-apply(matrix(beta_ib_pool, ncol=4, byrow=TRUE)[,-c(1,3)], 2, mean)-true_odds
bias_causal<-apply(matrix(beta_causal_pool, ncol=4, byrow=TRUE)[,-c(1,3)], 2, mean)-true_odds

(bias_causal^2+var_causal)/(bias_ib^2+var_ib)


