#############################################################
############################################################


## This is the code used for example in section 2. 

rm(list=ls())

Sigma_Adj=function(p,N,struc,a1,a2,b,c,s) # generate covariance matrix and adjcent matrix. 
{B=matrix(0,p,p)
for(i in 1:p)
{
  for(j in i:p)
  {
    if (abs(i-j)<p)
    {B[i,j]=0.7^(abs(i-j)/3)
    B[j,i]=B[i,j]}
  }
}
set.seed(s)
Paa=matrix(as.numeric(runif(N^2,0,1)<a1),N,N)
Pbb=matrix(as.numeric(runif(N^2,0,1)<a2),N,N)
Pcc=matrix(as.numeric(runif(N^2,0,1)<0),N,N)
Pab=matrix(as.numeric(runif(N^2,0,1)<b),N,N)


if(struc==1)
{
  W=rbind(cbind(Paa,Pab,c*Paa,c*Pab),cbind(Pab,Pbb,c*Pab,c*Pbb),cbind(c*Paa,c*Pab,Paa,Pab),cbind(c*Pab,c*Pbb,Pab,Pbb))
}
if (struc==2)
{
  W=rbind(cbind(Paa,Paa,Pab,Pab),cbind(Paa,Paa,Pab,Pab),cbind(Pab,Pab,Pbb,Pbb),cbind(Pab,Pab,Pbb,Pbb))
  
}
if (struc==3)
{
  W=rbind(cbind(Paa,Pab,Paa,Pab),cbind(Pab,Pbb,Pab,Pbb),cbind(Paa,Pab,c*Paa,c*Pab),cbind(Pab,Pbb,c*Pab,c*Pbb))
  
}
diag(W)=0


list(Sigma=B,Adj_mat=W)

}


library(MASS)
library(igraph)
library(e1071)

N=50; p=5; 
num_sim=100 # given an adjcent matrix repeat X for num_sim times. 
num_block=2
struc=2 #1:(2c); 2:(2a); 3:(2b)

mean_shift=1.6
rho=0.1 #rho=1,0.05
a1=0.5 
a2=0.1
c=0.1 ## 0.2,0.05 
b=seq(0,0.7,0.1)

SEED=seq(from=100,to=30000,by=33)  # random seed 
N_simu=10   # generate random adjcent matrix for N_simu times. 


timeused=matrix(nrow=N_simu,ncol=length(b))
Error<-matrix(nrow=N_simu,ncol=length(b))
################# 
for (xi in c(1:N_simu))
{
  Seed=SEED[xi]
  for(r in 1:length(b)){
    S_A=Sigma_Adj(p,N/2,struc,rho*a1,rho*a2,rho*b[r],c,Seed)
    B=S_A$Sigma
    W=S_A$Adj_mat
    
    #true label
    if(num_block==2){
      True_label=c(rep(1,N),rep(2,N));
    }else{
      True_label=c(rep(1,N/2),rep(2,N/2),rep(3,N/2),rep(4,N/2));
    }
    
    ############### covariates
    mu1=c(mean_shift,rep(0,p-1)); mu2=-mu1;
    mu3=c(0,mean_shift,rep(0,p-2))
    mu4=-mu3
    
    if(num_block==2){
      X1 <- mvrnorm(n=N, mu=mu1, Sigma=B)
      X2 <- mvrnorm(n=N, mu=mu2, Sigma=B)
      X=rbind(X1,X2)
    }else{
      X1 <- mvrnorm(n=N/2, mu=mu1, Sigma=B)
      X2 <- mvrnorm(n=N/2, mu=mu2, Sigma=B)
      X3 <- mvrnorm(n=N/2, mu=mu3, Sigma=B)
      X4 <- mvrnorm(n=N/2, mu=mu4, Sigma=B)
      X=rbind(X1,X2,X3,X4)
    }
    
    result=main(X,W,p,N,True_label,num_block,Iteration=FALSE)
    Error[xi,r]=result[[1]]
    timeused[xi,r]=result[[2]]
  }
  
  
  
} 

result=list('Ave_err_NDR-KM'=apply(Error,2,mean),'Ave_timeused1'=apply(timeused,2,mean))

print(result)


