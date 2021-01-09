
#the matrix of  hat{G} in the main text.
G_mat=function(X,S,N,p)  
{
  n=0
  G=matrix(0,p,p)
  for (i in 1:(2*N))
  { for (j in 1:(2*N))
  { 
    if (S[i,j]!=0)
    {
      X_ij=matrix(X[i,]-X[j,],nrow=1)
      G=G+S[i,j]*t(X_ij)%*%X_ij
    }
    n=n+1
  }
  }
  G=G/(n-2*N)
}

# compute the clasification error for the NDR based K-means algorithm.  
error_min=function(Est_label,True_label,N) 
{ 
  #library(e1071)
  Unique=unique(Est_label)
  error=rep(0,factorial(length(Unique))) 
  V=matrix(0,length(Est_label),length(Unique))
  for (i in 1:length(Est_label))
  {Id_i=which(Unique==Est_label[i])
  V[i,Id_i]=1
  }
  value_perm=t(permutations(length(Unique)))
  
  Label_perm=V%*%value_perm
  
  for (i in 1:factorial(length(Unique)))
  {
    error[i]=length(which(True_label!=Label_perm[,i]))/(num_block*N)
  } 
  error=min(error)
  
}


############################### main function

main=function(X,W,p,N,True_label,num_block,Iteration=FALSE)
{
  ##################
 
      S=1-W # adjacency matrix W
      G=G_mat(X,S,N,p)
      A=cov(X)
    
      ################################## NDR(N-KM)
      starttime=Sys.time()

      X_proj=vector();
      error=vector()
      if(Iteration==FALSE){
        A_eig=eigen(A)
        A_sqrt=A_eig$vectors%*%diag(1/sqrt(A_eig$values))%*%t(A_eig$vectors)
        
        E_lda=eigen(A_sqrt%*%G%*%A_sqrt)$values
        rank_r=c()
        
        if(num_block==2){
          r=1
        }else{
          for (i in 1:3) {
            rank_r[i]=(E_lda[i]-E_lda[i+1])/(E_lda[i]+E_lda[i+1])
          }
          r=which.max(rank_r) 
        }
        
        Id_eigv=c(1:r)
        
        V_lda=A_sqrt%*%(eigen(A_sqrt%*%G%*%A_sqrt)$vectors[,Id_eigv])
        
        X_proj=X%*%V_lda
        
        kc=kmeans(X_proj,num_block)
        SumSquare=kc$tot.withinss/kc$totss
        
        error=error_min(kc$cluster,True_label,N) 
     
      }else{
    
        n_iter=0 
        SumSquare=1000; SumSquare_old=SumSquare+1;
        
        while(SumSquare <SumSquare_old-0.05) # Iteration
        {SumSquare_old=SumSquare
        X_proj_old=X_proj
        error_old=error
        
        A_eig=eigen(A)
        A_sqrt=A_eig$vectors%*%diag(1/sqrt(A_eig$values))%*%t(A_eig$vectors)
        
        E_lda=eigen(A_sqrt%*%G%*%A_sqrt)$values
        rank_r=c()
        
        
        if(num_block==2){
          r=1
        }else{
          for (i in 1:3) {
            rank_r[i]=(E_lda[i]-E_lda[i+1])/(E_lda[i]+E_lda[i+1])
          }
          r=which.max(rank_r) 
        }
        
        Id_eigv=c(1:r)
        V_lda=A_sqrt%*%(eigen(A_sqrt%*%G%*%A_sqrt)$vectors[,Id_eigv])
        
        ####
        
        X_proj=X%*%V_lda
        
        kc=kmeans(X_proj,num_block)
        SumSquare=kc$tot.withinss/kc$totss
        
        error=error_min(kc$cluster,True_label,N)
        
        if (SumSquare_old<SumSquare) # if the results get worse than the previous round, we use the previous one. 
        {error=error_old
        X_proj=X_proj_old
        }
        
        ####### estimating the covariance matrix 
        I=c();COV=0
        for (k in 1:num_block) {
          I=which(kc$cluster==k)
          if(length(I)==1){
            COV=COV+var(X[I,]) 
          }
          else{
            COV=COV+cov(X[I,]) 
          }
          
        }
        hatA=COV/num_block
        A=hatA #
        n_iter=n_iter+1
        }# corresponding to while

      }
      
      endtime=Sys.time()
      timeused=endtime-starttime
      cat("error",error,"\n")
      cat("timeused",timeused,"\n")

  list("Error"=error,"timeused"=timeused)
  
}


