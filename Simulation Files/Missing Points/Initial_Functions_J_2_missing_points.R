
# %%first compare depth functions
# install.packages("EMMIXskew")
library(sn)
library(fda)
library(roahd)
library(mvtnorm)
library(fdasrvf)
library(doParallel)
library(MFHD)
library(mrfDepth)
library(doRNG)
####kernels and simulating data

#beta will be the scale parameter for any kernel
#different kernels
k_linear<-function(t1,t2,beta,c,inte){inte+beta*(t1-c)*(t2-c)}
k_exp<-function(t1,t2,beta,alpha,dummy){beta*exp(-(t1-t2)^2/(2*alpha^2))}
# k_exp2<-function(t1,t2,beta,alpha){alpha*exp(-beta*abs(t1-t2))}
k_periodic<-function(t1,t2,beta,alpha,p){beta*exp(-2/alpha^2*(sin(pi*(t1-t2)/p)^2))}


##get covariance over a grid, makes pd if not bc rounding
getCov<-function(grid,k){
  return(corpcor::make.positive.definite(outer(grid,grid,k)))
}


gfd<-function(N1=10,N2=10,
              res=1e2, 
              Cov1,
              Cov2,
              dist="t",
              del=0.9){
  
  
  grid = seq( 0, 1, length.out = res )
  
  
  if(dist=="N"){
    Data1 = generate_gauss_fdata( N1,  rep(0,res),Cov = Cov1 )
    Data2 = generate_gauss_fdata( N2, rep(0,res), Cov = Cov2 )
  }
  else if(dist=="t"){
    Data1 = rmvt( N1,  sigma =  Cov1/3,df=3 )
    Data2 = rmvt( N2,  sigma =  Cov2/3,df=3 )
  }
  else{
    cp1=list(mean=rep(0,nrow(Cov1)), var.cov=Cov1, gamma1=rep(del,nrow(Cov1))/nrow(Cov1))
    cp2=list(mean=rep(0,nrow(Cov2)), var.cov=Cov2, gamma1=rep(del,nrow(Cov2))/nrow(Cov2))
    dp1=cp2dp(cp1, "SN")
    dp2=cp2dp(cp2, "SN")
    
    Data1 = rmsn(N1, dp=dp1)
    Data2 = rmsn(N1, dp=dp2)
  }
  
  return(list(argvals =  grid, mdata=rbind(Data1,Data2)))
  
}


modify_curve=function(curve,amount_to_delete=0.2){
  
  curve[sample(2:(length(curve)-1),floor(amount_to_delete*length(curve)))]=rep(NA,floor(amount_to_delete*length(curve)))

  return(curve)
}

runMVSim<-function(N1,N2,c1,c2,grid,num_runs,fileName,dist="t",del=0.9){
  

    no_cores<-detectCores()-1
    no_cores<-50
    
    cl <- makeCluster(no_cores,type="FORK")   
    registerDoParallel(cl) 
    registerDoRNG(seed = 440)
    
    depth_values=try({foreach(i=1:num_runs) %dopar% {simulate_mv_depths(N1,N2,c1,c2,res=100,del=del,dist=dist,grid=grid)}})
    
    errorsp=inherits(depth_values, "try-error")
    if(errorsp){
      print("there was an error!")
      
    }

    stopCluster(cl)
    registerDoSEQ()
    closeAllConnections()
    
    dirr="/u/k3ramsay/ResearchDocuments/output/Functional Data Covariance Files/"
    save(depth_values,file=paste(dirr,fileName,sep=""))
  # }
}


num_runs=200
grid=seq( 0, 1, length.out = 100 )

