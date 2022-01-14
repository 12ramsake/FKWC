
#You need the code from the other authors for this one

# install.packages("EMMIXskew")
# library(EMMIXskew)
library(sn)
library(fda)
library(roahd)
library(mvtnorm)
library(fdasrvf)
library(doParallel)
library(MFHD)
library(mrfDepth)
library(doRNG)

library(ECF)
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

###generate a set of functional data with specific kernel, 
##and N is sample size
#generate data from spec. covariance kernel

###generate a set of functional data with specific kernel, 
##and N is sample size
#generate data from spec. covariance kernel
gfd_k_3<-function(N1=10,N2=10,N3=10,
                  res=1e2, 
                  Cov1,
                  Cov2,
                  Cov3,
                  dist="t",
                  del=0.9){
  
  
  grid = seq( 0, 1, length.out = res )
  
  
  if(dist=="N"){
    Data1 = generate_gauss_fdata( N1,  rep(0,res),Cov = Cov1 )
    Data2 = generate_gauss_fdata( N2, rep(0,res), Cov = Cov2 )
    Data3 = generate_gauss_fdata( N3, rep(0,res), Cov = Cov3 )
  }
  else if(dist=="t"){
    Data1 = rmvt( N1, sigma =  Cov1/3,df=3 )
    Data2 = rmvt( N2, sigma =  Cov2/3,df=3 )
    Data3 = rmvt( N3, sigma =  Cov3/3,df=3 )
  }
  else{
    # Data1 = rdmsn(N1, res, rep(0,res),Cov1, del=rep(del,res))
    # Data2 = rdmsn(N2, res, rep(0,res),Cov2, del=rep(del,res))
    # Data3 = rdmsn(N3, res, rep(0,res),Cov3, del=rep(del,res))
    cp1=list(mean=rep(0,nrow(Cov1)), var.cov=Cov1, gamma1=rep(del,nrow(Cov1))/nrow(Cov1))
    cp2=list(mean=rep(0,nrow(Cov2)), var.cov=Cov2, gamma1=rep(del,nrow(Cov2))/nrow(Cov2))
    cp3=list(mean=rep(0,nrow(Cov3)), var.cov=Cov3, gamma1=rep(del,nrow(Cov3))/nrow(Cov3))
    dp1=cp2dp(cp1, "SN")
    dp2=cp2dp(cp2, "SN")
    dp3=cp2dp(cp3, "SN")
    
    Data1 = rmsn(N1, dp=dp1)
    Data2 = rmsn(N2, dp=dp2)
    Data3 = rmsn(N3, dp=dp3)
    
  }
  
  return(list(argvals =  grid, mdata=rbind(Data1,Data2,Data3)))
  
}


oneRun<-function(N1,N2,N3,c1,c2,c3,res=100,
                 grid=seq( 0, 1, length.out = res ),
                 dist="t",
                 del=0.9){
  
  dat=gfd_k_3(N1,N2,N3,Cov1=c1,Cov2=c2,Cov3=c3,res=res,dist=dist,del=del)
  #boente
   BNT=testoper3P(grid,dat$mdata[1:N1,],dat$mdata[(N1+1):(N2+N1),],dat$mdata[(N2+N1+1):(N2+N1+N3),],10,5000)$pvalor
  #guo paper
  # require("ECF");
  sample_n=c(N1,N2,N3)
  
  datamx <- cbind(t(dat$mdata[1:N1,]),t(dat$mdata[(N1+1):(N2+N1),]),t(dat$mdata[(N2+N1+1):(N2+N1+N3),]));
  # listver=list(t(dat$mdata[1:N1,]),t(dat$mdata[(N1+1):(N2+N1),]));
  PV_L2nv <- ECFL2Rcpp(datamx,sample_n,0,centered=1);#L2nv
  PV_L2br <- ECFL2Rcpp(datamx,sample_n,1,centered=1);#L2br
  PV_Fnv <- ECFFnaiveRcpp(datamx,sample_n,centered=1);#GPFnv
  PV_rps <- rcpparma_ECFRPall(datamx,sample_n,1000,1);#L2rp,Tmax,GPFrp,Fmax
  #P-values of Boente L2nv,L2br,L2rp,Tmax,GPFnv,GPFrp,Fmax
  # KSCovsup(listver)
  pvals=c(BNT,PV_L2nv[2], PV_L2br[2], PV_rps[1], PV_rps[2], PV_Fnv[2], PV_rps[3], PV_rps[4])
  names(pvals)=c("Boen","L2nv","L2br","L2rp","Tmax","GPFnv","GPFrp","Fmax")
  # pvals=c(PV_L2nv[2], PV_L2br[2], PV_rps[1], PV_rps[2], PV_Fnv[2], PV_rps[3], PV_rps[4])
  # names(pvals)=c("L2nv","L2br","L2rp","Tmax","GPFnv","GPFrp","Fmax")
  
  return(pvals)
  
  
}
runMVSim<-function(N1,N2,N3,c1,c2,c3,grid,num_runs,fileName,dist="t",del=0.9){
  
  
  
  no_cores<-100
  
  cl <- makeCluster(no_cores,type="FORK")
  # cl <- makeCluster(no_cores)
  registerDoParallel(cl)
  registerDoRNG(seed = 440)
  # set.seed(440)
  p_values=try({foreach(i=1:num_runs) %dopar% {oneRun(N1,N2,N3,c1,c2,c3,res=100,del=del,dist=dist,grid=grid)}})
  
  errorsp=inherits(p_values, "try-error")
  if(errorsp){
    print("there was an error!")
    
  }
  
  stopCluster(cl)
  registerDoSEQ()
  closeAllConnections()
  dirr="~/OTHER_TEST_RESULTS/"
  save(p_values,file=paste(dirr,fileName,sep=""))
  # dirr="Functional Data Covariance Files/Two_Sample_Depths_Cov_Kernels/"
  # save(p_values,file=paste(dirr,fileName,sep=""))
}


##
N1=N2=N3=50
N1=N2=N3=100
num_runs=200
grid=seq( 0, 1, length.out = 100 )

for(N1 in c(50,100)){
  N2=N1
  print("sim 1")
  #changes in beta
  count=0
  k1=function(t1,t2){k_exp(t1,t2,.5,.05)}
  c1=getCov(grid,k1)
  c3=getCov(grid,k1)
  tmp=seq(0.5,1,length.out = 15)
  for(dist in c("N","t","sn")){
  # for(dist in c("t")){
    for(betai in 2:15){
      beta=tmp[betai]
      fileName=paste0("dist_",dist,"_other_k_exp_alpha1_p05_p05_p05_beta_p5_p",betai,"_p5_N_",N1,"_num_runs_200.Rda",sep="")
      k2=function(t1,t2){k_exp(t1,t2,beta,.05)}
      c2=getCov(grid,k2)
      runMVSim(N1,N2,N3,c1,c2,c3,grid,num_runs,fileName,dist=dist)
      count=count+1
      print(count/(42))
    }
  }
  
  

  print("sim 2")
  #changes in alpha
  count=0
  k1=function(t1,t2){k_exp(t1,t2,.5,.05)}
  c1=getCov(grid,k1)
  c3=getCov(grid,k1)
  tmp=seq(0.05,0.1,length.out=15)
  for(dist in c("N","t","sn")){
  # for(dist in c("t")){
    for(alphai in 2:15){
      alpha=tmp[alphai]
      fileName=paste0("dist_",dist,"_other_k_exp_alpha1_p05_",alphai,"_p05_beta_p5_p5_p5_N_",N1,"_num_runs_200.Rda",sep="")
      k2=function(t1,t2){k_exp(t1,t2,0.5,alpha)}
      c2=getCov(grid,k2)
      runMVSim(N1,N2,N3,c1,c2,c3,grid,num_runs,fileName,dist=dist)
      count=count+1
      print(count/(42))
    }
  }
}




for(N1 in c(50,100)){
  N2=N1
  print("sim 1")
  #changes in beta
  count=0
  k1=function(t1,t2){k_exp(t1,t2,.5,.05)}
  c1=getCov(grid,k1)
  c3=getCov(grid,k1)
  tmp=seq(0.5,1,length.out = 15)
  for(dist in c("N","t","sn")){
    # for(dist in c("t")){
    for(betai in 1){
      beta=tmp[betai]
      fileName=paste0("dist_",dist,"_other_k_exp_alpha1_p05_p05_p05_beta_p5_p",betai,"_p5_N_",N1,"_num_runs_200.Rda",sep="")
      k2=function(t1,t2){k_exp(t1,t2,beta,.05)}
      c2=getCov(grid,k2)
      runMVSim(N1,N2,N3,c1,c2,c3,grid,num_runs,fileName,dist=dist)
      count=count+1
      print(count/(42))
    }
  }
  
  

  print("sim 2")
  #changes in alpha
  count=0
  k1=function(t1,t2){k_exp(t1,t2,.5,.05)}
  c1=getCov(grid,k1)
  c3=getCov(grid,k1)
  tmp=seq(0.05,0.1,length.out=15)
  for(dist in c("N","t","sn")){
  # for(dist in c("t")){
    for(alphai in 1){
      alpha=tmp[alphai]
      fileName=paste0("dist_",dist,"_other_k_exp_alpha1_p05_",alphai,"_p05_beta_p5_p5_p5_N_",N1,"_num_runs_200.Rda",sep="")
      k2=function(t1,t2){k_exp(t1,t2,0.5,alpha)}
      c2=getCov(grid,k2)
      runMVSim(N1,N2,N3,c1,c2,c3,grid,num_runs,fileName,dist=dist)
      count=count+1
      print(count/(42))
    }
  }
}
