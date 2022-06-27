

# %%Basic -- Gaussian 
# %%heavy tails -- t
# %%skewed distributions -- sn

# install.packages("EMMIXskew")

setwd("/u/k3ramsay/ResearchDocuments/output/Functional Data Covariance Files/")
# install.packages(c("Guo/ECF_1.2.tar.gz","fChange_0.2.1.tar.gz"),type="source",repos=NULL)
require(ECF)
library(EMMIXskew)
library(fda)
library(roahd)
library(mvtnorm)
library(fdasrvf)
library(doParallel)
library(MFHD)
library(mrfDepth)
library(doRNG)
library(fChange)




gen_eigen_model<-function(N1,N2,nbasis=55,b1=rep(1,nbasis),b2=rep(1,nbasis)){
  
  
  dat1=fun_IID(N1, nbasis,Sigma=b1)
  dat2=fun_IID(N2, nbasis,Sigma=b2)
  dat=dat1
  dat$coefs=cbind(dat1$coefs,dat2$coefs)
  dat$fdnames$reps=c(dat$fdnames$reps,paste(rep("reps ",N2),as.character(((N1+1):(N1+N2)),collapse ="")))
  
  return(dat)
  
}


oneRun<-function(both,out,contamination_level,N1,N2,b1,b2,grid=seq( 0, 1, length.out = 100 )){
  
  
  
  fdat=gen_eigen_model(N1,N2,length(b1),b1,b2)
  dat=t(eval.fd(grid,fdat))
  
  dat1=list()
  dat1$argvals=grid
  dat1$mdata=dat
  
  # fderivs=deriv.fd(fdat,1)
  #
  dat=generate_outliers(dat1,N1,N2,contamination_level,outlier_type=out,both_samples=both)
  # matplot(t(dat$mdata),type='l')
  
  fdat=Data2fd(dat$argvals,t(dat$mdata))
  fderivs=deriv.fd(fdat,1)
  
  dat$arg=grid
  
  
  # dat=list("arg"=grid,"mdata"=t(eval.fd(grid,fdat)))
  #boente
  # source('~/research/PhD Thesis/Functional Data Covariance Files/Boente_MOd.R')
  BNT=testoper(grid,dat$mdata[1:N1,],dat$mdata[(N1+1):(N2+N1),],10,5000)$pvalor
  #guo paper
  # require("ECF");
  sample_n=c(N1,N2)
  
  datamx <- cbind(t(dat$mdata[1:N1,]),t(dat$mdata[(N1+1):(N2+N1),]));
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

runMVSim<-function(both,out,contamination_level,N1,N2,b1,b2,grid,num_runs,fileName){
  
  
  no_cores<-detectCores()-1
  no_cores<-50
  
  cl <- makeCluster(no_cores,type="FORK")   
  registerDoParallel(cl) 
  registerDoRNG(seed = 440)
  # set.seed(440)
  p_values=try({foreach(i=1:num_runs) %dopar% {oneRun(both,out,contamination_level,N1,N2,b1,b2,grid)}})
  
  # p_values=try({replicate(num_runs, oneRun(both,out,contamination_level,N1,N2,b1,b2,grid))})
  
  errorsp=inherits(p_values, "try-error")
  if(errorsp){
    print("there was an error!")
    
  }
  
  stopCluster(cl)
  registerDoSEQ()
  closeAllConnections()
  
  dirr="/u/k3ramsay/ResearchDocuments/output/Functional Data Covariance Files/EIGENVALUE/other/"
  save(p_values,file=paste(dirr,fileName,sep=""))
  return(  p_values)
}


N1=N2=25
num_runs=200
grid=seq( 0, 1, length.out = 100 )
contamination_level=0.1
cl=10

# fileName=paste0("other_dist_N_k_eigen_samesum_N_",N1,"_num_runs_200.Rda",sep="")

# b1=c(1,2,3)
# b2=c(3,2,1)
# runMVSim(both,out,contamination_level,N1,N2,b1,b2,grid,num_runs,fileName)
# 1:36

both=F
out=3

for(both in c(T,F)){
  for(out in 1:4){

    N1=N2=50
    
    
    grid=seq( 0, 1, length.out = 100 )
    
    
    
    b1=c(1,2,3)
    b2=c(3,2,1)
    
    fileName=paste0("cl_",cl,"_both_",both,"_out_",out,"_other_dist_N_k_eigen_1_samesum_N_",N1,"_num_runs_200.Rda",sep="")
    runMVSim(both,out,contamination_level,N1,N2,b1,b2,grid,num_runs,fileName)
    
    b1=1:11
    b2=11:1
    fileName=paste0("cl_",cl,"_both_",both,"_out_",out,"_other_dist_N_k_eigen_2_samesum_N_",N1,"_num_runs_200.Rda",sep="")
    runMVSim(both,out,contamination_level,N1,N2,b1,b2,grid,num_runs,fileName)
    
    b1=2^b1
    b2=2^b2
    fileName=paste0("cl_",cl,"_both_",both,"_out_",out,"_other_dist_N_k_eigen_3_samesum_N_",N1,"_num_runs_200.Rda",sep="")
    runMVSim(both,out,contamination_level,N1,N2,b1,b2,grid,num_runs,fileName)
    
    
    
    N1=N2=100
    
    b1=c(1,2,3)
    b2=c(3,2,1)
    
    fileName=paste0("cl_",cl,"_both_",both,"_out_",out,"_other_dist_N_k_eigen_1_samesum_N_",N1,"_num_runs_200.Rda",sep="")
    runMVSim(both,out,contamination_level,N1,N2,b1,b2,grid,num_runs,fileName)
    
    b1=1:11
    b2=11:1
    fileName=paste0("cl_",cl,"_both_",both,"_out_",out,"_other_dist_N_k_eigen_2_samesum_N_",N1,"_num_runs_200.Rda",sep="")
    runMVSim(both,out,contamination_level,N1,N2,b1,b2,grid,num_runs,fileName)
    
    b1=2^b1
    b2=2^b2
    fileName=paste0("cl_",cl,"_both_",both,"_out_",out,"_other_dist_N_k_eigen_3_samesum_N_",N1,"_num_runs_200.Rda",sep="")
    runMVSim(both,out,contamination_level,N1,N2,b1,b2,grid,num_runs,fileName)
    
    
    N1=N2=50
    
    b1=c(1,2,3)
    b2=b1*1.5
    
    fileName=paste0("cl_",cl,"_both_",both,"_out_",out,"_other_dist_N_depths_k_eigen_21_samesum_N_",N1,"_num_runs_200.Rda",sep="")
    runMVSim(both,out,contamination_level,N1,N2,b1,b2,grid,num_runs,fileName)
    
    b1=1:11
    b2=b1*1.5
    fileName=paste0("cl_",cl,"_both_",both,"_out_",out,"_other_dist_N_depths_k_eigen_22_samesum_N_",N1,"_num_runs_200.Rda",sep="")
    runMVSim(both,out,contamination_level,N1,N2,b1,b2,grid,num_runs,fileName)
    
    b1=2^b1
    b2=b1*1.5
    fileName=paste0("cl_",cl,"_both_",both,"_out_",out,"_other_dist_N_depths_k_eigen_23_samesum_N_",N1,"_num_runs_200.Rda",sep="")
    runMVSim(both,out,contamination_level,N1,N2,b1,b2,grid,num_runs,fileName)
    
    
    
    N1=N2=100
    
    b1=c(1,2,3)
    b2=b1*1.5
    
    fileName=paste0("cl_",cl,"_both_",both,"_out_",out,"_other_dist_N_depths_k_eigen_21_samesum_N_",N1,"_num_runs_200.Rda",sep="")
    runMVSim(both,out,contamination_level,N1,N2,b1,b2,grid,num_runs,fileName)
    
    b1=1:11
    b2=b1*1.5
    fileName=paste0("cl_",cl,"_both_",both,"_out_",out,"_other_dist_N_depths_k_eigen_22_samesum_N_",N1,"_num_runs_200.Rda",sep="")
    runMVSim(both,out,contamination_level,N1,N2,b1,b2,grid,num_runs,fileName)
    
    b1=2^b1
    b2=b1*1.5
    fileName=paste0("cl_",cl,"_both_",both,"_out_",out,"_other_dist_N_depths_k_eigen_23_samesum_N_",N1,"_num_runs_200.Rda",sep="")
    runMVSim(both,out,contamination_level,N1,N2,b1,b2,grid,num_runs,fileName)
    
    print(out/4)

  }
}

