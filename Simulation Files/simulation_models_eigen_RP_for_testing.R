
library(fChange)
library(EMMIXskew)
library(fda)
library(roahd)
library(mvtnorm)
library(fdasrvf)
library(doParallel)
library(MFHD)
library(mrfDepth)
library(doRNG)

gen_eigen_model<-function(N1,N2,nbasis=55,b1=rep(1,nbasis),b2=rep(1,nbasis)){
  
  
  dat1=fun_IID(N1, length(b1),Sigma=b1)
  dat2=fun_IID(N2, nbasis,Sigma=b2)
  dat=dat1
  dat$coefs=cbind(dat1$coefs,dat2$coefs)
  dat$fdnames$reps=c(dat$fdnames$reps,paste(rep("reps ",N2),as.character(((N1+1):(N1+N2)),collapse ="")))
  
  return(dat)
  
}

#simulate data and the resulting depths using the noderivative
#simulat
simulate_mv_depths<-function(N1,N2,b1,b2){
  
  fdat=gen_eigen_model(N1,N2,nbasis=length(b1),b1=b1,b2=b2)
  # fderivs=deriv.fd(fdat,1)
  #
  grid=seq(0,1,length.out=100)
  dat=t(eval.fd(grid,fdat))
  # ddat=t(eval.fd(grid, fderivs))
  # fdat<-fdata(dat$mdata,dat$argvals)
  # noderivs<-noderivcurves(dat$mdata)
  # fnoderivs<-fdata(noderivs,dat$argvals)
  
  RPD_depth=depth.RP(fdat,nproj=20,dfunc="Liu1")$dep

  
  depths=RPD_depth
  
  return(depths)
  
}




runMVSim<-function(N1,N2,b1,b2,grid,num_runs,fileName){
  
  
  no_cores<-detectCores()-1
  
  cl <- makeCluster(no_cores,type="FORK")   
  registerDoParallel(cl) 
  registerDoRNG(seed = 440)
  
  depth_values=try({foreach(i=1:num_runs) %dopar% {simulate_mv_depths(N1,N2,b1,b2)}})
  
  errorsp=inherits(depth_values, "try-error")
  if(errorsp){
    print("there was an error!")
    
  }
  
  stopCluster(cl)
  registerDoSEQ()
  closeAllConnections()
  
  dirr="/u/k3ramsay/ResearchDocuments/output/Functional Data Covariance Files/Two_Sample_Depths_Cov_Kernels/"
  save(depth_values,file=paste(dirr,fileName,sep=""))
  
}



N1=N2=50

num_runs=50
grid=seq( 0, 1, length.out = 100 )


N1=N2=2000
plot(simulate_mv_depths(N1,N2,b1=20:1,b2=1:20),type='l')

b1=c(1,2,3,4,5,6,7,8,9)
b2=c(9,8,7,6,5,4,3,2,1)
plot(simulate_mv_depths(N1,N2,b1,b2),type='l')
# fdat=gen_eigen_model(N1,N2,20,)



b1=c(1,2,3)
b2=c(3,2,1)



fileName=paste0("RPD_MBD/RP/RP_dist_N_noderiv_depths_k_eigen_1_samesum_N_",N1,"_ot_0_num_runs_50.Rda",sep="")
runMVSim(N1,N2,b1,b2,grid,num_runs,fileName)

b1=1:11
b2=11:1
fileName=paste0("RPD_MBD/RP/RP_dist_N_noderiv_depths_k_eigen_2_samesum_N_",N1,"_ot_0_num_runs_50.Rda",sep="")
runMVSim(N1,N2,b1,b2,grid,num_runs,fileName)

b1=2^b1
b2=2^b2
fileName=paste0("RPD_MBD/RP/RP_dist_N_noderiv_depths_k_eigen_3_samesum_N_",N1,"_ot_0_num_runs_50.Rda",sep="")
runMVSim(N1,N2,b1,b2,grid,num_runs,fileName)



N1=N2=100

b1=c(1,2,3)
b2=c(3,2,1)

fileName=paste0("RPD_MBD/RP/RP_dist_N_noderiv_depths_k_eigen_1_samesum_N_",N1,"_ot_0_num_runs_50.Rda",sep="")
runMVSim(N1,N2,b1,b2,grid,num_runs,fileName)

b1=1:11
b2=11:1
fileName=paste0("RPD_MBD/RP/RP_dist_N_noderiv_depths_k_eigen_2_samesum_N_",N1,"_ot_0_num_runs_50.Rda",sep="")
runMVSim(N1,N2,b1,b2,grid,num_runs,fileName)

b1=2^b1
b2=2^b2
fileName=paste0("RPD_MBD/RP/RP_dist_N_noderiv_depths_k_eigen_3_samesum_N_",N1,"_ot_0_num_runs_50.Rda",sep="")
runMVSim(N1,N2,b1,b2,grid,num_runs,fileName)

N1=N2=50

b1=c(1,2,3)
b2=b1*1.5

fileName=paste0("RPD_MBD/RP/RP_dist_N_noderiv_depths_k_eigen_21_samesum_N_",N1,"_ot_0_num_runs_50.Rda",sep="")
runMVSim(N1,N2,b1,b2,grid,num_runs,fileName)

b1=1:11
b2=b1*1.5
fileName=paste0("RPD_MBD/RP/RP_dist_N_noderiv_depths_k_eigen_22_samesum_N_",N1,"_ot_0_num_runs_50.Rda",sep="")
runMVSim(N1,N2,b1,b2,grid,num_runs,fileName)

b1=2^b1
b2=b1*1.5
fileName=paste0("RPD_MBD/RP/RP_dist_N_noderiv_depths_k_eigen_23_samesum_N_",N1,"_ot_0_num_runs_50.Rda",sep="")
runMVSim(N1,N2,b1,b2,grid,num_runs,fileName)



N1=N2=100

b1=c(1,2,3)
b2=b1*1.5

fileName=paste0("RPD_MBD/RP/RP_dist_N_noderiv_depths_k_eigen_21_samesum_N_",N1,"_ot_0_num_runs_50.Rda",sep="")
runMVSim(N1,N2,b1,b2,grid,num_runs,fileName)

b1=1:11
b2=b1*1.5
fileName=paste0("RPD_MBD/RP/RP_dist_N_noderiv_depths_k_eigen_22_samesum_N_",N1,"_ot_0_num_runs_50.Rda",sep="")
runMVSim(N1,N2,b1,b2,grid,num_runs,fileName)

b1=2^b1
b2=b1*1.5
fileName=paste0("RPD_MBD/RP/RP_dist_N_noderiv_depths_k_eigen_23_samesum_N_",N1,"_ot_0_num_runs_50.Rda",sep="")
runMVSim(N1,N2,b1,b2,grid,num_runs,fileName)



# 
# 
# 
# N1=N2=50
# 
# num_runs=50
# grid=seq( 0, 1, length.out = 100 )
# 
# 
# 
# b1=c(1,2,3)
# b2=c(3,2,1)
# 
# fileName=paste0("RPD_MBD/RP/RP_dist_N_noderiv_depths_k_eigen_1_samesum_N_",N1,"_ot_0_num_runs_50.Rda",sep="")
# runMVSim(N1,N2,b1,b2,grid,num_runs,fileName)
# 
# b1=1:11
# b2=11:1
# fileName=paste0("RPD_MBD/RP/RP_dist_N_noderiv_depths_k_eigen_2_samesum_N_",N1,"_ot_0_num_runs_50.Rda",sep="")
# runMVSim(N1,N2,b1,b2,grid,num_runs,fileName)
# 
# b1=2^b1
# b2=2^b2
# fileName=paste0("RPD_MBD/RP/RP_dist_N_noderiv_depths_k_eigen_3_samesum_N_",N1,"_ot_0_num_runs_50.Rda",sep="")
# runMVSim(N1,N2,b1,b2,grid,num_runs,fileName)
# 
# 
# 
# N1=N2=100
# 
# b1=c(1,2,3)
# b2=c(3,3,3)
# 
# fileName=paste0("RPD_MBD/RP/RP_dist_N_noderiv_depths_k_eigen_1_samesum_N_",N1,"_ot_0_num_runs_50.Rda",sep="")
# runMVSim(N1,N2,b1,b2,grid,num_runs,fileName)
# 
# b1=1:11
# b2=11:1
# fileName=paste0("RPD_MBD/RP/RP_dist_N_noderiv_depths_k_eigen_2_samesum_N_",N1,"_ot_0_num_runs_50.Rda",sep="")
# runMVSim(N1,N2,b1,b2,grid,num_runs,fileName)
# 
# b1=2^b1
# b2=2^b2
# fileName=paste0("RPD_MBD/RP/RP_dist_N_noderiv_depths_k_eigen_3_samesum_N_",N1,"_ot_0_num_runs_50.Rda",sep="")
# runMVSim(N1,N2,b1,b2,grid,num_runs,fileName)
# 
# 
# 
