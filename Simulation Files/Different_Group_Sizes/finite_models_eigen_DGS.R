# install.packages("fChange_0.2.1.tar.gz",type="source",repos=NULL)
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
  
  
  dat1=fun_IID(N1, nbasis,Sigma=b1)
  dat2=fun_IID(N2, nbasis,Sigma=b2)
  dat=dat1
  dat$coefs=cbind(dat1$coefs,dat2$coefs)
  dat$fdnames$reps=c(dat$fdnames$reps,paste(rep("reps ",N2),as.character(((N1+1):(N1+N2)),collapse ="")))
  
  return(dat)
  
}

#simulate data and the resulting depths using the derivative
#simulat
simulate_mv_depths<-function(N1,N2,b1,b2){
  
  fdat=gen_eigen_model(N1,N2,length(b1),b1,b2)
  fderivs=deriv.fd(fdat,1)
  #

  grid=seq(0,1,length.out=100)
  dat=t(eval.fd(grid,fdat))
  ddat=t(eval.fd(grid, fderivs))
  # fdat<-fdata(dat$mdata,dat$argvals)
  # derivs<-derivcurves(dat$mdata)
  # fderivs<-fdata(derivs,dat$argvals)
  

  
  fderivs<-fdata(ddat,grid)
  fdat=fdata(dat,grid)
  
  data_norms=cbind(c(norm.fdata(fdat)),c(norm.fdata(fderivs)))
  data_norms=t(apply(data_norms,1,'/',y=apply(data_norms,2,mad)))
  # data_norms
  LTR_depth=t(apply(data_norms,2,function(x){1/1+x}))
  LTR_depth=c(LTR_depth[1,]+LTR_depth[2,])
  kruskal.test(LTR_depth,g=c(rep(1,N1),rep(2,N2)))
  
  #fariman and muniz, really slow
  RPD_depth_like=depth.RPD(fdat)$dep
  RPD_depth_SD=depth.RPD(fdat,dfunc2='mdepth.SD')$dep
  FMp_depth=depth.FMp(list("f"=fdat,"f"=fderivs), dfunc = "mdepth.HS")$dep
  comb=mfData(grid,list(dat,ddat))
  MBD_depth=multiMBD(   comb, weights = 'uniform', manage_ties = FALSE )

  depths=data.frame(FMp_depth,RPD_depth_SD,MBD_depth,LTR_depth,RPD_depth_like)
  
  
  
  
  return(depths)
  
}




runMVSim<-function(N1,N2,b1,b2,grid,num_runs,fileName){
  

    no_cores<-detectCores()-1
    no_cores<-50
    
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
    
    dirr="/u/k3ramsay/ResearchDocuments/output/Functional Data Covariance Files/EIGENVALUE/"
    save(depth_values,file=paste(dirr,fileName,sep=""))
  
}



# N1=N2=50
num_runs=200
grid=seq( 0, 1, length.out = 100 )

for(N in c(50,100)*2){
  for(ratio in c(0.2,0.3,0.4)){
    N1=floor(ratio*N)
    N2=N-N1




    
    b1=c(1,2,3)
    b2=c(3,2,1)
    
    fileName=paste0("allD_dist_N_deriv_depths_k_eigen_1_samesum_N1_",N1,"_N2_",N2,"_num_runs_200.Rda",sep="")
    runMVSim(N1,N2,b1,b2,grid,num_runs,fileName)
    
    b1=1:11
    b2=11:1
    fileName=paste0("allD_dist_N_deriv_depths_k_eigen_2_samesum_N1_",N1,"_N2_",N2,"_num_runs_200.Rda",sep="")
    runMVSim(N1,N2,b1,b2,grid,num_runs,fileName)
    
    b1=2^b1
    b2=2^b2
    fileName=paste0("allD_dist_N_deriv_depths_k_eigen_3_samesum_N1_",N1,"_N2_",N2,"_num_runs_200.Rda",sep="")
    runMVSim(N1,N2,b1,b2,grid,num_runs,fileName)
    
    

    b1=c(1,2,3)
    b2=b1*1.5
    
    fileName=paste0("allD_dist_N_deriv_depths_k_eigen_21_samesum_N1_",N1,"_N2_",N2,"_num_runs_200.Rda",sep="")
    runMVSim(N1,N2,b1,b2,grid,num_runs,fileName)
    
    b1=1:11
    b2=b1*1.5
    fileName=paste0("allD_dist_N_deriv_depths_k_eigen_22_samesum_N1_",N1,"_N2_",N2,"_num_runs_200.Rda",sep="")
    runMVSim(N1,N2,b1,b2,grid,num_runs,fileName)
    
    b1=2^b1
    b2=b1*1.5
    fileName=paste0("allD_dist_N_deriv_depths_k_eigen_23_samesum_N1_",N1,"_N2_",N2,"_num_runs_200.Rda",sep="")
    runMVSim(N1,N2,b1,b2,grid,num_runs,fileName)

  }
  }
# 
# N1=N2=100
# 
# b1=c(1,2,3)
# b2=b1*1.5
# 
# fileName=paste0("allD_dist_N_deriv_depths_k_eigen_21_samesum__N1_",N1,"_N2_",N2,"_num_runs_200.Rda",sep="")
# runMVSim(N1,N2,b1,b2,grid,num_runs,fileName)
# 
# b1=1:11
# b2=b1*1.5
# fileName=paste0("allD_dist_N_deriv_depths_k_eigen_22_samesum__N1_",N1,"_N2_",N2,"_num_runs_200.Rda",sep="")
# runMVSim(N1,N2,b1,b2,grid,num_runs,fileName)
# 
# b1=2^b1
# b2=b1*1.5
# fileName=paste0("allD_dist_N_deriv_depths_k_eigen_23_samesum__N1_",N1,"_N2_",N2,"_num_runs_200.Rda",sep="")
# runMVSim(N1,N2,b1,b2,grid,num_runs,fileName)


