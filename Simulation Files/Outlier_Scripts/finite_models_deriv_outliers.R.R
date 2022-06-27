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

generate_outliers=function(dat,N1,N2,contamination_level=0.1,outlier_type=1,both_samples=T){
  
  con1=sample(1:N1,contamination_level*(N1+N2))
  if(both_samples)
    con2=sample((N1+1):(N1+N2),contamination_level*(N1+N2))
  else
    con2=NULL
  orig=dat$mdata[c(con1,con2),]
  
  if(outlier_type==1)
    dat$mdata[c(con1,con2),]=t(apply(t(orig)+2*grid,2,rev))
  
  ##model 2
  if(outlier_type==2){
    if(both)
      mod_2=replicate(contamination_level*(N1+N2)*2,0.2*dnorm(grid,runif(1,.25,.75),.2))
    else
      mod_2=replicate(contamination_level*(N1+N2),0.2*dnorm(grid,runif(1,.25,.75),.2))
    dat$mdata[c(con1,con2),]=t(mod_2)+orig
  }
  
  ##model 3
  if(outlier_type==3){
    if(both)
      mod_3=replicate(contamination_level*(N1+N2)*2,sin(8*(grid + runif(1,.25,.75))*pi))
    else
      mod_3=replicate(contamination_level*(N1+N2),sin(8*(grid + runif(1,.25,.75))*pi))
    dat$mdata[c(con1,con2),]=t(mod_3)+orig*.2
  }
  ##model 4
  if(outlier_type==4){
    if(both)
      mod_4=replicate(contamination_level*(N1+N2)*2,.25*dnorm(grid,runif(1,.25,.75),.01))
    else 
      mod_4=replicate(contamination_level*(N1+N2),.25*dnorm(grid,runif(1,.25,.75),.01))
    dat$mdata[c(con1,con2),]=t(mod_4)+orig*0.5
  }
  return(dat)
}
#simulate data and the resulting depths using the derivative
#simulat
simulate_mv_depths<-function(both,out,contamination_level,N1,N2,b1,b2){
  
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
  #
  # FSD_depth=spatial_depth(fdat,N1+N2)+spatial_depth(fderivs,N1+N2)
  # KFSD_depth=K_spatial_depth(fdat,N1+N2)+K_spatial_depth(fderivs,N1+N2)
  
  grid=seq(0,1,length.out=100)
  dat=t(eval.fd(grid,fdat))
  ddat=t(eval.fd(grid, fderivs))
  matplot(t(dat),type='l')
  # fdat<-fdata(dat$mdata,dat$argvals)
  # derivs<-derivcurves(dat$mdata)
  # fderivs<-fdata(derivs,dat$argvals)
  

  
  fderivs<-fdata(ddat,grid)
  fdat=fdata(dat,grid)
  #fariman and muniz, really slow
  RPD_depth_SD=depth.RPD(fdat,dfunc2='mdepth.SD')$dep
  RPD_depth_Like=depth.RPD(fdat)$dep
  FMp_depth=depth.FMp(list("f"=fdat,"f"=fderivs), dfunc = "mdepth.HS")$dep
  comb=mfData(grid,list(dat,ddat))
  MBD_depth=multiMBD( comb, weights = 'uniform', manage_ties = FALSE )
  data_norms=cbind(c(norm.fdata(fdat)),c(norm.fdata(fderivs)))
  data_norms=t(apply(data_norms,1,'/',y=apply(data_norms,2,mad)))
  # data_norms
  LTR_depth=t(apply(data_norms,2,function(x){1/1+x}))
  LTR_depth=c(LTR_depth[1,]+LTR_depth[2,])

  depths=data.frame(FMp_depth,RPD_depth_SD,MBD_depth,LTR_depth,RPD_depth_Like)
  
  return(depths)
  
}




runMVSim<-function(both,out,contamination_level,N1,N2,b1,b2,grid,num_runs,fileName){
  

    no_cores<-detectCores()-1
    # no_cores<-40
    no_cores<-60
    
    cl <- makeCluster(no_cores,type="FORK")   
    registerDoParallel(cl) 
    registerDoRNG(seed = 440)
    
    depth_values=try({foreach(i=1:num_runs) %dopar% {simulate_mv_depths(both,out,contamination_level,N1,N2,b1,b2)}})
    
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

contamination_level=0.1
cl=contamination_level*100



num_runs=200
grid=seq( 0, 1, length.out = 100 )

for(both in c(T,F)){
  for(out in 1:4){

    N1=N2=50
    
    b1=c(1,2,3)
    b2=c(3,2,1)
    
    fileName=paste0("cl_",cl,"_both_",both,"_out_",out,"_allD_dist_N_deriv_depths_k_eigen_1_samesum_N_",N1,"_num_runs_200.Rda",sep="")
    runMVSim(both,out,contamination_level,N1,N2,b1,b2,grid,num_runs,fileName)
    
    b1=1:11
    b2=11:1
    fileName=paste0("cl_",cl,"_both_",both,"_out_",out,"_allD_dist_N_deriv_depths_k_eigen_2_samesum_N_",N1,"_num_runs_200.Rda",sep="")
    runMVSim(both,out,contamination_level,N1,N2,b1,b2,grid,num_runs,fileName)
    
    b1=2^b1
    b2=2^b2
    fileName=paste0("cl_",cl,"_both_",both,"_out_",out,"_allD_dist_N_deriv_depths_k_eigen_3_samesum_N_",N1,"_num_runs_200.Rda",sep="")
    runMVSim(both,out,contamination_level,N1,N2,b1,b2,grid,num_runs,fileName)
    
    
    
    N1=N2=100
    
    b1=c(1,2,3)
    b2=c(3,2,1)
    
    fileName=paste0("cl_",cl,"_both_",both,"_out_",out,"_allD_dist_N_deriv_depths_k_eigen_1_samesum_N_",N1,"_num_runs_200.Rda",sep="")
    runMVSim(both,out,contamination_level,N1,N2,b1,b2,grid,num_runs,fileName)
    
    b1=1:11
    b2=11:1
    fileName=paste0("cl_",cl,"_both_",both,"_out_",out,"_allD_dist_N_deriv_depths_k_eigen_2_samesum_N_",N1,"_num_runs_200.Rda",sep="")
    runMVSim(both,out,contamination_level,N1,N2,b1,b2,grid,num_runs,fileName)
    
    b1=2^b1
    b2=2^b2
    fileName=paste0("cl_",cl,"_both_",both,"_out_",out,"_allD_dist_N_deriv_depths_k_eigen_3_samesum_N_",N1,"_num_runs_200.Rda",sep="")
    runMVSim(both,out,contamination_level,N1,N2,b1,b2,grid,num_runs,fileName)
    
    
    N1=N2=50
    
    b1=c(1,2,3)
    b2=b1*1.5
    
    fileName=paste0("cl_",cl,"_both_",both,"_out_",out,"_allD_dist_N_deriv_depths_k_eigen_21_samesum_N_",N1,"_num_runs_200.Rda",sep="")
    runMVSim(both,out,contamination_level,N1,N2,b1,b2,grid,num_runs,fileName)
    
    b1=1:11
    b2=b1*1.5
    fileName=paste0("cl_",cl,"_both_",both,"_out_",out,"_allD_dist_N_deriv_depths_k_eigen_22_samesum_N_",N1,"_num_runs_200.Rda",sep="")
    runMVSim(both,out,contamination_level,N1,N2,b1,b2,grid,num_runs,fileName)
    
    b1=2^b1
    b2=b1*1.5
    fileName=paste0("cl_",cl,"_both_",both,"_out_",out,"_allD_dist_N_deriv_depths_k_eigen_23_samesum_N_",N1,"_num_runs_200.Rda",sep="")
    runMVSim(both,out,contamination_level,N1,N2,b1,b2,grid,num_runs,fileName)
    
    
    
    N1=N2=100
    
    b1=c(1,2,3)
    b2=b1*1.5
    
    fileName=paste0("cl_",cl,"_both_",both,"_out_",out,"_allD_dist_N_deriv_depths_k_eigen_21_samesum_N_",N1,"_num_runs_200.Rda",sep="")
    runMVSim(both,out,contamination_level,N1,N2,b1,b2,grid,num_runs,fileName)
    
    b1=1:11
    b2=b1*1.5
    fileName=paste0("cl_",cl,"_both_",both,"_out_",out,"_allD_dist_N_deriv_depths_k_eigen_22_samesum_N_",N1,"_num_runs_200.Rda",sep="")
    runMVSim(both,out,contamination_level,N1,N2,b1,b2,grid,num_runs,fileName)
    
    b1=2^b1
    b2=b1*1.5
    fileName=paste0("cl_",cl,"_both_",both,"_out_",out,"_allD_dist_N_deriv_depths_k_eigen_23_samesum_N_",N1,"_num_runs_200.Rda",sep="")
    runMVSim(both,out,contamination_level,N1,N2,b1,b2,grid,num_runs,fileName)

  }
}
