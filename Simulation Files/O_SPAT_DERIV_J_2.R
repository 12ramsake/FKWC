#Source the initial and compute spatial scripts
#I am only doing 50 runs and 7 here. 

simulate_mv_depths<-function(N1,N2,c1,c2,res=100,
                             grid=seq( 0, 1, length.out = res ),
                             dist="t",
                             del=0.9){
  
  dat=gfd(N1,N2,Cov1=c1,Cov2=c2,res=res,dist=dist,del=del)
  
  
  fdat<-fdata(dat$mdata,dat$argvals)
  fdat<-fdata2fd(fdat)
  derivs<-derivcurves(dat$mdata)
  fderivs<-fdata(derivs,dat$argvals)
  fderivs<-fdata2fd(fderivs)
  
  FSD_depth=spatial_depth(fdat,N1+N2)+spatial_depth(fderivs,N1+N2)
  
  
  
  return(FSD_depth)
  
}


N1=N2=250
N1=N2=100
N1=N2=50
num_runs=50

#ready to go
#Simulation with the kernelized spatial depth
for(N1 in c(50,100,250)){
  N2=N1
  print("sim 1")
  #changes in beta
  count=0
  k1=function(t1,t2){k_exp(t1,t2,.5,.05)}
  c1=getCov(grid,k1)
  tmp=seq(0.5,1,length.out = 6)
  # for(dist in c("t")){
  for(dist in c("N","t","sn")){
    for(betai in 1:6){
      beta=tmp[betai]
      fileName=paste0("FKWCResults/ospat_dist_",dist,"_deriv_depths_k_exp_alpha1_p05_p05_beta_p",betai,"_p5_N_",N1,"_num_runs_",num_runs,".Rda",sep="")
      k2=function(t1,t2){k_exp(t1,t2,beta,.05)}
      c2=getCov(grid,k2)
      if(!file.exists(paste("",fileName,sep="")))
        if(!file.exists(paste("OSPAT_D/",fileName,sep="")))
          runMVSim(N1,N2,c1,c2,grid,num_runs,fileName,dist=dist)
      count=count+1
      print(count/(6*3))
    }
  }
  
  
  
  print("sim 2")
  #changes in alpha
  count=0
  k1=function(t1,t2){k_exp(t1,t2,.5,.05)}
  c1=getCov(grid,k1)
  tmp=seq(0.05,0.1,length.out=6)
  # for(dist in c("t")){
    for(dist in c("N","t","sn")){
    for(alphai in 2:6){
      alpha=tmp[alphai]
      fileName=paste0("FKWCResults/ospat_dist_",dist,"_deriv_depths_k_exp_alpha_p05_",alphai,
                      "_beta_p5_p5_N_",N1,"_num_runs_",num_runs,".Rda",sep="")
      k2=function(t1,t2){k_exp(t1,t2,0.5,alpha)}
      c2=getCov(grid,k2)
      if(!file.exists(paste("",fileName,sep="")))
        if(!file.exists(paste("OSPAT_D/",fileName,sep="")))
          runMVSim(N1,N2,c1,c2,grid,num_runs,fileName,dist=dist)
      count=count+1
      print(count/(6*3))
      
    }
  }
  
}
