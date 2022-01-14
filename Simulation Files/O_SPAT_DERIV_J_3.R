#Source initial and spatial scripts


simulate_mv_depths<-function(N1=10,N2=10,N3=10,
                             Cov1,
                             Cov2,
                             Cov3,
                             res=1e2, 
                             grid=seq( 0, 1, length.out = res ),
                             dist="t",
                             del=0.9){
  
  
  dat=gfd_k_3(N1,N2,N3,Cov1=c1,Cov2=c2,Cov3=c3,res=res,dist=dist,del=del)
  
  fdat<-fdata(dat$mdata,dat$argvals)
  fdat<-fdata2fd(fdat)
  derivs<-derivcurves(dat$mdata)
  fderivs<-fdata(derivs,dat$argvals)
  fderivs<-fdata2fd(fderivs)
  
  FSD_depth=spatial_depth(fdat,N1+N2+N3)+spatial_depth(fderivs,N1+N2+N3)
  
  
  
  return(FSD_depth)
  
  
}



N1=N2=N3=50
N1=N2=N3=100

for(N1 in c(50,100)){
  N2=N1
  N3=N1
  print("sim 1")
  #changes in beta
  count=0
  k1=function(t1,t2){k_exp(t1,t2,.5,.05)}
  c1=getCov(grid,k1)
  c3=getCov(grid,k1)
  tmp=seq(0.5,1,length.out = 15)
  for(dist in c("N","t","sn")){
    for(betai in 1:15){
      beta=tmp[betai]
      fileName=paste0("FKWCResults/ospat_dist_",dist,"_deriv_depths_k_exp_alpha1_p05_p05_p05_beta_p5_p",betai,"_p5_N_",N1,"_num_runs_200.Rda",sep="")
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
    for(alphai in 1:15){
      alpha=tmp[alphai]
      fileName=paste0("FKWCResults/ospat_dist_",dist,"_deriv_depths_k_exp_alpha1_p05_",alphai,"_p05_beta_p5_p5_p5_N_",N1,"_num_runs_200.Rda",sep="")
      k2=function(t1,t2){k_exp(t1,t2,0.5,alpha)}
      c2=getCov(grid,k2)
      runMVSim(N1,N2,N3,c1,c2,c3,grid,num_runs,fileName,dist=dist)
      count=count+1
      print(count/(42))
    }
  }
  
}