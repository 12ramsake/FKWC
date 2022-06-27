#Source initial script

#simulate data and the resulting depths using the derivative
simulate_mv_depths<-function(N1,N2,N3,c1,c2,c3,res=100,
                             grid=seq( 0, 1, length.out = res ),
                             dist="t",
                             del=0.9){
  
  dat=gfd_k_3(N1,N2,N3,Cov1=c1,Cov2=c2,Cov3=c3,res=res,dist=dist,del=del)
  fdat<-fdata(dat$mdata,dat$argvals)

  # depth.RPD() implements a depth measure based on random projections possibly using several derivatives (see Cuevas et al. 2007).
  RPD_depth=depth.RP(fdat,dfunc = "Liu1")$dep
  #fraiman and muniz, really slow
  FM_depth=depth.FM(fdat)$dep
  #band depth
  MBD_depth=MBD(dat$mdata)
  LTR_depth=c(norm.fdata(fdat))
  
  depths=data.frame(FM_depth,RPD_depth,MBD_depth,LTR_depth)
  
  return(depths)

  
}

##
N1=N2=N3=50
N1=N2=N3=100

##
N1=N2=N3=50
N1=N2=N3=100

for(N1 in c(50,100)){
  for(ratio in c(0.4,0.6,0.8)){
  N2=floor(N1*ratio)
  N3=N1
  print("sim 1")
  #changes in beta
  count=0
  k1=function(t1,t2){k_exp(t1,t2,.5,.05)}
  c1=getCov(grid,k1)
  c3=getCov(grid,k1)
  tmp=seq(0.5,1,length.out = 15)
  for(dist in c("N","t","sn")){
  # for(dist in c("t")){
    for(betai in c(7,15)){
      beta=tmp[betai]
      fileName=paste0("FKWCResults/dist_",dist,"_noderiv_depths_k_exp_alpha1_p05_p05_p05_beta_p5_p",betai,"_p5_N1_",N1,"_N2_",N2,"_N3_",N3,"_num_runs_200.Rda",sep="")
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
    for(alphai in c(7,15)){
      alpha=tmp[alphai]
      fileName=paste0("FKWCResults/dist_",dist,"_noderiv_depths_k_exp_alpha1_p05_",alphai,"_p05_beta_p5_p5_p5_N1_",N1,"_N2_",N2,"_N3_",N3,"_num_runs_200.Rda",sep="")
      k2=function(t1,t2){k_exp(t1,t2,0.5,alpha)}
      c2=getCov(grid,k2)
      runMVSim(N1,N2,N3,c1,c2,c3,grid,num_runs,fileName,dist=dist)
      count=count+1
      print(count/(42))
  }
  }
  }
}









for(N1 in c(50,100)){
  for(ratio in c(0.4,0.6,0.8)){
    N2=floor(N1*ratio)
    N3=N1
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
        fileName=paste0("FKWCResults/dist_",dist,"_noderiv_depths_k_exp_alpha1_p05_p05_p05_beta_p5_p",betai,"_p5_N1_",N1,"_N2_",N2,"_N3_",N3,"_num_runs_200.Rda",sep="")
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
        fileName=paste0("FKWCResults/dist_",dist,"_noderiv_depths_k_exp_alpha1_p05_",alphai,"_p05_beta_p5_p5_p5_N1_",N1,"_N2_",N2,"_N3_",N3,"_num_runs_200.Rda",sep="")
        k2=function(t1,t2){k_exp(t1,t2,0.5,alpha)}
        c2=getCov(grid,k2)
        runMVSim(N1,N2,N3,c1,c2,c3,grid,num_runs,fileName,dist=dist)
        count=count+1
        print(count/(42))
      }
    }
  }
}













