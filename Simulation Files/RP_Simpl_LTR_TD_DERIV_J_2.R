#Source the initial script
#simulate data and the resulting depths using the derivative
simulate_mv_depths<-function(N1,N2,c1,c2,res=100,
                             grid=seq( 0, 1, length.out = res ),
                             dist="t",
                             del=0.9){
  
  dat=gfd(N1,N2,Cov1=c1,Cov2=c2,res=res,dist=dist,del=del)
  
  
  fdat<-fdata(dat$mdata,dat$argvals)
  derivs<-derivcurves(dat$mdata)
  fderivs<-fdata(derivs,dat$argvals)
  # depth.RPD() implements a depth measure based on random projections possibly using several derivatives (see Cuevas et al. 2007).
  RPD_depth_SD=depth.RPD(fdat,dfunc2='mdepth.SD')$dep
  #fraiman and muniz, really slow
  #LTR
  data_norms=cbind(c(norm.fdata(fdat)),c(norm.fdata(fderivs)))
  data_norms=t(apply(data_norms,1,'/',y=apply(data_norms,2,mad)))
  # data_norms
  LTR_depth=t(apply(data_norms,2,function(x){1/1+x}))
  LTR_depth=c(LTR_depth[1,]+LTR_depth[2,])
  
  depths=data.frame(RPD_depth_SD,LTR_depth_HD)
  
  return(depths)
  
}



# N1=N2=250
# N1=N2=100
# N1=N2=50


for(N1 in c(50,100,250)){
  N2=N1
  print("sim 1")
  #changes in beta
  count=0
  k1=function(t1,t2){k_exp(t1,t2,.5,.05)}
  c1=getCov(grid,k1)
  tmp=seq(0.5,1,length.out = 15)
  for(dist in c("N","t","sn")){
    for(betai in 1:15){
      beta=tmp[betai]
      fileName=paste0("FKWCResults/dist_",dist,"_deriv_SIM_TD_k_exp_alpha1_p05_p05_beta_p",betai,"_p5_N_",N1,"_num_runs_200.Rda",sep="")
      k2=function(t1,t2){k_exp(t1,t2,beta,.05)}
      c2=getCov(grid,k2)
      runMVSim(N1,N2,c1,c2,grid,num_runs,fileName,dist=dist)
      count=count+1
      print(count/(14*3))
    }
  }
  
  
  
  print("sim 2")
  #changes in alpha
  count=0
  k1=function(t1,t2){k_exp(t1,t2,.5,.05)}
  c1=getCov(grid,k1)
  tmp=seq(0.05,0.1,length.out=15)
  for(dist in c("N","t","sn")){
    for(alphai in 1:15){
      alpha=tmp[alphai]
      fileName=paste0("FKWCResults/dist_",dist,"_deriv_SIM_TD_k_exp_alpha_p05_",alphai,"_beta_p5_p5_N_",N1,"_num_runs_200.Rda",sep="")
      k2=function(t1,t2){k_exp(t1,t2,0.5,alpha)}
      c2=getCov(grid,k2)
      runMVSim(N1,N2,c1,c2,grid,num_runs,fileName,dist=dist)
      count=count+1
      print(count/(14*3))
      
    }
  }
}







# 
# 
# 
# for(N1 in c(50,100,250)){
#   N2=N1
#   print("sim 1")
#   #changes in beta
#   count=0
#   k1=function(t1,t2){k_exp(t1,t2,.5,.05)}
#   c1=getCov(grid,k1)
#   tmp=seq(0.5,1,length.out = 15)
#   for(dist in c("N","t","sn")){
#     for(betai in 1:1){
#       beta=tmp[betai]
#       fileName=paste0("FKWCResults/dist_",dist,"_deriv_SIM_TD_k_exp_alpha1_p05_p05_beta_p",betai,"_p5_N_",N1,"_num_runs_200.Rda",sep="")
#       k2=function(t1,t2){k_exp(t1,t2,beta,.05)}
#       c2=getCov(grid,k2)
#       runMVSim(N1,N2,c1,c2,grid,num_runs,fileName,dist=dist)
#       count=count+1
#       print(count/(14*3))
#     }
#   }
#   
#   
#   
#   print("sim 2")
#   #changes in alpha
#   count=0
#   k1=function(t1,t2){k_exp(t1,t2,.5,.05)}
#   c1=getCov(grid,k1)
#   tmp=seq(0.05,0.1,length.out=15)
#   for(dist in c("N","t","sn")){
#     for(alphai in 1:1){
#       alpha=tmp[alphai]
#       fileName=paste0("FKWCResults/dist_",dist,"_deriv_SIM_TD_k_exp_alpha_p05_",alphai,"_beta_p5_p5_N_",N1,"_num_runs_200.Rda",sep="")
#       k2=function(t1,t2){k_exp(t1,t2,0.5,alpha)}
#       c2=getCov(grid,k2)
#       runMVSim(N1,N2,c1,c2,grid,num_runs,fileName,dist=dist)
#       count=count+1
#       print(count/(14*3))
#       
#     }
#   }
# }














