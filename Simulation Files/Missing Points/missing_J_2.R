#Source the initial script
#simulate data and the resulting depths using the derivative
simulate_mv_depths<-function(N1,N2,c1,c2,res=100,
                             grid=seq( 0, 1, length.out = res ),
                             dist="t",
                             del=0.9){
  
  dat=gfd(N1,N2,Cov1=c1,Cov2=c2,res=res,dist=dist,del=del)
  dat$mdata=apply(dat$mdata,2,modify_curve,amount_to_delete=0.2)
  # matplot(t(dat$mdata)[,1:5],type='l')

  aa=apply(dat$mdata, 1, zoo::na.spline)
  dat$mdata=t(aa)
  
  
  
  # b1 <- create.bspline.basis(rangeval=c(0,1),nbasis=11, norder=4)
  # tm=smooth.basis.sparse(grid,t(dat$mdata),fdPar(fdobj=b1, Lfdobj=1, lambda=1))
  # plot( dat$mdata[2,],type='l')
  # plot( zoo::na.approx(dat$mdata[2,]),type='l')
  # filled_in=apply(dat$mdata,1  ,zoo::na.approx)
  # if (is.null(dim( filled_in)))
  #   dat$mdata=do.call(cbind,filled_in)
  # plot( filled_in[2,])
  # matplot(t(dat$mdata),type='l')
  # colSums(apply(filled_in,2,is.na))
  # dat$mdata=filled_in
  # matplot(dat$mdata[,1:5],type='l')
  #NOW SMOOTH AND REDISC
  # b1 <- create.bspline.basis(nbasis=11, norder=4)
  # smoothed=Data2fd(grid,t(dat$mdata),basisobj = b1)
  
  # dat$mdata=t(fda::eval.fd(grid,smoothed))
  
  fdat<-fdata(dat$mdata,dat$argvals)
  # fdata.deriv( fdat)
  derivs<-derivcurves(dat$mdata)
  fderivs<-fdata(derivs,dat$argvals)
  # depth.RPD() implements a depth measure based on random projections possibly using several derivatives (see Cuevas et al. 2007).
  RPD_depth=depth.RPD(fdat)$dep
  RPD_depth_SD=depth.RPD(fdat, dfunc2='mdepth.SD')$dep
  #fraiman and muniz, really slow
  FMp_depth=depth.FMp(list(fdat,fderivs), dfunc = "mdepth.HS")$dep
  MBD_depth=multiMBD(list(dat$mdata,derivs), weights = 'uniform', manage_ties = FALSE )
  #LTR
  data_norms=cbind(c(norm.fdata(fdat)),c(norm.fdata(fderivs)))
  data_norms=t(apply(data_norms,1,'/',y=apply(data_norms,2,mad)))
  # data_norms
  LTR_depth=t(apply(data_norms,2,function(x){1/1+x}))
  LTR_depth=c(LTR_depth[1,]+LTR_depth[2,])
  
  depths=data.frame(FMp_depth,RPD_depth_SD,MBD_depth,LTR_depth,RPD_depth)
  
  return(depths)
  
}



# N1=N2=250
# N1=N2=100
# N1=N2=50


N1=100
  N2=N1
  print("sim 1")
  #changes in beta
  count=0
  k1=function(t1,t2){k_exp(t1,t2,.5,.05)}
  c1=getCov(grid,k1)
  tmp=seq(0.5,1,length.out = 15)
  for(dist in c("N","t")){
    for(betai in 1:15){
      beta=tmp[betai]
      fileName=paste0("FKWCResults/dist_",dist,"_missing_k_exp_alpha1_p05_p05_beta_p",betai,"_p5_N_",N1,"_num_runs_200.Rda",sep="")
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
  for(dist in c("N","t")){
    for(alphai in 1:15){
      alpha=tmp[alphai]
      fileName=paste0("FKWCResults/dist_",dist,"_missing_k_exp_alpha_p05_",alphai,"_beta_p5_p5_N_",N1,"_num_runs_200.Rda",sep="")
      k2=function(t1,t2){k_exp(t1,t2,0.5,alpha)}
      c2=getCov(grid,k2)
      runMVSim(N1,N2,c1,c2,grid,num_runs,fileName,dist=dist)
      count=count+1
      print(count/(14*3))
      
    }
  }



  
  
  
  
  
  
  
  
  
  # 
  # 
  # dat=gfd(N1,N2,Cov1=c1,Cov2=c2,res=res,dist=dist,del=del)
  # dat$mdata=apply(dat$mdata,2,modify_curve,amount_to_delete=0.2)
  # matplot(t(dat$mdata)[,1:5],type='l')
  # aa=apply(dat$mdata, 1, zoo::na.spline)
  # zoo::na.spline(dat$mdata[2,])
  # 
  # b1 <- create.bspline.basis(rangeval=c(0,1),nbasis=11, norder=4)
  # Data2fd(grid,t(dat$mdata))
  # tm=  smooth.basis.sparse(grid,t(dat$mdata),fdPar(fdobj=b1, Lfdobj=1, lambda=1))
  # 
  # 
 
  


