
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

#load in Guo and Boente methods

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
gfd<-function(N1=10,N2=10,
              res=1e2, 
              Cov1,
              Cov2,
              dist="t",
              del=0.9){
  
  
  grid = seq( 0, 1, length.out = res )
  
  
  if(dist=="N"){
    Data1 = generate_gauss_fdata( N1,  rep(0,res),Cov = Cov1 )
    Data2 = generate_gauss_fdata( N2, rep(0,res), Cov = Cov2 )
  }
  else if(dist=="t"){
    Data1 = rmvt( N1, sigma =  Cov1/3,df=3 )
    Data2 = rmvt( N2, sigma =  Cov2/3,df=3 )
  }
  else{
    cp1=list(mean=rep(0,nrow(Cov1)), var.cov=Cov1, gamma1=rep(del,nrow(Cov1))/nrow(Cov1))
    cp2=list(mean=rep(0,nrow(Cov2)), var.cov=Cov2, gamma1=rep(del,nrow(Cov2))/nrow(Cov2))
    dp1=cp2dp(cp1, "SN")
    dp2=cp2dp(cp2, "SN")
    
    Data1 = rmsn(N1, dp=dp1)
    Data2 = rmsn(N2, dp=dp2)
  }
  
  return(list(argvals =  grid, mdata=rbind(Data1,Data2)))
  
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


oneRun<-function(N1,N2,c1,c2,res=100,
                 grid=seq( 0, 1, length.out = res ),
                 dist="t",
                 del=0.9,
                 outlier_type,both,contamination_level){
  
  dat=gfd(N1,N2,Cov1=c1,Cov2=c2,res=res,dist=dist,del=del)
  dat=generate_outliers(dat,N1,N2,contamination_level=contamination_level,outlier_type=outlier_type,both_samples=both)
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

runMVSim<-function(both,out,contamination_level,N1,N2,c1,c2,grid=seq( 0, 1, length.out = 100 ),num_runs=200,fileName,dist="t",del=0.9){
  
  
  
  # no_cores<-detectCores()-1
  # no_cores<-floor(detectCores()*0.8)
  no_cores<-40
  
  
  cl <- makeCluster(no_cores,type="FORK")
  # cl <- makeCluster(no_cores)
  registerDoParallel(cl)
  registerDoRNG(seed = 440)
  # set.seed(440)
  p_values=try({foreach(i=1:num_runs) %dopar% {oneRun(N1,N2,c1,c2,res=100,del=del,dist=dist,grid=grid,
                                                      both=both,outlier_type=out,contamination_level=contamination_level)}})
  
  errorsp=inherits(p_values, "try-error")
  if(errorsp){
    print("there was an error!")
    
  }
  
  stopCluster(cl)
  registerDoSEQ()
  closeAllConnections()
  
  # save(p_values,file=fileName)
  dirr="/u/k3ramsay/ResearchDocuments/output/Functional Data Covariance Files/OTHER_TEST_RESULTS/"
  save(p_values,file=paste(dirr,fileName,sep=""))
  # save(p_values,file=paste(dirr,fileName,sep=""))
}



num_runs=200
grid=seq( 0, 1, length.out = 100 )


contamination_level=0.1
cl=contamination_level*100

contamination_level=0.05
cl=contamination_level*100





N1=100
both=FALSE
both=TRUE

for(contamination_level in c(0.01,0.025,0.05)){
  cl=contamination_level*100

  for(N1 in 100){
    for(both in c(T,F)){
      for(out in 4){
        N2=N1
        print("sim 1")
        #changes in beta
        count=0
        k1=function(t1,t2){k_exp(t1,t2,.5,.05)}
        c1=getCov(grid,k1)
        tmp=seq(0.5,1,length.out = 15)
        # for(dist in c("t")){
        for(dist in c("N")){
          for(betai in 1:15){
            beta=tmp[betai]
            fileName=paste0("cl_",cl,"_both_",both,"_out_",out,"dist_",dist,"_other_k_exp_alpha1_p05_p05_beta_p",betai,"_p5_N_",N1,"_num_runs_200.Rda",sep="")
            k2=function(t1,t2){k_exp(t1,t2,beta,.05)}
            c2=getCov(grid,k2)
            runMVSim(both,out,contamination_level,N1,N2,c1,c2,grid,num_runs,fileName,dist=dist)
            count=count+1
            print(count/(14*3))
          }
        }
        
        
        ##
        # N1=N2=100
        # N1=N2=50
        print("sim 2")
        #changes in alpha
        count=0
        k1=function(t1,t2){k_exp(t1,t2,.5,.05)}
        c1=getCov(grid,k1)
        tmp=seq(0.05,0.1,length.out=15)
        # for(dist in c("t")){
        for(dist in c("N")){
          for(alphai in 1:15){
            alpha=tmp[alphai]
            fileName=paste0("cl_",cl,"_both_",both,"_out_",out,"dist_",dist,"_other_k_exp_alpha1_p05_",
                            alphai,"_beta_p5_p5_N_",N1,"_num_runs_200.Rda",sep="")
            k2=function(t1,t2){k_exp(t1,t2,0.5,alpha)}
            grid=seq( 0, 1, length.out = 100 )
            c2=getCov(grid,k2)
            runMVSim(both,out,contamination_level,N1,N2,c1,c2,grid,num_runs,fileName,dist=dist)
            count=count+1
            print(count/(14*3))
          }
        }
        #just guo
      }
    }
  } 
  
  
  
}
  
  
  
  
  

for(N1 in c(50,100,250)){
  for(both in c(T,F)){
    for(out in 1:4){
    N2=N1
    print("sim 1")
    #changes in beta
    count=0
    k1=function(t1,t2){k_exp(t1,t2,.5,.05)}
    c1=getCov(grid,k1)
    tmp=seq(0.5,1,length.out = 15)
    # for(dist in c("t")){
    for(dist in c("N")){
      for(betai in 1:15){
        beta=tmp[betai]
        fileName=paste0("cl_",cl,"_both_",both,"_out_",out,"dist_",dist,"_other_k_exp_alpha1_p05_p05_beta_p",betai,"_p5_N_",N1,"_num_runs_200.Rda",sep="")
        k2=function(t1,t2){k_exp(t1,t2,beta,.05)}
        c2=getCov(grid,k2)
        runMVSim(both,out,contamination_level,N1,N2,c1,c2,grid,num_runs,fileName,dist=dist)
        count=count+1
        print(count/(14*3))
      }
    }
  
  
  ##
  # N1=N2=100
  # N1=N2=50
    print("sim 2")
    #changes in alpha
    count=0
    k1=function(t1,t2){k_exp(t1,t2,.5,.05)}
    c1=getCov(grid,k1)
    tmp=seq(0.05,0.1,length.out=15)
    # for(dist in c("t")){
    for(dist in c("N")){
      for(alphai in 1:15){
        alpha=tmp[alphai]
        fileName=paste0("cl_",cl,"_both_",both,"_out_",out,"dist_",dist,"_other_k_exp_alpha1_p05_",
                        alphai,"_beta_p5_p5_N_",N1,"_num_runs_200.Rda",sep="")
        k2=function(t1,t2){k_exp(t1,t2,0.5,alpha)}
        grid=seq( 0, 1, length.out = 100 )
        c2=getCov(grid,k2)
        runMVSim(both,out,contamination_level,N1,N2,c1,c2,grid,num_runs,fileName,dist=dist)
        count=count+1
        print(count/(14*3))
      }
    }
  #just guo
    }
    }
}







both=F


for(N1 in c(100)){
  # for(both in c(T,F)){
    # for(out in 1:4){
      for(out in 3:4){
  N2=N1
  print("sim 1")
  #changes in beta
  count=0
  k1=function(t1,t2){k_exp(t1,t2,.5,.05)}
  c1=getCov(grid,k1)
  tmp=seq(0.5,1,length.out = 15)
  # for(dist in c("t")){
  for(dist in c("N")){
    for(betai in 2:15){
      beta=tmp[betai]
      fileName=paste0("cl_",cl,"_both_",both,"_out_",out,
                      "dist_",dist,"_other_k_exp_alpha1_p05_p05_beta_p",betai,"_p5_N_",N1,"_num_runs_200.Rda",sep="")
      k2=function(t1,t2){k_exp(t1,t2,beta,.05)}
      c2=getCov(grid,k2)
      runMVSim(both,out,contamination_level,N1,N2,c1,c2,grid,num_runs,fileName,dist=dist)
      count=count+1
      print(count/(14*3))
    }
  }
  
  
  ##
  # N1=N2=100
  # N1=N2=50
  print("sim 2")
  #changes in alpha
  count=0
  k1=function(t1,t2){k_exp(t1,t2,.5,.05)}
  c1=getCov(grid,k1)
  tmp=seq(0.05,0.1,length.out=15)
  # for(dist in c("t")){
  for(dist in c("N")){
    for(alphai in 5){
      alpha=tmp[alphai]
      fileName=paste0("cl_",cl,"_both_",both,"_out_",out,"dist_",dist,
                      "_other_k_exp_alpha1_p05_",alphai,"_beta_p5_p5_N_",N1,"_num_runs_200.Rda",sep="")
      k2=function(t1,t2){k_exp(t1,t2,0.5,alpha)}
      grid=seq( 0, 1, length.out = 100 )
      c2=getCov(grid,k2)
      runMVSim(both,out,contamination_level,N1,N2,c1,c2,grid,num_runs,fileName,dist=dist)
      count=count+1
      print(count/(14*3))
    }
    }
  #just guo
  }
}
# }




















for(N1 in c(50,100,250)){
 both=F
 for(out in 2:4){
      N2=N1
      print("sim 1")
      #changes in beta
      count=0
      k1=function(t1,t2){k_exp(t1,t2,.5,.05)}
      c1=getCov(grid,k1)
      tmp=seq(0.5,1,length.out = 15)
      # for(dist in c("t")){
      dist='N'
        for(betai in 2:15){
          beta=tmp[betai]
          fileName=paste0("cl_",cl,"_both_",both,"_out_",out,"dist_",dist,"_other_k_exp_alpha1_p05_p05_beta_p",betai,"_p5_N_",N1,"_num_runs_200.Rda",sep="")
          k2=function(t1,t2){k_exp(t1,t2,beta,.05)}
          c2=getCov(grid,k2)
          runMVSim(both,out,contamination_level,N1,N2,c1,c2,grid,num_runs,fileName,dist=dist)
          count=count+1
          print(count/(14*3))
        }
      
      
      
      ##
      # N1=N2=100
      # N1=N2=50
      print("sim 2")
      #changes in alpha
      count=0
      k1=function(t1,t2){k_exp(t1,t2,.5,.05)}
      c1=getCov(grid,k1)
      tmp=seq(0.05,0.1,length.out=15)
      # for(dist in c("t")){
      dist='N'
        for(alphai in 1:15){
          alpha=tmp[alphai]
          fileName=paste0("cl_",cl,"_both_",both,"_out_",out,"dist_",dist,"_other_k_exp_alpha1_p05_",
                          alphai,"_beta_p5_p5_N_",N1,"_num_runs_200.Rda",sep="")
          k2=function(t1,t2){k_exp(t1,t2,0.5,alpha)}
          grid=seq( 0, 1, length.out = 100 )
          c2=getCov(grid,k2)
          runMVSim(both,out,contamination_level,N1,N2,c1,c2,grid,num_runs,fileName,dist=dist)
          count=count+1
          print(count/(14*3))
        }
      
      #just guo
    }
  
}








for(N1 in c(250)){
  both=T
  for(out in 1:4){
    N2=N1
    print("sim 1")
    #changes in beta
    count=0
    k1=function(t1,t2){k_exp(t1,t2,.5,.05)}
    c1=getCov(grid,k1)
    tmp=seq(0.5,1,length.out = 15)
    # for(dist in c("t")){
    dist='N'
    for(betai in 2:15){
      beta=tmp[betai]
      fileName=paste0("cl_",cl,"_both_",both,"_out_",out,"dist_",dist,"_other_k_exp_alpha1_p05_p05_beta_p",betai,"_p5_N_",N1,"_num_runs_200.Rda",sep="")
      k2=function(t1,t2){k_exp(t1,t2,beta,.05)}
      c2=getCov(grid,k2)
      runMVSim(both,out,contamination_level,N1,N2,c1,c2,grid,num_runs,fileName,dist=dist)
      count=count+1
      print(count/(14*3))
    }
    
    
    
    ##
    # N1=N2=100
    # N1=N2=50
    print("sim 2")
    #changes in alpha
    count=0
    k1=function(t1,t2){k_exp(t1,t2,.5,.05)}
    c1=getCov(grid,k1)
    tmp=seq(0.05,0.1,length.out=15)
    # for(dist in c("t")){
    dist='N'
    for(alphai in 1:15){
      alpha=tmp[alphai]
      fileName=paste0("cl_",cl,"_both_",both,"_out_",out,"dist_",dist,"_other_k_exp_alpha1_p05_",
                      alphai,"_beta_p5_p5_N_",N1,"_num_runs_200.Rda",sep="")
      k2=function(t1,t2){k_exp(t1,t2,0.5,alpha)}
      grid=seq( 0, 1, length.out = 100 )
      c2=getCov(grid,k2)
      runMVSim(both,out,contamination_level,N1,N2,c1,c2,grid,num_runs,fileName,dist=dist)
      count=count+1
      print(count/(14*3))
    }
    
    #just guo
  }

}


for(N1 in c(250)){
  both=F
  for(out in 1){
    N2=N1
    print("sim 1")
    #changes in beta
    count=0
    k1=function(t1,t2){k_exp(t1,t2,.5,.05)}
    c1=getCov(grid,k1)
    tmp=seq(0.5,1,length.out = 15)
    # for(dist in c("t")){
    dist='N'
    for(betai in 2:15){
      beta=tmp[betai]
      fileName=paste0("cl_",cl,"_both_",both,"_out_",out,"dist_",dist,"_other_k_exp_alpha1_p05_p05_beta_p",betai,"_p5_N_",N1,"_num_runs_200.Rda",sep="")
      k2=function(t1,t2){k_exp(t1,t2,beta,.05)}
      c2=getCov(grid,k2)
      runMVSim(both,out,contamination_level,N1,N2,c1,c2,grid,num_runs,fileName,dist=dist)
      count=count+1
      print(count/(14*3))
    }
    
    
    
    ##
    # N1=N2=100
    # N1=N2=50
    print("sim 2")
    #changes in alpha
    count=0
    k1=function(t1,t2){k_exp(t1,t2,.5,.05)}
    c1=getCov(grid,k1)
    tmp=seq(0.05,0.1,length.out=15)
    # for(dist in c("t")){
    dist='N'
    for(alphai in 1:15){
      alpha=tmp[alphai]
      fileName=paste0("cl_",cl,"_both_",both,"_out_",out,"dist_",dist,"_other_k_exp_alpha1_p05_",
                      alphai,"_beta_p5_p5_N_",N1,"_num_runs_200.Rda",sep="")
      k2=function(t1,t2){k_exp(t1,t2,0.5,alpha)}
      grid=seq( 0, 1, length.out = 100 )
      c2=getCov(grid,k2)
      runMVSim(both,out,contamination_level,N1,N2,c1,c2,grid,num_runs,fileName,dist=dist)
      count=count+1
      print(count/(14*3))
    }
    
    #just guo
  }
  
}










