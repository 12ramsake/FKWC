################################ Graphs actually in the manuscript


# install.packages("EMMIXskew")
# library(EMMIXskew)
library(sn)
library(fda)
library(roahd)
library(mvtnorm)
library(fdasrvf)

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
  con1=1:floor(contamination_level*(N1+N2))
  if(both_samples)
    con2=N1+1:floor(contamination_level*(N1+N2))
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












#####################


wdd="C:/Users/12RAM/OneDrive/Documents/research/PhD Thesis/Functional Data Covariance Files/Nov_2020_sim_results/Power Files/"
graph_folder="C:/Users/12RAM/OneDrive/Documents/research/PhD Thesis/Functional Data Covariance Files/Nov_2020_sim_results/POWER_GRAPHS_REVISION/"

pink_grad=colorRampPalette(c("#470f32","#821e5c", "#d63696","#ffd1ee"))
cols=c(RColorBrewer::brewer.pal(9, "Blues")[5],"#ffd1ee")
library(RColorBrewer)
cols=c("black",brewer.pal(9, "Blues")[7])


inds=c(1:5,19,21:23,25:26)

make_Plots_2=function(xs=bet,leg=T,b=T){
  #depth plots
  colnames(prop_rejs)[1]="MFHD"
  colnames(prop_rejs)[5]="RPD-L"
  # colnames(prop_rejs)[5]="KSD"
  indz=c(1:5,19,21:23,25:26)
  
  matplot(xs[,indz],prop_rejs[,indz],type='l',
          col=c(rep(cols[1],5),rep(cols[2],6)),
          lty=c(1:5,1:6),
          lw=1,bty="n",
          cex.axis=2,main="",
          ylab="",
          xlab="",yaxt="n",
          xaxt="n")
  if(b)
    mtext(expression(beta),1,line=1.5,cex=2)
  else
    mtext(expression(alpha),1,line=1.5,cex=2)
  mtext("Power",2,line=1,cex=2)
  axis(2,line=-1.5,cex.axis=2)
  axis(1,line=0,cex.axis=2)
  if(leg)
    legend("bottomright",legend=colnames(prop_rejs)[indz],lty=c(1:5,1:6),
           col=c(rep(cols[1],5),rep(cols[2],6)),cex=3,lwd=1,bty = "n")
  # legend("bottomright",legend=c("W","Other"),col=cols,lty=1)
  
}


###############################################################



# We should show gaussian for shape and scale, and t with one
make_plot=function(fn,xs=bet,leg=F,b=T){
  Cairo::CairoPDF(fn,height=10,width=14)
  make_Plots_2(xs,leg,b); 
  dev.off()
}
load(paste(wdd,"dist_N_All_power_rev_k_exp_beta_N_100J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots_2(leg=T); 
make_plot(paste0(graph_folder,"dist_N_All_power_rev_k_exp_beta_N_100J_2.pdf"),leg=T)

load(paste(wdd,"dist_N_All_power_rev_k_exp_alpha_N_100J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots_2(xs=alp,b=F); 
make_plot(paste0(graph_folder,"dist_N_All_power_rev_k_exp_alpha_N_100J_2.pdf"),xs=alp,b=F); 
# make_plot(paste0(graph_folder,"dist_N_All_power_rev_k_exp_beta_N_100J_2.pdf"))
load(paste(wdd,"dist_t_All_power_rev_k_exp_beta_N_100J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots_2(leg=F); 
make_plot(paste0(graph_folder,"dist_t_All_power_rev_k_exp_beta_N_100J_2.pdf"))
load(paste(wdd,"dist_t_All_power_rev_k_exp_alpha_N_100J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots_2(xs=alp,leg=F,b=F); 
make_plot(paste0(graph_folder,"dist_t_All_power_rev_k_exp_alpha_N_100J_2.pdf"))


# Size for each scenario is fine
b_sizes=a_sizes=NULL
tmp=seq(0.05,0.1,length.out=15)
for(dist in c("N","t","sn")){
for(N1 in c(50,100,250)){
    FN=paste0("dist_",dist,"_All_power_rev_k_exp_beta_N_",N1,"J_2_num_runs_200.Rda",sep="")
    load(paste0(wdd,FN))
    b_sizes=cbind(b_sizes,prop_rejs[1,])
    
   
    FN=paste0("dist_",dist,"_All_power_rev_k_exp_alpha_N_",N1,"J_2_num_runs_200.Rda",sep="")
    load(paste0(wdd,FN))
    a_sizes=cbind(a_sizes,prop_rejs[1,])
  }
}

#Size table
rn=c("MFHD","RP$","$MBD$" ,"LTR","RPD^*"   ,rownames(a_sizes)[6:nrow(a_sizes)] )
rownames(a_sizes)=rn
xtable::xtable(a_sizes[inds,],auto=F,label='tab::null',caption='Empirical sizes for $J=2$ for different tests under the infinite dimensional models. The first row indicates the underlying process and the second row indicates the sample size of each group.')



# we should show the eigenvalue table
#Eigenvalue Decay
pr=NULL
for(N1 in c(50,100)){
  dist="N"
  FN=paste0(wdd,"eigen_",dist,"_All_power_rev_N_",N1,"_num_runs_200.Rda",sep="")
  load(FN)
  pr=rbind(pr,prop_rejs[,inds])
}

beepr::beep()
pr=t(pr)
pr2=pr[1:6,7:12]
rownames(pr2)=c("MFHD","RP","MBD" ,"LTR","RPD*","Competing" )
xtable(pr2)


#################################################### Outliers ####################################################

dist="N"
grid=seq( 0, 1, length.out = 100 )
k1=function(t1,t2){k_exp(t1,t2,.5,.05)}
c1=getCov(grid,k1)
k2=k1
c2=getCov(grid,k2)

#Type 1
set.seed(440)
dat=gfd(5,5,Cov1=c1,Cov2=c2,res=100,dist=dist,del=del)
dat=generate_outliers(dat,5,5,contamination_level=0.1,outlier_type=1,both_samples=F)
matplot(matrix(dat$argvals,nrow=ncol(dat$mdata),ncol=nrow(dat$mdata))[,1:5], t(dat$mdata)[,1:5],
        col=c("magenta",rep("black",100)),lty=c(2,rep(1,100)),lwd=c(5,rep(3,100)),main="",type='l',bty="n",xlab="",ylab="")
fn=paste0(graph_folder,"outlier_1.pdf")
Cairo::CairoPDF(fn,height=10,width=12)
matplot(matrix(dat$argvals,nrow=ncol(dat$mdata),ncol=nrow(dat$mdata))[,1:5], t(dat$mdata)[,1:5],
        col=c("magenta",rep("black",100)),lty=c(2,rep(1,100)),lwd=c(5,rep(3,100)),cex.main=4,main="",type='l',bty="n",xlab="",ylab="")
dev.off()

#Type 2
set.seed(440)
dat=gfd(5,5,Cov1=c1,Cov2=c2,res=100,dist=dist,del=del)
dat=generate_outliers(dat,5,5,contamination_level=0.1,outlier_type=2,both_samples=F)
matplot(matrix(dat$argvals,nrow=ncol(dat$mdata),ncol=nrow(dat$mdata))[,1:5], t(dat$mdata)[,1:5],
        col=c("magenta",rep("black",100)),lty=c(2,rep(1,100)),lwd=c(5,rep(3,100)),main="",type='l',bty="n",xlab="",ylab="")

#Type 3
set.seed(440)
dat=gfd(5,5,Cov1=c1,Cov2=c2,res=100,dist=dist,del=del)
dat=generate_outliers(dat,5,5,contamination_level=0.1,outlier_type=3,both_samples=F)
matplot(matrix(dat$argvals,nrow=ncol(dat$mdata),ncol=nrow(dat$mdata))[,1:5], t(dat$mdata)[,1:5],
        col=c("magenta",rep("black",100)),lty=c(2,rep(1,100)),lwd=c(5,rep(3,100)),main="",type='l',bty="n",xlab="",ylab="")
fn=paste0(graph_folder,"outlier_2.pdf")
Cairo::CairoPDF(fn,height=10,width=12)
matplot(matrix(dat$argvals,nrow=ncol(dat$mdata),ncol=nrow(dat$mdata))[,1:5], t(dat$mdata)[,1:5],
        col=c("magenta",rep("black",100)),lty=c(2,rep(1,100)),lwd=c(5,rep(3,100)),main="",type='l',bty="n",xlab="",ylab="")
dev.off()

#Type 4
set.seed(440)
dat=gfd(5,5,Cov1=c1,Cov2=c2,res=100,dist=dist,del=del)
dat=generate_outliers(dat,5,5,contamination_level=0.1,outlier_type=4,both_samples=F)
matplot(matrix(dat$argvals,nrow=ncol(dat$mdata),ncol=nrow(dat$mdata))[,1:5], t(dat$mdata)[,1:5],
        col=c("magenta",rep("black",100)),lty=c(2,rep(1,100)),lwd=c(5,rep(3,100)),main="",type='l',bty="n",xlab="",ylab="")
fn=paste0(graph_folder,"outlier_3.pdf")
Cairo::CairoPDF(fn,height=10,width=12)
matplot(matrix(dat$argvals,nrow=ncol(dat$mdata),ncol=nrow(dat$mdata))[,1:5], t(dat$mdata)[,1:5],
        col=c("magenta",rep("black",100)),lty=c(2,rep(1,100)),lwd=c(5,rep(3,100)),main="",type='l',bty="n",xlab="",ylab="")
dev.off()

############################# BOTH POWER #################################

make_Plots_2=function(xs=bet,leg=T,b=T,main=""){
  #depth plots
  colnames(prop_rejs)[1]="MFHD"
  colnames(prop_rejs)[5]="RPD-L"
  # colnames(prop_rejs)[9]="LTR_PM"
  # colnames(prop_rejs)[5]="KSD"
  
  in2=c(1:5,11,13:15,17:18)
  matplot(xs[,in2],prop_rejs[,in2],type='l',
          col=c(rep(cols[1],5),rep(cols[2],6)),
          lty=c(1:5,1:6),
          lw=1,bty="n",
          cex.axis=2,main=main,cex.main=5,
          ylab="",
          xlab="",yaxt="n",
          xaxt="n")
  if(b)
    mtext(expression(beta),1,line=1.5,cex=2)
  else
    mtext(expression(alpha),1,line=1.5,cex=2)
  mtext("Power",2,line=1,cex=2)
  axis(2,line=-1.5,cex.axis=2)
  axis(1,line=0,cex.axis=2)
  if(leg)
    legend("bottomright",legend=colnames(prop_rejs)[in2],lty=c(1:5,1:6),
           col=c(rep(cols[1],5),rep(cols[2],6)),cex=3,lwd=4,bty = "n")
  # legend("bottomright",legend=c("W","Other"),col=cols,lty=1)
  
}

make_plot=function(fn,xs=bet,leg=F,b=T,main=""){
  Cairo::CairoPDF(fn,height=10,width=14)
  make_Plots_2(xs,leg,b,main=main); 
  dev.off()
}
########################################## 
cl=2.5
# cl=5
both=TRUE
######## Outlier type 1 , 10% observation a 2x upward tilt is added
out=1
N1=100
load(paste(wdd,"cl_",cl,"_both_",both,"_out_",out,"dist_N_All_power_rev_k_exp_beta_N_",N1,"J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots_2()
make_plot(paste0(graph_folder,"cl_2p5_both_",both,"_out_",out,"dist_N_All_power_rev_k_exp_beta_N_",N1,"J_2.pdf"),b=T,main="Drift Outliers"); 
sizes=prop_rejs[1,]

######## Outlier type 2 , 10% observation arandom measurement error
out=2

load(paste(wdd,"cl_",cl,"_both_",both,"_out_",out,"dist_N_All_power_rev_k_exp_beta_N_",N1,"J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots_2()
# make_plot(paste0(graph_folder,"cl_2p5_both_",both,"_out_",out,"dist_N_All_power_rev_k_exp_beta_N_",N1,"J_2.pdf"),b=T); 
######## Outlier type 3 , 10% observation sine waves with random starting point on x axis
out=3

load(paste(wdd,"cl_",cl,"_both_",both,"_out_",out,"dist_N_All_power_rev_k_exp_beta_N_",N1,"J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots_2(xs=bet)
make_plot(paste0(graph_folder,"cl_2p5_both_",both,"_out_",out,"dist_N_All_power_rev_k_exp_beta_N_",N1,"J_2.pdf"),b=T,main="Wavy Outliers"); 
sizes=rbind(sizes,prop_rejs[1,])

######## Outlier type 4 , 10% observation spike outliers
out=4


load(paste(wdd,"cl_",cl,"_both_",both,"_out_",out,"dist_N_All_power_rev_k_exp_beta_N_",N1,"J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots_2()
make_plot(paste0(graph_folder,"cl_2p5_both_",both,"_out_",out,"dist_N_All_power_rev_k_exp_beta_N_",N1,"J_2.pdf"),b=T,main="Spike Outliers");  
sizes=rbind(sizes,prop_rejs[1,])


sizes





############################# FALSE SIZE #################################
########################################## 
cl=2.5
both=FALSE
######## Outlier type 1 , 10% observation a 2x upward tilt is added
out=1
N1=100
load(paste(wdd,"cl_",cl,"_both_",both,"_out_",out,"dist_N_All_power_rev_k_exp_beta_N_",N1,"J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots_2()
make_plot(paste0(graph_folder,"cl_2p5_both_",both,"_out_",out,"dist_N_All_power_rev_k_exp_beta_N_",N1,"J_2.pdf"),b=T,main="Drift Outliers"); 
sizes=prop_rejs[1,]

######## Outlier type 2 , 10% observation arandom measurement error
out=2

load(paste(wdd,"cl_",cl,"_both_",both,"_out_",out,"dist_N_All_power_rev_k_exp_beta_N_",N1,"J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots_2()
# make_plot(paste0(graph_folder,"cl_2p5_both_",both,"_out_",out,"dist_N_All_power_rev_k_exp_beta_N_",N1,"J_2.pdf"),b=T); 
######## Outlier type 3 , 10% observation sine waves with random starting point on x axis
out=3

load(paste(wdd,"cl_",cl,"_both_",both,"_out_",out,"dist_N_All_power_rev_k_exp_beta_N_",N1,"J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots_2(xs=bet)
make_plot(paste0(graph_folder,"cl_2p5_both_",both,"_out_",out,"dist_N_All_power_rev_k_exp_beta_N_",N1,"J_2.pdf"),b=T,main="Wavy Outliers"); 
sizes=rbind(sizes,prop_rejs[1,])

######## Outlier type 4 , 10% observation spike outliers
out=4


load(paste(wdd,"cl_",cl,"_both_",both,"_out_",out,"dist_N_All_power_rev_k_exp_beta_N_",N1,"J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots_2()
make_plot(paste0(graph_folder,"cl_2p5_both_",both,"_out_",out,"dist_N_All_power_rev_k_exp_beta_N_",N1,"J_2.pdf"),b=T,main="Spike Outliers");  
sizes=rbind(sizes,prop_rejs[1,])


xtable::xtable(sizes[,c(1:5,13:15,17:18)])

#################################################### DGS ####################################################



#changes in beta
# in2=c(1:5,13:15,17:18)
tmp=seq(0.5,1,length.out = 15)
tmp[7]
t_pr_2= size=NULL
ratios= c(0.2,0.3,0.4)
N=500
# for(N in c(50,100,250)*2){

for(dist in c("N","t","sn")){
  for(ratio in ratios){
    N1=floor(ratio*N)
    N2=N-N1
    
    FN=paste0("Power Files/dist_",dist,"_All_power_rev_k_exp_beta_N1_",N1,"_N2_",N2,"J_2_num_runs_200.Rda",sep="")
    load(FN)
    size=rbind(size, prop_rejs[1,inds])
    t_pr_2=rbind(t_pr_2,prop_rejs[2,inds])
  }
}
# }
xtable::xtable(t(t_pr_2))
size


t_pr_2= size=NULL
ratios= c(0.2,0.3,0.4)
N=200
N=100
for(dist in c("N","t","sn")){
  for(ratio in ratios){
    N1=floor(ratio*N)
    N2=N-N1
    
    FN=paste0("Power Files/dist_",dist,"_All_power_rev_k_exp_beta_N1_",N1,"_N2_",N2,"J_2_num_runs_200.Rda",sep="")
    load(FN)
    size=rbind(size, prop_rejs[1,inds])
    t_pr_2=rbind(t_pr_2,prop_rejs[2,inds])
  }
}
t_pr_2
size
#L2 br andL2 nv have trouble with DGS

#changes in beta
t_pr_2= size=NULL
ratios= c(0.2,0.3,0.4)*2
N1=100
N1=50
# for(N in c(50,100,250)*2){
for(dist in c("N","t","sn")){
  
  for(ratio in ratios){
    N2=floor(N1*ratio)
    N3=N1
    
    FN=paste0("Power Files/dist_",dist,"_All_power_rev_k_exp_alpha_N1_",N1,"_N2_",N2,"_N3_",N3,"J_3_num_runs_200.Rda",sep="")
    load(FN)
    size=rbind(size, prop_rejs[1,inds])
    t_pr_2=rbind(t_pr_2,prop_rejs[2,inds])
  }
}
# }
t_pr_2
size
#Need to make sure the other tests are working properly


#changes in alpha
t_pr_2= size=NULL
ratios= c(0.2,0.3,0.4)
N=500
# for(N in c(50,100,250)*2){
for(dist in c("N","t","sn")){
  for(ratio in ratios){
    N1=floor(ratio*N)
    N2=N-N1
    
    FN=paste0("Power Files/dist_",dist,"_All_power_rev_k_exp_alpha_N1_",N1,"_N2_",N2,"J_2_num_runs_200.Rda",sep="")
    load(FN)
    size=rbind(size, prop_rejs[1,inds])
    t_pr_2=rbind(t_pr_2,prop_rejs[2,inds])
  }
}
# }
t_pr_2
xtable::xtable(t(t_pr_2))
size


t_pr_2= size=NULL
ratios= c(0.2,0.3,0.4)*2
N=200
N1=100
for(dist in c("N","t","sn")){
  for(ratio in ratios){
    N2=floor(N1*ratio)
    N3=N1
    
    FN=paste0("Power Files/dist_",dist,"_All_power_rev_k_exp_alpha_N1_",N1,"_N2_",N2,"_N3_",N3,"J_3_num_runs_200.Rda",sep="")
    load(FN)
    size=rbind(size, prop_rejs[1,c(1:5,19:26)])
    t_pr_2=rbind(t_pr_2,prop_rejs[2,c(1:5,19:26)])
  }
}
t_pr_2
size

ratios= c(0.2,0.3,0.4)
num_runs=200
library(xtable)
t_pr_2= size=NULL
dist='N'
N=200
for(ratio in ratios){
  N1=floor(ratio*N)
  N2=N-N1
  dist="N"
  FN=paste0(wdd,"eigen_",dist,"_All_power_rev_N1_",N1,"_N2_",N2,"_num_runs_200.Rda",sep="")
  load(FN)
  t_pr=t(prop_rejs)
  colnames(t_pr)=rep(toString(N1),ncol(t_pr))
  t_pr_2=cbind(t_pr_2, t_pr)
}


t_pr_2[inds,]
xtable::xtable(t(t_pr_2[inds,]))
