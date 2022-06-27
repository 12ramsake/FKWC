
#####################


wdd="C:/Users/12RAM/OneDrive/Documents/research/PhD Thesis/Functional Data Covariance Files/Nov_2020_sim_results/Power Files/"
graph_folder="C:/Users/12RAM/OneDrive/Documents/research/PhD Thesis/Functional Data Covariance Files/Nov_2020_sim_results/POWER_GRAPHS_REVISION"

pink_grad=colorRampPalette(c("#470f32","#821e5c", "#d63696","#ffd1ee"))
cols=c(RColorBrewer::brewer.pal(9, "Blues")[5],"#ffd1ee")
library(RColorBrewer)
cols=c("black",brewer.pal(9, "Blues")[7])

make_Plots=function(xs=bet,leg=T){
  
  matplot(xs[,19:26],prop_rejs[,19:26],type='l',
          col=rainbow(8),
          lty=1,
          lw=5,bty="n",
          cex.axis=2,main="",
          ylab="",
          xlab="")
  if(leg)
    legend("bottomright",legend=colnames(prop_rejs)[19:26],lty=1,col=rainbow(8),lwd=3,cex=2)
}


######################################## Regular Scenario  ################################################

######################################### 2 Groups #####################################

load(paste(wdd,"dist_N_All_power_rev_k_exp_beta_N_50J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots()
sizes=prop_rejs[1,19:26]


load(paste(wdd,"dist_N_All_power_rev_k_exp_beta_N_100J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots(); 
sizes=cbind(sizes,prop_rejs[1,19:26])


load(paste(wdd,"dist_N_All_power_rev_k_exp_beta_N_250J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots(); sizes=cbind(sizes,prop_rejs[1,19:26])
sizes=cbind(sizes,prop_rejs[1,19:26])

xtable::xtable(sizes)
# 1% higher size in normal



load(paste(wdd,"dist_N_All_power_rev_k_exp_alpha_N_50J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots(); 
sizes=prop_rejs[1,19:26]



load(paste(wdd,"dist_N_All_power_rev_k_exp_alpha_N_100J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots(); 
sizes=cbind(sizes,prop_rejs[1,19:26])



load(paste(wdd,"dist_N_All_power_rev_k_exp_alpha_N_250J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots(); 
sizes=cbind(sizes,prop_rejs[1,19:26])

#gonna save as 10x12

xtable::xtable(sizes)
# 1% higher size in normal


#T dist
load(paste(wdd,"dist_t_All_power_rev_k_exp_beta_N_50J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots()
sizes=prop_rejs[1,19:26]
# xtable::xtable(sizes)



load(paste(wdd,"dist_t_All_power_rev_k_exp_beta_N_100J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots(); sizes=cbind(sizes,prop_rejs[1,19:26])



load(paste(wdd,"dist_t_All_power_rev_k_exp_beta_N_250J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots(); sizes=cbind(sizes,prop_rejs[1,19:26])


xtable::xtable(sizes)

# same as for normal

load(paste(wdd,"dist_t_All_power_rev_k_exp_alpha_N_50J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots();  sizes=prop_rejs[1,19:26]


load(paste(wdd,"dist_t_All_power_rev_k_exp_alpha_N_100J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots(); sizes=cbind(sizes,prop_rejs[1,19:26])



load(paste(wdd,"dist_t_All_power_rev_k_exp_alpha_N_250J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots(); sizes=cbind(sizes,prop_rejs[1,19:26])
xtable::xtable(sizes)


##SG
load(paste(wdd,"dist_sn_All_power_rev_k_exp_beta_N_50J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots() ;  sizes=prop_rejs[1,19:26]



load(paste(wdd,"dist_sn_All_power_rev_k_exp_beta_N_100J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots(); sizes=cbind(sizes,prop_rejs[1,19:26])

load(paste(wdd,"dist_sn_All_power_rev_k_exp_beta_N_250J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots(); sizes=cbind(sizes,prop_rejs[1,19:26])



xtable::xtable(sizes)




load(paste(wdd,"dist_sn_All_power_rev_k_exp_alpha_N_50J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots();   sizes=prop_rejs[1,19:26]



load(paste(wdd,"dist_sn_All_power_rev_k_exp_alpha_N_100J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots(); sizes=cbind(sizes,prop_rejs[1,19:26])

load(paste(wdd,"dist_sn_All_power_rev_k_exp_alpha_N_250J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots(); sizes=cbind(sizes,prop_rejs[1,19:26])


xtable(sizes)

################################ 3 groups  ########################################



load(paste(wdd,"dist_N_All_power_rev_k_exp_beta_N_50J_3_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots(); sizes=prop_rejs[1,19:26]


load(paste(wdd,"dist_N_All_power_rev_k_exp_beta_N_100J_3_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots(); sizes=cbind(sizes,prop_rejs[1,19:26])


load(paste(wdd,"dist_N_All_power_rev_k_exp_alpha_N_50J_3_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots(); sizes=cbind(sizes,prop_rejs[1,19:26])


load(paste(wdd,"dist_N_All_power_rev_k_exp_alpha_N_100J_3_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots(); sizes=cbind(sizes,prop_rejs[1,19:26])



load(paste(wdd,"dist_t_All_power_rev_k_exp_beta_N_50J_3_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots(); sizes=cbind(sizes,prop_rejs[1,19:26])


load(paste(wdd,"dist_t_All_power_rev_k_exp_beta_N_100J_3_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots(); sizes=cbind(sizes,prop_rejs[1,19:26])


load(paste(wdd,"dist_t_All_power_rev_k_exp_alpha_N_50J_3_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots(); sizes=cbind(sizes,prop_rejs[1,19:26])


load(paste(wdd,"dist_t_All_power_rev_k_exp_alpha_N_100J_3_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots(); sizes=cbind(sizes,prop_rejs[1,19:26])



load(paste(wdd,"dist_sn_All_power_rev_k_exp_beta_N_50J_3_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots(); sizes=cbind(sizes,prop_rejs[1,19:26])


load(paste(wdd,"dist_sn_All_power_rev_k_exp_beta_N_100J_3_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots(); sizes=cbind(sizes,prop_rejs[1,19:26])


load(paste(wdd,"dist_sn_All_power_rev_k_exp_alpha_N_50J_3_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots(); sizes=cbind(sizes,prop_rejs[1,19:26])


load(paste(wdd,"dist_sn_All_power_rev_k_exp_alpha_N_100J_3_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots(); sizes=cbind(sizes,prop_rejs[1,19:26])

sizes
# xtable::xtable(sizes)

######################################## Eigenvalue Scenario  #########################################
num_runs=200
library(xtable)
t_pr_2= size=NULL
for(N1 in c(50,100)){
  dist="N"
  FN=paste0("Power Files/eigen_",dist,"_All_power_rev_N_",N1,"_num_runs_200.Rda",sep="")
  load(FN)
  t_pr=t(prop_rejs)
  colnames(t_pr)=rep(toString(N1),ncol(t_pr))
  t_pr_2=cbind(t_pr_2,t_pr[19:26,])
  
}

xtable(t_pr_2)



######################################## Different Group Sizes  ################################################
########################## Tables ####################################

#changes in beta
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
    size=rbind(size, prop_rejs[1,19:26])
    t_pr_2=rbind(t_pr_2,prop_rejs[2,19:26])
  }
}
# }
t_pr_2
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
    size=rbind(size, prop_rejs[1,19:26])
    t_pr_2=rbind(t_pr_2,prop_rejs[2,19:26])
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
    
    FN=paste0("Power Files/dist_",dist,"_All_power_rev_k_exp_beta_N1_",N1,"_N2_",N2,"_N3_",N3,"J_3_num_runs_200.Rda",sep="")
    load(FN)
    size=rbind(size, prop_rejs[1,19:26])
    t_pr_2=rbind(t_pr_2,prop_rejs[2,19:26])
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
    size=rbind(size, prop_rejs[1,19:26])
    t_pr_2=rbind(t_pr_2,prop_rejs[2,19:26])
  }
}
# }
t_pr_2
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
    size=rbind(size, prop_rejs[1,19:26])
    t_pr_2=rbind(t_pr_2,prop_rejs[2,19:26])
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
  t_pr_2=cbind(t_pr_2,t_pr[19:26,])
}


t_pr_2


#comparable, some dont work well with the smaller sample size



######################################## Outliers  ################################################
########################## In Both  ####################################
both=TRUE
cl=10
dist='N'
######################################### 2 Groups #####################################
######## Outlier type 1 , 10% observation a 2x upward tilt is added
out=1

N1=50
load(paste(wdd,"cl_",cl,"_both_",both,"_out_",out,"dist_N_All_power_rev_k_exp_beta_N_",N1,"J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots()
sizes=prop_rejs[1,]; sizes

N1=100
load(paste(wdd,"cl_",cl,"_both_",both,"_out_",out,"dist_N_All_power_rev_k_exp_beta_N_",N1,"J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots()
sizes=prop_rejs[1,]; sizes

N1=250
load(paste(wdd,"cl_",cl,"_both_",both,"_out_",out,"dist_N_All_power_rev_k_exp_beta_N_",N1,"J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots()
sizes=prop_rejs[1,]; sizes


N1=50
load(paste(wdd,"cl_",cl,"_both_",both,"_out_",out,"dist_N_All_power_rev_k_exp_alpha_N_",N1,"J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots(xs=alp,leg=T)
sizes=prop_rejs[1,]; sizes

N1=100
load(paste(wdd,"cl_",cl,"_both_",both,"_out_",out,"dist_N_All_power_rev_k_exp_alpha_N_",N1,"J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots(xs=alp)
sizes=prop_rejs[1,]; sizes

N1=250
load(paste(wdd,"cl_",cl,"_both_",both,"_out_",out,"dist_N_All_power_rev_k_exp_alpha_N_",N1,"J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots(xs=alp)
sizes=prop_rejs[1,]; sizes


######## Outlier type 2 , 10% observation arandom measurement error
out=2

N1=50
load(paste(wdd,"cl_",cl,"_both_",both,"_out_",out,"dist_N_All_power_rev_k_exp_beta_N_",N1,"J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots()
sizes=prop_rejs[1,]; sizes

N1=100
load(paste(wdd,"cl_",cl,"_both_",both,"_out_",out,"dist_N_All_power_rev_k_exp_beta_N_",N1,"J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots()
sizes=prop_rejs[1,]; sizes

N1=250
load(paste(wdd,"cl_",cl,"_both_",both,"_out_",out,"dist_N_All_power_rev_k_exp_beta_N_",N1,"J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots(xs=alp)
sizes=prop_rejs[1,]; sizes


N1=50
load(paste(wdd,"cl_",cl,"_both_",both,"_out_",out,"dist_N_All_power_rev_k_exp_alpha_N_",N1,"J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots(xs=alp,leg=T)
sizes=prop_rejs[1,]; sizes

N1=100
load(paste(wdd,"cl_",cl,"_both_",both,"_out_",out,"dist_N_All_power_rev_k_exp_alpha_N_",N1,"J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots(xs=alp)
sizes=prop_rejs[1,]; sizes

N1=250
load(paste(wdd,"cl_",cl,"_both_",both,"_out_",out,"dist_N_All_power_rev_k_exp_alpha_N_",N1,"J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots(xs=alp)
sizes=prop_rejs[1,]; sizes

######## Outlier type 3 , 10% observation sine waves with random starting point on x axis
out=3

N1=50
load(paste(wdd,"cl_",cl,"_both_",both,"_out_",out,"dist_N_All_power_rev_k_exp_beta_N_",N1,"J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots()
sizes=prop_rejs[1,]; sizes

N1=100
load(paste(wdd,"cl_",cl,"_both_",both,"_out_",out,"dist_N_All_power_rev_k_exp_beta_N_",N1,"J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots()
sizes=prop_rejs[1,]; sizes

N1=250
load(paste(wdd,"cl_",cl,"_both_",both,"_out_",out,"dist_N_All_power_rev_k_exp_beta_N_",N1,"J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots(xs=alp)
sizes=prop_rejs[1,]; sizes


N1=50
load(paste(wdd,"cl_",cl,"_both_",both,"_out_",out,"dist_N_All_power_rev_k_exp_alpha_N_",N1,"J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots(xs=alp,leg=T)
sizes=prop_rejs[1,]; sizes

N1=100
load(paste(wdd,"cl_",cl,"_both_",both,"_out_",out,"dist_N_All_power_rev_k_exp_alpha_N_",N1,"J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots(xs=alp)
sizes=prop_rejs[1,]; sizes

N1=250
load(paste(wdd,"cl_",cl,"_both_",both,"_out_",out,"dist_N_All_power_rev_k_exp_alpha_N_",N1,"J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots(xs=alp)
sizes=prop_rejs[1,]; sizes

#Derivs more comparable here


######## Outlier type 4 , 10% observation spike outliers
out=4

N1=50
load(paste(wdd,"cl_",cl,"_both_",both,"_out_",out,"dist_N_All_power_rev_k_exp_beta_N_",N1,"J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots()
sizes=prop_rejs[1,]; sizes

N1=100
load(paste(wdd,"cl_",cl,"_both_",both,"_out_",out,"dist_N_All_power_rev_k_exp_beta_N_",N1,"J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots()
sizes=prop_rejs[1,]; sizes

N1=250
load(paste(wdd,"cl_",cl,"_both_",both,"_out_",out,"dist_N_All_power_rev_k_exp_beta_N_",N1,"J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots(xs=alp)
sizes=prop_rejs[1,]; sizes


N1=50
load(paste(wdd,"cl_",cl,"_both_",both,"_out_",out,"dist_N_All_power_rev_k_exp_alpha_N_",N1,"J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots(xs=alp,leg=T)
sizes=prop_rejs[1,]; sizes

N1=100
load(paste(wdd,"cl_",cl,"_both_",both,"_out_",out,"dist_N_All_power_rev_k_exp_alpha_N_",N1,"J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots(xs=alp)
sizes=prop_rejs[1,]; sizes

N1=250
load(paste(wdd,"cl_",cl,"_both_",both,"_out_",out,"dist_N_All_power_rev_k_exp_alpha_N_",N1,"J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots(xs=alp)
sizes=prop_rejs[1,]; sizes



#EigenValues
num_runs=200

#Eigenvalue Decay
prop_rejs=NULL
N1=100
for(out in 1:4){
  dist="N"
  load(paste0("Power Files/cl_",cl,"_both_",both,"_out_",out,"eigen_",dist,"_All_power_rev_N_",N1,"_num_runs_200.Rda",sep=""))
  print(paste0("out ",out))
  # print(t(prop_rejs)[c(1:5),]-t(prop_rejs)[c(1:5),])
  print(t(prop_rejs)[19:26,])
  # print(xtable(t(prop_rejs)[c(1:4,11:14),]))
  
}
# THe spike outliers destroy LTR with the derivative



prop_rejs=NULL
N1=50
for(out in 1:4){
  dist="N"
  load(paste0("Power Files/cl_",cl,"_both_",both,"_out_",out,"eigen_",dist,"_All_power_rev_N_",N1,"_num_runs_200.Rda",sep=""))
  print(paste0("out ",out))
  print(t(prop_rejs)[19:26,])
  # print(t(prop_rejs)[c(1:4,11:14),])
  # print(xtable(t(prop_rejs)[c(1:4,11:14),]))
  
}
# THe sp




########################## In One sample  ####################################
both=FALSE
cl=10
#can adjust cl 
cl=1
dist='N'
######################################### 2 Groups #####################################
######## Outlier type 1 , 10% observation a 2x upward tilt is added
out=1

N1=50
load(paste(wdd,"cl_",cl,"_both_",both,"_out_",out,"dist_N_All_power_rev_k_exp_beta_N_",N1,"J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots()
sizes=prop_rejs[1,]; sizes

N1=100
load(paste(wdd,"cl_",cl,"_both_",both,"_out_",out,"dist_N_All_power_rev_k_exp_beta_N_",N1,"J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots()
sizes=prop_rejs[1,]; sizes

N1=250
load(paste(wdd,"cl_",cl,"_both_",both,"_out_",out,"dist_N_All_power_rev_k_exp_beta_N_",N1,"J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots(xs=alp)
sizes=prop_rejs[1,]; sizes


N1=50
load(paste(wdd,"cl_",cl,"_both_",both,"_out_",out,"dist_N_All_power_rev_k_exp_alpha_N_",N1,"J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots(xs=alp,leg=T)
sizes=prop_rejs[1,]; sizes

N1=100
load(paste(wdd,"cl_",cl,"_both_",both,"_out_",out,"dist_N_All_power_rev_k_exp_alpha_N_",N1,"J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots(xs=alp)
sizes=prop_rejs[1,]; sizes

N1=250
load(paste(wdd,"cl_",cl,"_both_",both,"_out_",out,"dist_N_All_power_rev_k_exp_alpha_N_",N1,"J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots(xs=alp)
sizes=prop_rejs[1,]; sizes


######## Outlier type 2 , 10% observation arandom measurement error
out=2

N1=50
load(paste(wdd,"cl_",cl,"_both_",both,"_out_",out,"dist_N_All_power_rev_k_exp_beta_N_",N1,"J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots()
sizes=prop_rejs[1,]; sizes

N1=100
load(paste(wdd,"cl_",cl,"_both_",both,"_out_",out,"dist_N_All_power_rev_k_exp_beta_N_",N1,"J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots()
sizes=prop_rejs[1,]; sizes

N1=250
load(paste(wdd,"cl_",cl,"_both_",both,"_out_",out,"dist_N_All_power_rev_k_exp_beta_N_",N1,"J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots(xs=alp)
sizes=prop_rejs[1,]; sizes


N1=50
load(paste(wdd,"cl_",cl,"_both_",both,"_out_",out,"dist_N_All_power_rev_k_exp_alpha_N_",N1,"J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots(xs=alp,leg=T)
sizes=prop_rejs[1,]; sizes

N1=100
load(paste(wdd,"cl_",cl,"_both_",both,"_out_",out,"dist_N_All_power_rev_k_exp_alpha_N_",N1,"J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots(xs=alp)
sizes=prop_rejs[1,]; sizes

N1=250
load(paste(wdd,"cl_",cl,"_both_",both,"_out_",out,"dist_N_All_power_rev_k_exp_alpha_N_",N1,"J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots(xs=alp)
sizes=prop_rejs[1,]; sizes

######## Outlier type 3 , 10% observation sine waves with random starting point on x axis
out=3

N1=50
load(paste(wdd,"cl_",cl,"_both_",both,"_out_",out,"dist_N_All_power_rev_k_exp_beta_N_",N1,"J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots()
sizes=prop_rejs[1,]; sizes

N1=100
load(paste(wdd,"cl_",cl,"_both_",both,"_out_",out,"dist_N_All_power_rev_k_exp_beta_N_",N1,"J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots()
sizes=prop_rejs[1,]; sizes

N1=250
load(paste(wdd,"cl_",cl,"_both_",both,"_out_",out,"dist_N_All_power_rev_k_exp_beta_N_",N1,"J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots(xs=alp)
sizes=prop_rejs[1,]; sizes


N1=50
load(paste(wdd,"cl_",cl,"_both_",both,"_out_",out,"dist_N_All_power_rev_k_exp_alpha_N_",N1,"J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots(xs=alp,leg=T)
sizes=prop_rejs[1,]; sizes

N1=100
load(paste(wdd,"cl_",cl,"_both_",both,"_out_",out,"dist_N_All_power_rev_k_exp_alpha_N_",N1,"J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots(xs=alp)
sizes=prop_rejs[1,]; sizes

N1=250
load(paste(wdd,"cl_",cl,"_both_",both,"_out_",out,"dist_N_All_power_rev_k_exp_alpha_N_",N1,"J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots(xs=alp)
sizes=prop_rejs[1,]; sizes

#Derivs more comparable here


######## Outlier type 4 , 10% observation spike outliers
out=4

N1=50
load(paste(wdd,"cl_",cl,"_both_",both,"_out_",out,"dist_N_All_power_rev_k_exp_beta_N_",N1,"J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots()
sizes=prop_rejs[1,]; sizes

N1=100
load(paste(wdd,"cl_",cl,"_both_",both,"_out_",out,"dist_N_All_power_rev_k_exp_beta_N_",N1,"J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots()
sizes=prop_rejs[1,]; sizes

N1=250
load(paste(wdd,"cl_",cl,"_both_",both,"_out_",out,"dist_N_All_power_rev_k_exp_beta_N_",N1,"J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots(xs=alp)
sizes=prop_rejs[1,]; sizes


N1=50
load(paste(wdd,"cl_",cl,"_both_",both,"_out_",out,"dist_N_All_power_rev_k_exp_alpha_N_",N1,"J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots(xs=alp,leg=T)
sizes=prop_rejs[1,]; sizes

N1=100
load(paste(wdd,"cl_",cl,"_both_",both,"_out_",out,"dist_N_All_power_rev_k_exp_alpha_N_",N1,"J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots(xs=alp)
sizes=prop_rejs[1,]; sizes

N1=250
load(paste(wdd,"cl_",cl,"_both_",both,"_out_",out,"dist_N_All_power_rev_k_exp_alpha_N_",N1,"J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots(xs=alp)
sizes=prop_rejs[1,]; sizes





# library(xtable)
# sizesT=t(sizes)
# xtable(sizes[c(1:4,11:14),])
# xtable(sizes)

# both+FALSE

prop_rejs=NULL
N1=100
for(out in 1:4){
  dist="N"
  load(paste0("Power Files/cl_",cl,"_both_",both,"_out_",out,"eigen_",dist,"_All_power_rev_N_",N1,"_num_runs_200.Rda",sep=""))
  print(paste0("out ",out))
  # print(t(prop_rejs)[c(1:4),]-t(prop_rejs)[c(11:14),])
  print(t(prop_rejs)[19:26,])
  # print(xtable(t(prop_rejs)[c(1:4,11:14),]))
  
}
# THe spike outliers destroy LTR with the derivative



prop_rejs=NULL
N1=50
for(out in 1:4){
  dist="N"
  load(paste0("Power Files/cl_",cl,"_both_",both,"_out_",out,"eigen_",dist,"_All_power_rev_N_",N1,"_num_runs_200.Rda",sep=""))
  print(paste0("out ",out))
  print(t(prop_rejs)[19:26,])
  # print(t(prop_rejs)[c(1:4,11:14),])
  # print(xtable(t(prop_rejs)[c(1:4,11:14),]))
  
}
# THe sp


########### ANALYZING ADDITIONAL OUTLIERS




make_Plots=function(xs=bet,leg=T){
  
  matplot(xs[,c(11:18)],prop_rejs[,c(11:18)],type='l',
          col=rainbow(8),
          lty=1,
          lw=5,bty="n",
          cex.axis=2,main="",
          ylab="",
          xlab="")
  if(leg)
    legend("bottomright",legend=colnames(prop_rejs)[c(11:18)],lty=1,col=rainbow(8),lwd=3,cex=2)
}



########################## In One sample  ####################################
both=FALSE
cl=10
#can adjust cl 
cl=1
dist='N'
######################################### 2 Groups #####################################
######## Outlier type 1 , 10% observation a 2x upward tilt is added
out=1
N1=100
load(paste(wdd,"cl_",cl,"_both_",both,"_out_",out,"dist_N_All_power_rev_k_exp_beta_N_",N1,"J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots()
sizes=prop_rejs[1,c(1:5,11:18)]; sizes


load(paste(wdd,"cl_",cl,"_both_",both,"_out_",out,"dist_N_All_power_rev_k_exp_alpha_N_",N1,"J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots(xs=alp)
sizes=prop_rejs[1,]; sizes


######## Outlier type 2 , 10% observation arandom measurement error
out=2

load(paste(wdd,"cl_",cl,"_both_",both,"_out_",out,"dist_N_All_power_rev_k_exp_beta_N_",N1,"J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots()
sizes=prop_rejs[1,]; sizes


load(paste(wdd,"cl_",cl,"_both_",both,"_out_",out,"dist_N_All_power_rev_k_exp_alpha_N_",N1,"J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots(xs=alp)
sizes=prop_rejs[1,]; sizes

######## Outlier type 3 , 10% observation sine waves with random starting point on x axis
out=3



load(paste(wdd,"cl_",cl,"_both_",both,"_out_",out,"dist_N_All_power_rev_k_exp_beta_N_",N1,"J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots()
sizes=prop_rejs[1,]; sizes


load(paste(wdd,"cl_",cl,"_both_",both,"_out_",out,"dist_N_All_power_rev_k_exp_alpha_N_",N1,"J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots(xs=alp)
sizes=prop_rejs[1,]; sizes



######## Outlier type 4 , 10% observation spike outliers
out=4


load(paste(wdd,"cl_",cl,"_both_",both,"_out_",out,"dist_N_All_power_rev_k_exp_beta_N_",N1,"J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots()
sizes=prop_rejs[1,]; sizes


load(paste(wdd,"cl_",cl,"_both_",both,"_out_",out,"dist_N_All_power_rev_k_exp_alpha_N_",N1,"J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots(xs=alp)
sizes=prop_rejs[1,]; sizes

########################################## 
cl=2.5
######## Outlier type 1 , 10% observation a 2x upward tilt is added
out=1
load(paste(wdd,"cl_",cl,"_both_",both,"_out_",out,"dist_N_All_power_rev_k_exp_beta_N_",N1,"J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots()
sizes=prop_rejs[1,]; sizes


load(paste(wdd,"cl_",cl,"_both_",both,"_out_",out,"dist_N_All_power_rev_k_exp_alpha_N_",N1,"J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots(xs=alp)
sizes=prop_rejs[1,]; sizes


######## Outlier type 2 , 10% observation arandom measurement error
out=2

load(paste(wdd,"cl_",cl,"_both_",both,"_out_",out,"dist_N_All_power_rev_k_exp_beta_N_",N1,"J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots()
sizes=prop_rejs[1,]; sizes


load(paste(wdd,"cl_",cl,"_both_",both,"_out_",out,"dist_N_All_power_rev_k_exp_alpha_N_",N1,"J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots(xs=alp)
sizes=prop_rejs[1,]; sizes

######## Outlier type 3 , 10% observation sine waves with random starting point on x axis
out=3



load(paste(wdd,"cl_",cl,"_both_",both,"_out_",out,"dist_N_All_power_rev_k_exp_beta_N_",N1,"J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots()
sizes=prop_rejs[1,]; sizes


load(paste(wdd,"cl_",cl,"_both_",both,"_out_",out,"dist_N_All_power_rev_k_exp_alpha_N_",N1,"J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots(xs=alp)
sizes=prop_rejs[1,]; sizes



######## Outlier type 4 , 10% observation spike outliers
out=4


load(paste(wdd,"cl_",cl,"_both_",both,"_out_",out,"dist_N_All_power_rev_k_exp_beta_N_",N1,"J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots()
sizes=prop_rejs[1,]; sizes


load(paste(wdd,"cl_",cl,"_both_",both,"_out_",out,"dist_N_All_power_rev_k_exp_alpha_N_",N1,"J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots(xs=alp)
sizes=prop_rejs[1,]; sizes















