wdd="C:/Users/12RAM/OneDrive/Documents/research/PhD Thesis/Functional Data Covariance Files/Nov_2020_sim_results/Power Files/"
graph_folder="C:/Users/12RAM/OneDrive/Documents/research/PhD Thesis/Functional Data Covariance Files/Nov_2020_sim_results/POWER_GRAPHS_REVISION"
colfunc <- colorRampPalette(c("orange","red","magenta", "blue"))
#J=2
#BETA
sizes=NULL
cols=c("darkgreen","black")
pink_grad=colorRampPalette(c("#470f32","#821e5c", "#d63696","#ffd1ee"))
cols=c(RColorBrewer::brewer.pal(9, "Blues")[5],"#ffd1ee")


#for saving
make_Plots=function(xs=bet,leg=T){
  #depth plots
  colnames(prop_rejs)[1]="MFHD"
  # colnames(prop_rejs)[8]="LTR_PM"
  # colnames(prop_rejs)[5]="KSD"
  num_d=4
  #Derivative vs No Derivative
  matplot(xs[,c(1:num_d,11:14)],prop_rejs[,c(1:num_d,11:14)],type='l',
          main="",
          col=rep(cols,each=4),lty=rep(1:num_d,2)+1,lw=3,bty="n",
          cex.axis=2,ylab="",xlab="")
  if(leg){
    legend(.7,.7,legend=colnames(prop_rejs)[1:num_d],col=1,lty=1:num_d+1,cex=3,lwd=3,bty = "n")
    legend(.7,.3,legend=c("Derivative","No Derivative"),col=cols,lty=1,cex=3,lwd=4,bty = "n")
  }
  
  
  
  num_d=5
  #Derivative vs No Derivative
  matplot(xs[,1:num_d],prop_rejs[,1:num_d],type='l',
          main="",
          col=rainbow(num_d),lty=1,lw=3,bty="n",
          cex.axis=2,ylab="",xlab="")
  if(leg){
    legend(.7,.7,legend=colnames(prop_rejs)[1:num_d],col=rainbow(num_d),lty=1,cex=3,lwd=3,bty = "n")
  }
  
  
  
  
  
  #THIS IS THE PM
  matplot(xs[,1:10],prop_rejs[,c(1:10)],
          main="",
          type='l',col=rep(cols,each=5),lty=rep(1:5,2),
          lw=5,bty="n",cex.axis=1.5,ylab="")
  legend("topright",legend=colnames(prop_rejs)[1:5],col=1,lty=1:5,cex=1.5)
  legend(xs[5,1],0.7,legend=c(expression(italic(W)),expression(italic(M))),col=cols,lty=1,cex=1.5)
  
  #This is Deriv vs Competing
  matplot(xs[,c(6:10,19:26)],prop_rejs[,c(6:10,19:26)],type='l',
          col=c(rep(cols[1],5),rep(cols[2],8)),
          lty=c(1:5,1:8),
          lw=5,bty="n",
          cex.axis=2,main="",
          ylab="",
          xlab="")
  if(leg)
    legend("bottomright",legend=colnames(prop_rejs)[c(1:5,19:26)],lty=c(1:5,1:8),col=c(rep(cols[1],5),rep(cols[2],8)),lwd=3,cex=2)
  # legend("bottomright",legend=c("W","Other"),col=cols,lty=1)
  
}

#deriv vs no derivative
make_Plots_1=function(xs=bet,leg=T,b=T){
  #depth plots
  colnames(prop_rejs)[1]="MFHD"
  # colnames(prop_rejs)[8]="LTR_PM"
  # colnames(prop_rejs)[5]="KSD"
  matplot(xs[,c(1:4,11:14)],prop_rejs[,c(1:4,11:14)],type='l',
          main="",
          col=rep(cols,each=4),lty=rep(1:4,2),lw=5,bty="n",yaxt="n",
          xaxt="n",ylab="",xlab="")
  if(b)
    mtext(expression(beta),1,line=1.5,cex=2)
  else
    mtext(expression(alpha),1,line=1.5,cex=2)
  mtext("Power",2,line=1,cex=2)
  axis(2,line=-1.5,cex.axis=2)
  axis(1,line=0,cex.axis=2)
  
  if(leg){
    legend(.7,.7,legend=colnames(prop_rejs)[1:4],col=1,lty=1:5,cex=3,lwd=4,bty = "n")
    legend(.7,.3,legend=c("Derivative","No Derivative"),col=cols,lty=1,cex=3,lwd=4,bty = "n")
  }
  
  
}


make_Plots_2=function(xs=bet,leg=T,b=T){
  #depth plots
  colnames(prop_rejs)[1]="MFHD"
  # colnames(prop_rejs)[9]="LTR_PM"
  # colnames(prop_rejs)[5]="KSD"
  
  
  matplot(xs[,c(1:5,19,23:26)],prop_rejs[,c(1:5,19,23:26)],type='l',
          col=c(rep(cols[1],5),rep(cols[2],5)),
          lty=c(1:5,1:5),
          lw=5,bty="n",
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
    legend("bottomright",legend=colnames(prop_rejs)[c(1:5,19,23:26)],lty=c(1:5,1:5),
           col=c(rep(cols[1],4),rep(cols[2],5)),cex=3,lwd=4,bty = "n")
  # legend("bottomright",legend=c("W","Other"),col=cols,lty=1)
  
}


# make_Plots=function(){
#         #depth plots
#         matplot(prop_rejs[,c(1:5,11:15)],type='l',
#                 main="Deriv vs No Deriv",
#                 col=rep(cols,each=5),lty=rep(1:5,2),lw=5,bty="n",cex.axis=1.5,ylab="")
#         legend("topright",legend=colnames(prop_rejs)[1:5],col=1,lty=1:5)
#         legend("bottomright",legend=c("Deriv","No Deriv"),col=cols,lty=1)
#         
#         
#         matplot(prop_rejs[,c(1:10)],
#                 main="KW vs PM",
#                 type='l',col=rep(cols,each=5),lty=rep(1:5,2),
#                 lw=5,bty="n",cex.axis=1.5,ylab="")
#         legend("topright",legend=colnames(prop_rejs)[1:5],col=1,lty=1:5)
#         legend("bottomright",legend=c("W","M"),col=cols,lty=1)
#         
#         
#         matplot(prop_rejs[,c(1:5,21:28)],type='l',col=c(rep(cols[1],5),rep(cols[2],8)),
#                 lty=c(1:5,1:8),lw=5,bty="n",cex.axis=1.5,main="KW vs Others",ylab="")
#         legend("topright",legend=colnames(prop_rejs)[c(1:5,21:28)],lty=c(1:5,1:8),col=c(rep(cols[1],5),rep(cols[2],8)))
#         # legend("bottomright",legend=c("W","Other"),col=cols,lty=1)
#         
# }
library(RColorBrewer)
cols=c("black",brewer.pal(9, "Blues")[7])

# load(paste(wdd,"PM_dist_N_SIM_TD_All_power_rev_k_exp_alpha_N_50J_2_num_runs_200.rda",sep="")); dim(prop_rejs)
# new_depth_def=prop_rejs[,15,]
load(paste(wdd,"dist_N_All_power_rev_k_exp_beta_N_50J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
# prop_rejs[,'RPD']=new_depth_def[1,]
# prop_rejs[,'LTR']=new_depth_def[2,]
make_Plots()
sizes=cbind(sizes,prop_rejs[1,])

pdf(file="DvsNDdist_N_beta_N_50J_2_rev.pdf",width=12,height=10)
make_Plots_1();
dev.off()
pdf(file="KvOdist_N_beta_N_50J_2_rev.pdf",width=12,height=10)
make_Plots_2();
dev.off()


load(paste(wdd,"dist_N_All_power_rev_k_exp_beta_N_100J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots(); 
sizes=cbind(sizes,prop_rejs[1,])
xtable::xtable(sizes)
pdf(file="DvsNDdist_N_beta_N_100J_2_rev.pdf",width=12,height=10)
make_Plots_1();
dev.off()
pdf(file="KvOdist_N_beta_N_100J_2_rev.pdf",width=12,height=10)
make_Plots_2();
dev.off()

load(paste(wdd,"dist_N_All_power_rev_k_exp_beta_N_250J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots(); sizes=cbind(sizes,prop_rejs[1,])

pdf(file="DvsNDdist_N_beta_N_250J_2_rev.pdf",width=12,height=10)
make_Plots_1();
dev.off()
pdf(file="KvOdist_N_beta_N_250J_2_rev.pdf",width=12,height=10)
make_Plots_2();
dev.off()





# DvsNDdist_N_beta_N_250J_2
# KvOdist_N_beta_N_250J_2 
#10x12

load(paste(wdd,"dist_N_All_power_rev_k_exp_alpha_N_50J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots(); 
sizes=cbind(sizes,prop_rejs[1,])
pdf(file="DvsNDdist_N_alpha_N_50J_2_rev.pdf",width=12,height=10)
make_Plots_1(xs=alp,leg=F,b=F);
dev.off()
pdf(file="KvOdist_N_alpha_N_50J_2_rev.pdf",width=12,height=10)
make_Plots_2(xs=alp,leg=F,b=F);
dev.off()


load(paste(wdd,"dist_N_All_power_rev_k_exp_alpha_N_100J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots(); 
sizes=cbind(sizes,prop_rejs[1,])
pdf(file="DvsNDdist_N_alpha_N_100J_2_rev.pdf",width=12,height=10)
make_Plots_1(xs=alp,leg=F,b=F);
dev.off()
pdf(file="KvOdist_N_alpha_N_100J_2_rev.pdf",width=12,height=10)
make_Plots_2(xs=alp,leg=F,b=F);
dev.off()


load(paste(wdd,"dist_N_All_power_rev_k_exp_alpha_N_250J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots(); 

sizes=cbind(sizes,prop_rejs[1,])
pdf(file="DvsNDdist_N_alpha_N_250J_2_rev.pdf",width=12,height=10)
make_Plots_1(xs=alp,leg=F,b=F);
dev.off()
pdf(file="KvOdist_N_alpha_N_250J_2_rev.pdf",width=12,height=10)
make_Plots_2(xs=alp,leg=F,b=F);
dev.off()

#gonna save as 10x12
# KvOdist_N_alpha_N_250J_2 
# DvsNDdist_N_alpha_N_250J_2


load(paste(wdd,"dist_N_All_power_rev_k_exp_beta_N_50J_3_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots(); #sizes=cbind(sizes,prop_rejs[1,])
pdf(file="DvsNDdist_N_beta_N_50J_3_rev.pdf",width=12,height=10)
make_Plots_1();
dev.off()
pdf(file="KvOdist_N_beta_N_50J_3_rev.pdf",width=12,height=10)
make_Plots_2();
dev.off()

load(paste(wdd,"dist_N_All_power_rev_k_exp_beta_N_100J_3_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots(); #sizes=cbind(sizes,prop_rejs[1,])
pdf(file="DvsNDdist_N_beta_N_100J_3_rev.pdf",width=12,height=10)
make_Plots_1();
dev.off()
pdf(file="KvOdist_N_beta_N_100J_3_rev.pdf",width=12,height=10)
make_Plots_2();
dev.off()



load(paste(wdd,"dist_N_All_power_rev_k_exp_alpha_N_50J_3_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots();
pdf(file="DvsNDdist_N_alpha_N_50J_3_rev.pdf",width=12,height=10)
make_Plots_1(xs=alp,leg=F,b=F);
dev.off()
pdf(file="KvOdist_N_alpha_N_50J_3_rev.pdf",width=12,height=10)
make_Plots_2(xs=alp,leg=F,b=F);
dev.off()

load(paste(wdd,"dist_N_All_power_rev_k_exp_alpha_N_100J_3_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots(); 
pdf(file="DvsNDdist_N_alpha_N_100J_3_rev.pdf",width=12,height=10)
make_Plots_1(xs=alp,leg=F,b=F);
dev.off()
pdf(file="KvOdist_N_alpha_N_100J_3_rev.pdf",width=12,height=10)
make_Plots_2(xs=alp,leg=F,b=F);
dev.off()






#T dist
load(paste(wdd,"dist_t_All_power_rev_k_exp_beta_N_50J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots()
sizes=cbind(sizes,prop_rejs[1,])
xtable::xtable(sizes)
pdf(file="DvsNDdist_t_beta_N_50J_2_rev.pdf",width=12,height=10)
make_Plots_1();
dev.off()
pdf(file="KvOdist_t_beta_N_50J_2_rev.pdf",width=12,height=10)
make_Plots_2();
dev.off()


load(paste(wdd,"dist_t_All_power_rev_k_exp_beta_N_100J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots(); sizes=cbind(sizes,prop_rejs[1,])
xtable::xtable(sizes)
pdf(file="DvsNDdist_t_beta_N_100J_2_rev.pdf",width=12,height=10)
make_Plots_1();
dev.off()
pdf(file="KvOdist_t_beta_N_100J_2_rev.pdf",width=12,height=10)
make_Plots_2();
dev.off()

load(paste(wdd,"dist_t_All_power_rev_k_exp_beta_N_250J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots(); sizes=cbind(sizes,prop_rejs[1,])

pdf(file="DvsNDdist_t_beta_N_250J_2_rev.pdf",width=12,height=10)
make_Plots_1(leg=F);
dev.off()
pdf(file="KvOdist_t_beta_N_250J_2_rev.pdf",width=12,height=10)
make_Plots_2(leg=F);
dev.off()





# DvsNDdist_t_beta_N_250J_2
# KvOdist_t_beta_N_250J_2 
#10x12

load(paste(wdd,"dist_t_All_power_rev_k_exp_alpha_N_50J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots(); 
pdf(file="DvsNDdist_t_alpha_N_50J_2_rev.pdf",width=12,height=10)
make_Plots_1(xs=alp,leg=F,b=F);
dev.off()
pdf(file="KvOdist_t_alpha_N_50J_2_rev.pdf",width=12,height=10)
make_Plots_2(xs=alp,leg=F,b=F);
dev.off()


load(paste(wdd,"dist_t_All_power_rev_k_exp_alpha_N_100J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots();
pdf(file="DvsNDdist_t_alpha_N_100J_2_rev.pdf",width=12,height=10)
make_Plots_1(xs=alp,leg=F,b=F);
dev.off()
pdf(file="KvOdist_t_alpha_N_100J_2_rev.pdf",width=12,height=10)
make_Plots_2(xs=alp,leg=F,b=F);
dev.off()


load(paste(wdd,"dist_t_All_power_rev_k_exp_alpha_N_250J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots(); 


pdf(file="DvsNDdist_t_alpha_N_250J_2_rev.pdf",width=12,height=10)
make_Plots_1(xs=alp,leg=F,b=F);
dev.off()
pdf(file="KvOdist_t_alpha_N_250J_2_rev.pdf",width=12,height=10)
make_Plots_2(xs=alp,leg=F,b=F);
dev.off()

#gonna save as 10x12
# KvOdist_t_alpha_N_250J_2 
# DvsNDdist_t_alpha_N_250J_2


load(paste(wdd,"dist_t_All_power_rev_k_exp_beta_N_50J_3_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots(); #sizes=cbind(sizes,prop_rejs[1,])
pdf(file="DvsNDdist_t_beta_N_50J_3_rev.pdf",width=12,height=10)
make_Plots_1();
dev.off()
pdf(file="KvOdist_t_beta_N_50J_3_rev.pdf",width=12,height=10)
make_Plots_2();
dev.off()

load(paste(wdd,"dist_t_All_power_rev_k_exp_beta_N_100J_3_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots(); #sizes=cbind(sizes,prop_rejs[1,])
pdf(file="DvsNDdist_t_beta_N_100J_3_rev.pdf",width=12,height=10)
make_Plots_1();
dev.off()
pdf(file="KvOdist_t_beta_N_100J_3_rev.pdf",width=12,height=10)
make_Plots_2();
dev.off()



load(paste(wdd,"dist_t_All_power_rev_k_exp_alpha_N_50J_3_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots(); 
pdf(file="DvsNDdist_t_alpha_N_50J_3_rev.pdf",width=12,height=10)
make_Plots_1(xs=alp,leg=F,b=F);
dev.off()
pdf(file="KvOdist_t_alpha_N_50J_3_rev.pdf",width=12,height=10)
make_Plots_2(xs=alp,leg=F,b=F);
dev.off()

load(paste(wdd,"dist_t_All_power_rev_k_exp_alpha_N_100J_3_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots(); 
pdf(file="DvsNDdist_t_alpha_N_100J_3_rev.pdf",width=12,height=10)
make_Plots_1(xs=alp,leg=F,b=F);
dev.off()
pdf(file="KvOdist_t_alpha_N_100J_3_rev.pdf",width=12,height=10)
make_Plots_2(xs=alp,leg=F,b=F);
dev.off()

##SG
load(paste(wdd,"dist_sn_All_power_rev_k_exp_beta_N_50J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots()
sizes=cbind(sizes,prop_rejs[1,])
pdf(file="DvsNDdist_sg_beta_N_50J_2_rev.pdf",width=12,height=10)
make_Plots_1();
dev.off()
pdf(file="KvOdist_sg_beta_N_50J_2_rev.pdf",width=12,height=10)
make_Plots_2();
dev.off()


load(paste(wdd,"dist_sn_All_power_rev_k_exp_beta_N_100J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots(); sizes=cbind(sizes,prop_rejs[1,])

pdf(file="DvsNDdist_sg_beta_N_100J_2_rev.pdf",width=12,height=10)
make_Plots_1();
dev.off()
pdf(file="KvOdist_sg_beta_N_100J_2_rev.pdf",width=12,height=10)
make_Plots_2();
dev.off()

load(paste(wdd,"dist_sn_All_power_rev_k_exp_beta_N_250J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots(); sizes=cbind(sizes,prop_rejs[1,])

pdf(file="DvsNDdist_sg_beta_N_250J_2_rev.pdf",width=12,height=10)
make_Plots_1(leg=F);
dev.off()
pdf(file="KvOdist_sg_beta_N_250J_2_rev.pdf",width=12,height=10)
make_Plots_2(leg=F);
dev.off()





# DvsNDdist_sg_beta_N_250J_2
# KvOdist_sg_beta_N_250J_2 
#10x12

load(paste(wdd,"dist_sn_All_power_rev_k_exp_alpha_N_50J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots(); 
pdf(file="DvsNDdist_sg_alpha_N_50J_2_rev.pdf",width=12,height=10)
make_Plots_1(xs=alp,leg=F,b=F);
dev.off()
pdf(file="KvOdist_sg_alpha_N_50J_2_rev.pdf",width=12,height=10)
make_Plots_2(xs=alp,leg=F,b=F);
dev.off()


load(paste(wdd,"dist_sn_All_power_rev_k_exp_alpha_N_100J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots(); 
pdf(file="DvsNDdist_sg_alpha_N_100J_2_rev.pdf",width=12,height=10)
make_Plots_1(xs=alp,leg=F,b=F);
dev.off()
pdf(file="KvOdist_sg_alpha_N_100J_2_rev.pdf",width=12,height=10)
make_Plots_2(xs=alp,leg=F,b=F);
dev.off()


load(paste(wdd,"dist_sn_All_power_rev_k_exp_alpha_N_250J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots(); 


pdf(file="DvsNDdist_sg_alpha_N_250J_2_rev.pdf",width=12,height=10)
make_Plots_1(xs=alp,leg=F,b=F);
dev.off()
pdf(file="KvOdist_sg_alpha_N_250J_2_rev.pdf",width=12,height=10)
make_Plots_2(xs=alp,leg=F,b=F);
dev.off()

#gonna save as 10x12
# KvOdist_sg_alpha_N_250J_2 
# DvsNDdist_sg_alpha_N_250J_2


load(paste(wdd,"dist_sn_All_power_rev_k_exp_beta_N_50J_3_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots();#sizes=cbind(sizes,prop_rejs[1,])
pdf(file="DvsNDdist_sg_beta_N_50J_3_rev.pdf",width=12,height=10)
make_Plots_1();
dev.off()
pdf(file="KvOdist_sg_beta_N_50J_3_rev.pdf",width=12,height=10)
make_Plots_2();
dev.off()

load(paste(wdd,"dist_sn_All_power_rev_k_exp_beta_N_100J_3_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots(); #sizes=cbind(sizes,prop_rejs[1,])
pdf(file="DvsNDdist_sg_beta_N_100J_3_rev.pdf",width=12,height=10)
make_Plots_1();
dev.off()
pdf(file="KvOdist_sg_beta_N_100J_3_rev.pdf",width=12,height=10)
make_Plots_2();
dev.off()



load(paste(wdd,"dist_sn_All_power_rev_k_exp_alpha_N_50J_3_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots(); 
pdf(file="DvsNDdist_sg_alpha_N_50J_3_rev.pdf",width=12,height=10)
make_Plots_1(xs=alp,leg=F,b=F);
dev.off()
pdf(file="KvOdist_sg_alpha_N_50J_3_rev.pdf",width=12,height=10)
make_Plots_2(xs=alp,leg=F,b=F);
dev.off()

load(paste(wdd,"dist_sn_All_power_rev_k_exp_alpha_N_100J_3_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots(); 
pdf(file="DvsNDdist_sg_alpha_N_100J_3_rev.pdf",width=12,height=10)
make_Plots_1(xs=alp,leg=F,b=F);
dev.off()
pdf(file="KvOdist_sg_alpha_N_100J_3_rev.pdf",width=12,height=10)
make_Plots_2(xs=alp,leg=F,b=F);
dev.off()



library(xtable)
sizesT=t(sizes)
xtable(sizes[c(1:4,9:12,17:24),])

xtable(sizes)

