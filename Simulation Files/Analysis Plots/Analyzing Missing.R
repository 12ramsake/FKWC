
wdd="C:/Users/12RAM/OneDrive/Documents/research/PhD Thesis/Functional Data Covariance Files/Nov_2020_sim_results/Power Files/"
graph_folder="C:/Users/12RAM/OneDrive/Documents/research/PhD Thesis/Functional Data Covariance Files/Nov_2020_sim_results/POWER_GRAPHS_REVISION/"

pink_grad=colorRampPalette(c("#470f32","#821e5c", "#d63696","#ffd1ee"))
cols=c(RColorBrewer::brewer.pal(9, "Blues")[5],"#ffd1ee")
library(RColorBrewer)
cols=c("black",brewer.pal(9, "Blues")[7])

make_Plots_2=function(xs=bet,leg=T,b=T){
  #depth plots
  colnames(prop_rejs)[1]="MFHD"
  colnames(prop_rejs)[5]="RPD-L"
  # colnames(prop_rejs)[5]="KSD"
  indz=c(1:5,6,21:23-13,25:26-13)
  
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


make_plot=function(fn,xs=bet,leg=F,b=T){
  Cairo::CairoPDF(fn,height=10,width=17)
  make_Plots_2(xs,leg,b); 
  dev.off()
}
######################################## Regular Scenario  ################################################

######################################### 2 Groups #####################################
N1=100
dist='N'
load(paste(wdd,"missing_dist_",dist,"_All_power_rev_k_exp_beta_N_",N1,"J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots_2()
make_plot(paste(wdd,"missing_dist_",dist,"_All_power_rev_k_exp_beta_N_",N1,"J_2_num_runs_200.pdf",sep=""),b=T)
sizes=prop_rejs[1,]

dist='t'
load(paste(wdd,"missing_dist_",dist,"_All_power_rev_k_exp_beta_N_",N1,"J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots_2()
make_plot(paste(wdd,"missing_dist_",dist,"_All_power_rev_k_exp_beta_N_",N1,"J_2_num_runs_200.pdf",sep=""),b=T)
sizes=rbind(prop_rejs[1,],sizes)

dist='N'
load(paste(wdd,"missing_dist_",dist,"_All_power_rev_k_exp_alpha_N_",N1,"J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots_2()
make_plot(paste(wdd,"missing_dist_",dist,"_All_power_rev_k_exp_alpha_N_",N1,"J_2_num_runs_200.pdf",sep=""),b=F,xs=alp)
sizes=rbind(prop_rejs[1,],sizes)

dist='t'
load(paste(wdd,"missing_dist_",dist,"_All_power_rev_k_exp_alpha_N_",N1,"J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
make_Plots_2()
make_plot(paste(wdd,"missing_dist_",dist,"_All_power_rev_k_exp_alpha_N_",N1,"J_2_num_runs_200.pdf",sep=""),b=F,xs=alp)
sizes=rbind(prop_rejs[1,],sizes)







