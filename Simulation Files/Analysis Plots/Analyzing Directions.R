
wdd="C:/Users/12RAM/OneDrive/Documents/research/PhD Thesis/Functional Data Covariance Files/Nov_2020_sim_results/Power Files/"
graph_folder="C:/Users/12RAM/OneDrive/Documents/research/PhD Thesis/Functional Data Covariance Files/Nov_2020_sim_results/POWER_GRAPHS_REVISION/"

pink_grad=colorRampPalette(c("#470f32","#821e5c", "#d63696","#ffd1ee"))
cols=c(RColorBrewer::brewer.pal(9, "Blues")[5],"#ffd1ee")
library(RColorBrewer)
cols=c("black",brewer.pal(9, "Blues")[7])

make_Plots_2=function(xs=bet,leg=T,b=T){
  #depth plots


  matplot(xs,prop_rejs,type='l',
          col=rep(rainbow(4),each=2),
          lty=rep(c(1,2),4),
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
  if(leg){
    legend(xs[9,1],0.28,legend=paste0(c(20,50,75,100)," directions"),lty=1,col=rep(rainbow(4),each=1),cex=3,lwd=1,bty = "n")
    legend(xs[9,1],0.45,legend=c("Simplicial","Likelihood"),lty=1:2,cex=3,lwd=1,bty = "n")
  }
    # legend("bottomright",legend=colnames(prop_rejs),lty=rep(c(1,2),4),
    #        col=rep(rainbow(4),each=2),cex=3,lwd=1,bty = "n")
  # legend("bottomright",legend=c("W","Other"),col=cols,lty=1)
  
}
make_Plots_2(xs=alp)



load_stuff_beta=function(dist,N1){
  load(paste(wdd,"dist_",dist,"_All_power_rev_k_exp_beta_N_",N1,"J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
  pr1=prop_rejs[,c(2,5)]
  load(paste(wdd,"Directions_dist_",dist,"_All_power_rev_k_exp_beta_N_",N1,"J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
  prop_rejs=cbind(pr1,prop_rejs)
  colnames(prop_rejs)[1:2]=c('RPD_20','RPD_Like_20')
  return(prop_rejs)
}
load_stuff_alpha=function(dist,N1){
  load(paste(wdd,"dist_",dist,"_All_power_rev_k_exp_alpha_N_",N1,"J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
  pr1=prop_rejs[,c(2,5)]
  load(paste(wdd,"Directions_dist_",dist,"_All_power_rev_k_exp_alpha_N_",N1,"J_2_num_runs_200.Rda",sep="")); dim(prop_rejs)
  prop_rejs=cbind(pr1,prop_rejs)
  colnames(prop_rejs)[1:2]=c('RPD_20','RPD_Like_20')
  return(prop_rejs)
}


make_plot=function(fn,xs=bet,leg=F,b=T){
  Cairo::CairoPDF(fn,height=10,width=14)
  make_Plots_2(xs,leg,b); 
  dev.off()
}
######################################## Regular Scenario  ################################################

######################################### 2 Groups #####################################
N1=50
dist='N'
prop_rejs=load_stuff_beta('N',50)
make_Plots_2(b=T)
make_plot(paste(wdd,"Directions_dist_",dist,"_All_power_rev_k_exp_beta_N_",N1,"_.pdf",sep=""),leg=T,b=T)
sizes=prop_rejs[1,]

N1=50
dist='t'
prop_rejs=load_stuff_beta('t',50)
make_Plots_2(b=T)
make_plot(paste(wdd,"Directions_dist_",dist,"_All_power_rev_k_exp_beta_N_",N1,"_.pdf",sep=""),b=T)
sizes=prop_rejs[1,]


N1=50
dist='sn'
prop_rejs=load_stuff_beta('sn',50)
make_Plots_2(b=T)
make_plot(paste(wdd,"Directions_dist_",dist,"_All_power_rev_k_exp_beta_N_",N1,"_.pdf",sep=""),b=T)
sizes=prop_rejs[1,]




N1=50
dist='N'
prop_rejs=load_stuff_alpha('N',50)
make_Plots_2(xs=alp)
make_plot(paste(wdd,"Directions_dist_",dist,"_All_power_rev_k_exp_alpha_N_",N1,"_.pdf",sep=""),b=F,xs=alp)
sizes=prop_rejs[1,]

N1=50
dist='t'
prop_rejs=load_stuff_alpha('t',50)
make_Plots_2(xs=alp)
make_plot(paste(wdd,"Directions_dist_",dist,"_All_power_rev_k_exp_alpha_N_",N1,"_.pdf",sep=""),b=F,xs=alp)
sizes=prop_rejs[1,]


N1=50
dist='sn'
prop_rejs=load_stuff_alpha('sn',50)
make_Plots_2(xs=alp)
make_plot(paste(wdd,"Directions_dist_",dist,"_All_power_rev_k_exp_alpha_N_",N1,"_.pdf",sep=""),b=F,xs=alp)
sizes=prop_rejs[1,]





