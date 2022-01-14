#
#as in CHenouri paper
getTj<-function(j,groups,ranks,nprime){
  ranks=ranks[groups[j]:(groups[j+1]-1)]
  sum(nprime-ranks[ranks<=nprime]+1)
}

getKj<-function(Tj,nj,nprime,N){
  ETj=nj*nprime*(nprime+1)/(2*N)
  VTj=(nj*(N-nj)*nprime*(nprime+1)*(2*N*(2*nprime+1)-
                                      3*nprime*(nprime+1)))/(12*N^2*(N-1))
  Kj=(Tj-ETj)^2/VTj
  return(Kj)
}

#groups is the group indices at which each sample starts and N+1 is at the end
getKW_PM<-function(depths,groups,r=0.8,k=2,N){
  nprime=floor(r*N)
  ranks=rank(depths,ties.method = "average")
  Tjs<-sapply(1:k,getTj,groups=c(groups,N+1),ranks=ranks,nprime=nprime)
  njs=c(groups[2:k],N+1)-groups
  Kjs=mapply(getKj,Tjs,njs,MoreArgs = list(nprime=nprime,N=N))
  return(sum((1-njs/N)*Kjs))
}

#takes matrix of depths, and then gets ts for each
getKWPM_all_depths<-function(depths,groups,r=0.8,k=2,N){
  apply(depths,2,getKW_PM,groups,r=r,k=2,N=N)
}
# par(mfrow=c(1,5))

prop_rejs=NULL
N1=N2=250
###MVsims
# pdf(file=paste0("Simulation Study Power Graph Results 2020/mv_power_graphs.pdf",sep=""),height=10,width=23)
cols=c("red","black","blue","green","purple")
print("sim 1")
#changes in beta
count=0
tmp=seq(0.5,1,length.out = 15)
for(dist in c("N","t","sn")){
  for(betai in 1:15){
    fileName=paste0("RPD_MBD/dist_",dist,"_deriv_depths_k_exp_alpha1_p05_p05_beta_p",beta,"_p5_N_",N1,"_ot_0_num_runs_50.Rda",sep="")
    load(paste("",fileName,sep=""))
    depth_values1=depth_values
    fileName=paste0("RPD_MBD/fm_dist_",dist,"_deriv_depths_k_exp_alpha1_p05_p05_beta_p",beta,"_p5_N_",N1,"_ot_0_num_runs_50.Rda",sep="")
    load(paste("",fileName,sep=""))
    for(i in 1:length(depth_values))
      depth_values[[i]]=cbind(depth_values1[[i]],depth_values[[i]])
    depth_values1=depth_values
    fileName=paste0("RPD_MBD/spat_dist_",dist,"_deriv_depths_k_exp_alpha1_p05_p05_beta_p",beta,"_p5_N_",N1,"_ot_0_num_runs_50.Rda",sep="")
    load(paste("",fileName,sep=""))
    
    
    for(i in 1:length(depth_values))
      depth_values[[i]]=cbind(depth_values1[[i]],depth_values[[i]])

    
    test_stats=lapply(depth_values,getKWPM_all_depths,groups=c(1,N1+1),r=0.8,k=2,N=N1+N2)

    crit=qchisq(.95,1)
    
    prop_rej=colMeans(do.call(rbind,lapply(test_stats,'>',crit)))
    prop_rej
    names(prop_rej)=c("RPD","MBD","FMD","SPAT","KSPAT")
    prop_rejs=rbind(prop_rejs,prop_rej)
  }
  matplot(matrix(tmp,ncol=ncol(prop_rejs),nrow=nrow(prop_rejs)),prop_rejs,type='l',lwd=2,lty=rep(c(2,3),each=8),xlim=c(min(tmp),max(tmp)+.1),col=cols,xlab="Beta",ylab="prop null rej",main=paste0("dist " ,dist))
  legend(x="topright",y=NULL,lty=rep(c(2,3),each=8),col=cols,legend = names(prop_rej))
  # legend(x="topright",y=NULL,lty=rep(c(2,3),each=8),col=1,legend = names(prop_rej))
  # dev.off()
  print(  prop_rejs)
  prop_rejs=NULL
}
dev.off()

print("sim 2")
#changes in alpha

tmp=seq(0.05,0.1,length.out=7)
for(dist in c("N","t","sn")){
  for(alpha in tmp){
    fileName=paste0("RPD_MBD/dist_",dist,"_deriv_depths_k_exp_alpha1_p05_",alpha,"_beta_p5_p5_N_",N1,"_ot_0_num_runs_50.Rda",sep="")
    load(paste("Two_Sample_Depths_Cov_Kernels/",fileName,sep=""))
    depth_values1=depth_values
    fileName=paste0("RPD_MBD/fm_dist_",dist,"_deriv_depths_k_exp_alpha1_p05_",alpha,"_beta_p5_p5_N_",N1,"_ot_0_num_runs_50.Rda",sep="")
    load(paste("Two_Sample_Depths_Cov_Kernels/",fileName,sep=""))
    
    for(i in 1:length(depth_values))
      depth_values[[i]]=cbind(depth_values[[i]],depth_values1[[i]])

    
    test_stats=lapply(depth_values,getKWPM_all_depths,groups=c(1,N1+1),r=0.85,k=2,N=N1+N2)
    
    crit=qchisq(.95,1)
    
    prop_rej=colMeans(do.call(rbind,lapply(test_stats,'>',crit)))
    prop_rej
    names(prop_rej)=c("FMD","RPD","MBD")
    
    
    test_stats=lapply(depth_values,getKWPM_all_depths,groups=c(1,N1+1),r=0.8,k=2,N=N1+N2)
    
    crit=qchisq(.95,1)
    
    prop_rej=c(prop_rej,colMeans(do.call(rbind,lapply(test_stats,'>',crit))))
    prop_rej
    
    test_stats=lapply(depth_values,getKWPM_all_depths,groups=c(1,N1+1),r=0.3,k=2,N=N1+N2)
    
    crit=qchisq(.95,1)
    
    prop_rej=c(prop_rej,colMeans(do.call(rbind,lapply(test_stats,'>',crit))))
    prop_rej
    

    names(prop_rej)=c("FMD8","RPD8","MBD8","FMD5","RPD5","MBD5","FMD3","RPD3","MBD3")
    
    prop_rejs=rbind(prop_rejs,prop_rej)
  }
  matplot(matrix(tmp,ncol=ncol(prop_rejs),nrow=nrow(prop_rejs)),prop_rejs,type='l',lwd=2,
          lty=rep(c(2,3),each=8),xlim=c(min(tmp),max(tmp)+.01),col=c(rep(1,3),rep(2,3),rep(3,3)),xlab="Alpha",ylab="prop null rej",main=paste0("dist " ,dist))
  legend(x="topright",y=NULL,lty=rep(c(2,3),each=8),col=cols,legend = names(prop_rej)) # legend(x="topright",y=NULL,lty=rep(c(2,3),each=8),col=1,legend = names(prop_rej))
  # dev.off()
  prop_rejs=NULL
}
    
    
dev.off()