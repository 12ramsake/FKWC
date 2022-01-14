#
#as in Chenouri paper
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
  apply(depths,2,getKW_PM,groups,r=r,k=k,N=N)
}

getKW<-function(depths,N1){
  ranks=rank(depths,ties.method = "average")
  rbar=mean(ranks)
  return((2*N1-1)*(N1*(mean(ranks[1:N1])-rbar)^2+N1*(mean(ranks[(N1+1):length(ranks)])-rbar)^2)/sum((ranks-rbar)^2))
}

#takes matrix of depths, and then gets ts for each
getKW_all_depths<-function(depths,N1){
  apply(depths,2,getKW,N1)
}

setwd("Functional Data Covariance Files")



#J=2 
# Beta Changes
rs=seq(0.25,1,length.out = 15)

print("sim 1")
#changes in beta
prop_rejs=NULL
prop_rejs=array(0,dim=c(4,15,15))
count=0
tmp=seq(0.5,1,length.out = 15)
for(N1 in c(50,100,250)){
  for(dist in c("N","t","sn")){
    for(betai in 1:15){
      #for N1=100 I ran FM seperately
      beta=tmp[betai]
      fileName=paste0("FKWCResults/RPDMBD/dist_",dist,"_deriv_depths_k_exp_alpha1_p05_p05_beta_p",betai,"_p5_N_",N1,"_num_runs_200.Rda",sep="")
      load(fileName)
      depth_values1=depth_values
      
      
      crit=qchisq(.95,1)
      ####PM
      #dim = depth, rs, parameter

      
      for(j in 1:length(rs)){
        test_stats=lapply(depth_values,getKWPM_all_depths,groups=c(1,N1+1),r=rs[j],k=2,N=N1+N1)
        prop_rej=colMeans(do.call(rbind,lapply(test_stats,'>',crit)))
        prop_rejs[1:4,j,which(tmp==beta)]=prop_rej
      }
      # prop_rejs
      
      # names(prop_rej)=c("FMD","RPD","MBD","SPAT","KSPAT","FMD_PM","RPD_PM","MBD_PM","SPAT_PM","KSPAT_PM")
      # 
    }
    print(prop_rejs)
    FN=paste0("Power Files/PM_dist_",dist,"_All_power_k_exp_beta_N_",N1,"J_2_num_runs_50.Rda",sep="")
    save(rs,bet,prop_rejs,file=FN)
  }
}

#
##J=2

# prop_rejs=array(0,dim=c(10,10,7))
library(viridis)
count=0
tmp=seq(0.5,1,length.out = 15)
N1=250
dist="N"
for(N1 in c(50,100,250)){
  for(dist in c("N","t","sn")){
    
    
    FN=paste0("Power Files/PM_dist_",dist,"_All_power_k_exp_beta_N_",N1,"J_2_num_runs_50.Rda",sep="")
    load(file=FN)
    bet=matrix(tmp,ncol=15,nrow=15)
    
    matplot(bet,t(prop_rejs[2,,]),
            type='l',col=viridis(15),
            bty="n",
            lwd=4,
            xlim=c(0.5,0.75),
  #          main=paste0(dist,N1),
            ylab="",yaxt="n",xaxt="n",
            cex.lab=2,
            xlab="")
    mtext("Power",2,line=1,cex=2)
    axis(2,line=-1.5,cex.axis=2)
    axis(1,line=0,cex.axis=2)
    mtext(expression(beta),1,line=1.5,cex=2)
    # prop_rejs=NULL
    legend_image <- as.raster(matrix(rev(viridis(15)), ncol=1))
    rasterImage(legend_image, 0.7, 0.5, 0.71,0.9,cex=2)
    text(x=0.72, y = seq(0.535,.875,l=8), labels = c(round(rs[seq(1,15,2)],2)),cex=1.5,col="black")
    text(x=0.71,labels = expression(italic(r)),y=0.93,cex=2)
  }
}


count=0
tmp=seq(0.05,0.1,length.out=7)

for(N1 in c(50,100,250)){
  for(dist in c("N","t","sn")){
    FN=paste0("Two_Sample_Depths_Cov_Kernels/RPD_MBD/Power Files/PM_",dist,"_All_power_k_exp_alpha_N_",N1,"J_2_num_runs_50.Rda",sep="")
    load(file=FN)
    alp=matrix(tmp,ncol=10,nrow=7)

    matplot(alp,t(prop_rejs[9,,]),type='l',col=colfunc(10),main=paste0(dist,N1))
    # prop_rejs=NULL
  }
}







# ALPHA changes

prop_rejs=array(0,dim=c(10,10,7))
count=0
tmp=seq(0.05,0.1,length.out=7)
for(N1 in c(50,100,250)){
  for(dist in c("N","t","sn")){
    for(alpha in tmp){
      if(N1==100){
        fileName=paste0("RPD_MBD/Regular/dist_",dist,"_deriv_depths_k_exp_alpha1_p05_",alpha,"_beta_p5_p5_N_",N1,"_ot_0_num_runs_50.Rda",sep="")
        load(paste("Two_Sample_Depths_Cov_Kernels/",fileName,sep=""))
        depth_values1=depth_values
        fileName=paste0("RPD_MBD/Just FM/fm_dist_",dist,"_deriv_depths_k_exp_alpha1_p05_",alpha,"_beta_p5_p5_N_",N1,"_ot_0_num_runs_50.Rda",sep="")
        load(paste("Two_Sample_Depths_Cov_Kernels/",fileName,sep=""))
        for(i in 1:length(depth_values))
          depth_values[[i]]=cbind(depth_values[[i]],depth_values1[[i]])
      }
      else{
        fileName=paste0("RPD_MBD/Regular/dist_",dist,"_deriv_depths_k_exp_alpha1_p05_",alpha,"_beta_p5_p5_N_",N1,"_ot_0_num_runs_50.Rda",sep="")
        load(paste("Two_Sample_Depths_Cov_Kernels/",fileName,sep=""))
        
      }
      
      
      
      #spatial
      depth_values1=depth_values
      fileName=paste0("RPD_MBD/Spatial_Added/o_spat_dist_",dist,"_deriv_depths_k_exp_alpha1_p05_",alpha,"_beta_p5_p5_N_",N1,"_ot_0_num_runs_50.Rda",sep="")
      load(paste("Two_Sample_Depths_Cov_Kernels/",fileName,sep=""))
      for(i in 1:length(depth_values))
        depth_values[[i]]=cbind(depth_values1[[i]],depth_values[[i]])
      
      depth_values1=depth_values
      fileName=paste0("RPD_MBD/Spatial_Added/KSPAT/k_spat_dist_",dist,"_deriv_depths_k_exp_alpha1_p05_",alpha,"_beta_p5_p5_N_",N1,"_ot_0_num_runs_50.Rda",sep="")
      load(paste("Two_Sample_Depths_Cov_Kernels/",fileName,sep=""))
      for(i in 1:length(depth_values))
        depth_values[[i]]=cbind(depth_values1[[i]],depth_values[[i]])
      
      
      
      crit=qchisq(.95,1)
      
      for(j in 1:length(rs)){
        test_stats=lapply(depth_values,getKWPM_all_depths,groups=c(1,N1+1),r=rs[j],k=2,N=N1+N1)
        prop_rej=colMeans(do.call(rbind,lapply(test_stats,'>',crit)))
        prop_rejs[1:5,j,which(tmp==alpha)]=prop_rej
      }
      #NO DERIVATIVE
      
      fileName=paste0("RPD_MBD/No derivative/dist_",dist,"_noderiv_depths_k_exp_alpha1_p05_p05_beta_p",beta,"_p5_N_",N1,"_ot_0_num_runs_50.Rda",sep="")
      load(paste("Two_Sample_Depths_Cov_Kernels/",fileName,sep=""))
      
      
      #RP done seperate with simplicial depth
      depth_values1=depth_values
      fileName=paste0("RPD_MBD/RP/RP_dist_",dist,"_noderiv_depths_k_exp_alpha1_p05_",alpha,"_beta_p5_p5_N_",N1,"_ot_0_num_runs_50.Rda",sep="")
      load(paste("Two_Sample_Depths_Cov_Kernels/",fileName,sep=""))
      for(i in 1:length(depth_values)){
        depth_values1[[i]][,2]=depth_values[[i]]
      }
      depth_values=depth_values1
      
      
      #spatial
      depth_values1=depth_values
      fileName=paste0("RPD_MBD/Spatial No Deriv/o_spat_dist_",dist,"_noderiv_depths_k_exp_alpha1_p05_",alpha,"_beta_p5_p5_N_",N1,"_ot_0_num_runs_50.Rda",sep="")
      load(paste("Two_Sample_Depths_Cov_Kernels/",fileName,sep=""))
      for(i in 1:length(depth_values))
        depth_values[[i]]=cbind(depth_values1[[i]],depth_values[[i]])
      
      depth_values1=depth_values
      fileName=paste0("RPD_MBD/Spatial No Deriv/KSPAT/k_spat_dist_",dist,"_noderiv_depths_k_exp_alpha1_p05_",alpha,"_beta_p5_p5_N_",N1,"_ot_0_num_runs_50.Rda",sep="")
      load(paste("Two_Sample_Depths_Cov_Kernels/",fileName,sep=""))
      for(i in 1:length(depth_values))
        depth_values[[i]]=cbind(depth_values1[[i]],depth_values[[i]])
      
      #KW      crit=qchisq(.95,1)
      
      for(j in 1:length(rs)){
        test_stats=lapply(depth_values,getKWPM_all_depths,groups=c(1,N1+1),r=rs[j],k=2,N=N1+N1)
        prop_rej=colMeans(do.call(rbind,lapply(test_stats,'>',crit)))
        prop_rejs[6:10,j,which(tmp==alpha)]=prop_rej
      }
      
      
    }
    alp=matrix(tmp,ncol=ncol(prop_rejs),nrow=nrow(prop_rejs))
    print(prop_rejs)
    FN=paste0("Two_Sample_Depths_Cov_Kernels/RPD_MBD/Power Files/PM_",dist,"_All_power_k_exp_alpha_N_",N1,"J_2_num_runs_50.Rda",sep="")
    save(rs,prop_rejs,file=FN)
    # prop_rejs=NULL
  }
}









#
# 
# ##J=3
# # Beta Changes
# 
# 
# print("sim 1")
# #changes in beta
# prop_rejs=NULL
# count=0
# tmp=seq(0.5,1,length.out = 7)
# for(N1 in c(50,100)){
#   for(dist in c("N","t","sn")){
#     for(beta in tmp){
#       
#       fileName=paste0("RPD_MBD/Regular/dist_",dist,"_deriv_depths_k_exp_alpha1_p05_p05_p05_beta_p5_p",beta,"_p5_N_",N1,"_ot_0_num_runs_50.Rda",sep="")
#       load(paste("Two_Sample_Depths_Cov_Kernels/",fileName,sep=""))
#       
#       
#       
#       
#       
#       #Spatial Depth
#       depth_values1=depth_values
#       fileName=paste0("RPD_MBD/Spatial_Added/o_spat_dist_",dist,"_deriv_depths_k_exp_alpha1_p05_p05_p05_beta_p5_p",beta,"_p5_N_",N1,"_ot_0_num_runs_50.Rda",sep="")
#       load(paste("Two_Sample_Depths_Cov_Kernels/",fileName,sep=""))
#       for(i in 1:length(depth_values))
#         depth_values[[i]]=cbind(depth_values1[[i]],depth_values[[i]])
#       
#       
#       depth_values1=depth_values
#       fileName=paste0("RPD_MBD/Spatial_Added/KSPAT/k_spat_dist_",dist,"_deriv_depths_k_exp_alpha1_p05_p05_p05_beta_p5_p",beta,"_p5_N_",N1,"_ot_0_num_runs_50.Rda",sep="")
#       load(paste("Two_Sample_Depths_Cov_Kernels/",fileName,sep=""))
#       for(i in 1:length(depth_values))
#         depth_values[[i]]=cbind(depth_values1[[i]],depth_values[[i]])
#       
#       
#       
#       
#       #KW
#       test_stats=lapply(depth_values,getKW_all_depths,N1=N1)
#       
#       crit=qchisq(.95,1)
#       
#       prop_rej=colMeans(do.call(rbind,lapply(test_stats,'>',crit)));       prop_rej
#       
#       
#       
#       ####PM
#       test_stats=lapply(depth_values,getKWPM_all_depths,groups=c(1,N1+1),r=0.8,k=2,N=N1+N1)
#       
#       crit=qchisq(.95,1)
#       
#       prop_rej=c(prop_rej,colMeans(do.call(rbind,lapply(test_stats,'>',crit))))
#       
#       names(prop_rej)=c("FMD","RPD","MBD","SPAT","KSPAT","FMD_PM","RPD_PM","MBD_PM","SPAT_PM","KSPAT_PM")
#       
#       
#       
#       
#       fileName=paste0("RPD_MBD/No derivative/dist_",dist,"_noderiv_depths_k_exp_alpha1_p05_p05_beta_p5_p",beta,"_p5_N_",N1,"_ot_0_num_runs_50.Rda",sep="")
#       load(paste("Two_Sample_Depths_Cov_Kernels/",fileName,sep=""))
#       
#       #RP done seperate with simplicial depth
#       depth_values1=depth_values
#       fileName=paste0("RPD_MBD/RP/RP_dist_",dist,"_noderiv_depths_k_exp_alpha1_p05_p05_beta_p5_p",beta,"_p5_N_",N1,"_ot_0_num_runs_50.Rda",sep="")
#       load(paste("Two_Sample_Depths_Cov_Kernels/",fileName,sep=""))
#       for(i in 1:length(depth_values)){
#         depth_values1[[i]][,2]=depth_values[[i]]
#       }
#       depth_values=depth_values1
#       
#       
#       #Spatial Depth
#       #Spatial Depth
#       depth_values1=depth_values
#       fileName=paste0("RPD_MBD/Spatial No Deriv/o_spat_dist_",dist,"_noderiv_depths_k_exp_alpha1_p05_p05_beta_p5_p",beta,"_p5_N_",N1,"_ot_0_num_runs_50.Rda",sep="")
#       load(paste("Two_Sample_Depths_Cov_Kernels/",fileName,sep=""))
#       for(i in 1:length(depth_values))
#         depth_values[[i]]=cbind(depth_values1[[i]],depth_values[[i]])
#       
#       
#       depth_values1=depth_values
#       fileName=paste0("RPD_MBD/Spatial No Deriv/KSPAT/k_spat_dist_",dist,"_noderiv_depths_k_exp_alpha1_p05_p05_beta_p5_p",beta,"_p5_N_",N1,"_ot_0_num_runs_50.Rda",sep="")
#       load(paste("Two_Sample_Depths_Cov_Kernels/",fileName,sep=""))
#       for(i in 1:length(depth_values))
#         depth_values[[i]]=cbind(depth_values1[[i]],depth_values[[i]])
#       
#       
#       #KW
#       test_stats=lapply(depth_values,getKW_all_depths,N1=N1)
#       
#       crit=qchisq(.95,1)
#       
#       prop_rej2=colMeans(do.call(rbind,lapply(test_stats,'>',crit)))
#       prop_rej=c(prop_rej, prop_rej2)
#       prop_rej
#       
#       
#       ####PM
#       test_stats=lapply(depth_values,getKWPM_all_depths,groups=c(1,N1+1),r=0.8,k=2,N=N1+N1)
#       
#       crit=qchisq(.95,1)
#       
#       prop_rej=c(prop_rej,colMeans(do.call(rbind,lapply(test_stats,'>',crit))))
#       
#       names(prop_rej)=c(names(prop_rej)[1:10],
#                         paste(c("FMD","RPD","MBD","SPAT","KSPAT",
#                                 "FMD_PM","RPD_PM","MBD_PM","SPAT_PM","KSPAT_PM"),"_ND",sep=""))
#       
#       
#       
#       
#       ##other tests, had to run guo seperately here
#       fileName=paste0("Two_Sample_Depths_Cov_Kernels/RPD_MBD/Other Tests/dist_",dist,"_other_k_exp_alpha1_p05_p05_beta_p5_p",beta,"_p5_N_",N1,"_ot_0_num_runs_50.Rda",sep="")
#       
#       load(fileName)
#       pvalues=do.call(rbind, p_values)
#       
#       
#       prop_rej=c(prop_rej,apply(pvalues<0.05,2,mean))
#       prop_rejs=rbind(prop_rejs,prop_rej)
#       
#       
#     }
#     bet=matrix(tmp,ncol=ncol(prop_rejs),nrow=nrow(prop_rejs))
#     print(prop_rejs)
#     FN=paste0("Two_Sample_Depths_Cov_Kernels/RPD_MBD/Power Files/dist_",dist,"_All_power_k_exp_beta_N_",N1,"J_3_num_runs_50.Rda",sep="")
#     save(bet,prop_rejs,file=FN)
#     prop_rejs=NULL
#   }
# }
# 
# 
# 
# 
# 
# print("sim 2")
# #changes in alpha
# prop_rejs=NULL
# count=0
# tmp=seq(0.05,0.1,length.out=7)
# for(N1 in c(50,100)){
#   for(dist in c("N","t","sn")){
#     for(alpha in tmp){
#       
#       fileName=paste0("RPD_MBD/Regular/dist_",dist,"_deriv_depths_k_exp_alpha1_p05_",alpha,"_p05_beta_p5_p5_p5_N_",N1,"_ot_0_num_runs_50.Rda",sep="")
#       load(paste("Two_Sample_Depths_Cov_Kernels/",fileName,sep=""))
#       
#       
#       
#       
#       
#       #Spatial Depth
#       depth_values1=depth_values
#       fileName=paste0("RPD_MBD/Spatial_Added/o_spat_dist_",dist,"_deriv_depths_k_exp_alpha1_p05_",alpha,"_p05_beta_p5_p5_N_",N1,"_ot_0_num_runs_50.Rda",sep="")
#       load(paste("Two_Sample_Depths_Cov_Kernels/",fileName,sep=""))
#       for(i in 1:length(depth_values))
#         depth_values[[i]]=cbind(depth_values1[[i]],depth_values[[i]])
#       
#       
#       depth_values1=depth_values
#       fileName=paste0("RPD_MBD/Spatial_Added/KSPAT/k_spat_dist_",dist,"_deriv_depths_k_exp_alpha1_p05_",alpha,"_p05_beta_p5_p5_N_",N1,"_ot_0_num_runs_50.Rda",sep="")
#       load(paste("Two_Sample_Depths_Cov_Kernels/",fileName,sep=""))
#       for(i in 1:length(depth_values))
#         depth_values[[i]]=cbind(depth_values1[[i]],depth_values[[i]])
#       
#       
#       
#       
#       #KW
#       test_stats=lapply(depth_values,getKW_all_depths,N1=N1)
#       
#       crit=qchisq(.95,1)
#       
#       prop_rej=colMeans(do.call(rbind,lapply(test_stats,'>',crit)));       prop_rej
#       
#       
#       
#       ####PM
#       test_stats=lapply(depth_values,getKWPM_all_depths,groups=c(1,N1+1),r=0.8,k=2,N=N1+N1)
#       
#       crit=qchisq(.95,1)
#       
#       prop_rej=c(prop_rej,colMeans(do.call(rbind,lapply(test_stats,'>',crit))))
#       
#       names(prop_rej)=c("FMD","RPD","MBD","SPAT","KSPAT","FMD_PM","RPD_PM","MBD_PM","SPAT_PM","KSPAT_PM")
#       
#       
#       
#       
#       fileName=paste0("RPD_MBD/No derivative/dist_",dist,"_noderiv_depths_k_exp_alpha1_p05_",alpha,"_p05_beta_p5_p5_N_",N1,"_ot_0_num_runs_50.Rda",sep="")
#       load(paste("Two_Sample_Depths_Cov_Kernels/",fileName,sep=""))
#       
#       #RP done seperate with simplicial depth
#       depth_values1=depth_values
#       fileName=paste0("RPD_MBD/RP/RP_dist_",dist,"_noderiv_depths_k_exp_alpha1_p05_",alpha,"_p05_beta_p5_p5_N_",N1,"_ot_0_num_runs_50.Rda",sep="")
#       load(paste("Two_Sample_Depths_Cov_Kernels/",fileName,sep=""))
#       for(i in 1:length(depth_values)){
#         depth_values1[[i]][,2]=depth_values[[i]]
#       }
#       depth_values=depth_values1
#       
#       
#       #Spatial Depth
#       #Spatial Depth
#       depth_values1=depth_values
#       fileName=paste0("RPD_MBD/Spatial No Deriv/o_spat_dist_",dist,"_noderiv_depths_k_exp_alpha1_p05_",alpha,"_p05_beta_p5_p5_N_",N1,"_ot_0_num_runs_50.Rda",sep="")
#       load(paste("Two_Sample_Depths_Cov_Kernels/",fileName,sep=""))
#       for(i in 1:length(depth_values))
#         depth_values[[i]]=cbind(depth_values1[[i]],depth_values[[i]])
#       
#       
#       depth_values1=depth_values
#       fileName=paste0("RPD_MBD/Spatial No Deriv/KSPAT/k_spat_dist_",dist,"_noderiv_depths_k_exp_alpha1_p05_",alpha,"_p05_beta_p5_p5_N_",N1,"_ot_0_num_runs_50.Rda",sep="")
#       load(paste("Two_Sample_Depths_Cov_Kernels/",fileName,sep=""))
#       for(i in 1:length(depth_values))
#         depth_values[[i]]=cbind(depth_values1[[i]],depth_values[[i]])
#       
#       
#       #KW
#       test_stats=lapply(depth_values,getKW_all_depths,N1=N1)
#       
#       crit=qchisq(.95,1)
#       
#       prop_rej2=colMeans(do.call(rbind,lapply(test_stats,'>',crit)))
#       prop_rej=c(prop_rej, prop_rej2)
#       prop_rej
#       
#       
#       ####PM
#       test_stats=lapply(depth_values,getKWPM_all_depths,groups=c(1,N1+1),r=0.8,k=2,N=N1+N1)
#       
#       crit=qchisq(.95,1)
#       
#       prop_rej=c(prop_rej,colMeans(do.call(rbind,lapply(test_stats,'>',crit))))
#       
#       names(prop_rej)=c(names(prop_rej)[1:10],
#                         paste(c("FMD","RPD","MBD","SPAT","KSPAT",
#                                 "FMD_PM","RPD_PM","MBD_PM","SPAT_PM","KSPAT_PM"),"_ND",sep=""))
#       
#       
#       
#       
#       ##other tests, had to run guo seperately here
#       fileName=paste0("Two_Sample_Depths_Cov_Kernels/RPD_MBD/Other Tests/dist_",dist,"_other_k_exp_alpha1_p05_",alpha,"_p05_beta_p5_p5_N_",N1,"_ot_0_num_runs_50.Rda",sep="")
#       
#       load(fileName)
#       pvalues=do.call(rbind, p_values)
#       
#       
#       prop_rej=c(prop_rej,apply(pvalues<0.05,2,mean))
#       prop_rejs=rbind(prop_rejs,prop_rej)
#       
#       
#     }
#     bet=matrix(tmp,ncol=ncol(prop_rejs),nrow=nrow(prop_rejs))
#     print(prop_rejs)
#     FN=paste0("Two_Sample_Depths_Cov_Kernels/RPD_MBD/Power Files/dist_",dist,"_All_power_k_exp_alpha_N_",N1,"J_3_num_runs_50.Rda",sep="")
#     save(bet,prop_rejs,file=FN)
#     prop_rejs=NULL
#   }
# }
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# library(xtable)
# #Eigenvalue Decay
# prop_rejs=NULL
# scens=c("1","2","3","21","22","23")
# for(N1 in c(50,100)){
#   dist="N"
#   for(scen in scens){
#     
#     
#     fileName=paste0("Two_Sample_Depths_Cov_Kernels/RPD_MBD/Eigen Scenario/allD_dist_N_deriv_depths_k_eigen_",scen,"_samesum_N_",N1,"_ot_0_num_runs_50.Rda",sep="")
#     load(fileName)
#     
#     
#     #KW
#     test_stats=lapply(depth_values,getKW_all_depths,N1=N1)
#     
#     crit=qchisq(.95,1)
#     
#     prop_rej=colMeans(do.call(rbind,lapply(test_stats,'>',crit)));       prop_rej
#     
#     
#     
#     ####PM
#     test_stats=lapply(depth_values,getKWPM_all_depths,groups=c(1,N1+1),r=0.8,k=2,N=N1+N1)
#     
#     crit=qchisq(.95,1)
#     
#     prop_rej=c(prop_rej,colMeans(do.call(rbind,lapply(test_stats,'>',crit))))
#     
#     names(prop_rej)=c("FMD","RPD","MBD","SPAT","KSPAT","FMD_PM","RPD_PM","MBD_PM","SPAT_PM","KSPAT_PM")
#     
#     
#     fileName=paste0("Two_Sample_Depths_Cov_Kernels/RPD_MBD/Eigen Scenario/RP_dist_N_noderiv_depths_k_eigen_",scen,"_samesum_N_",N1,"_ot_0_num_runs_50.Rda",sep="")
#     load(fileName)
#     rp=depth_values
#     fileName=paste0("Two_Sample_Depths_Cov_Kernels/RPD_MBD/Eigen Scenario/allD_dist_N_noderiv_depths_k_eigen_",scen,"_samesum_N_",N1,"_ot_0_num_runs_50.Rda",sep="")
#     load(fileName)
#     for(i in 1:length(depth_values)){
#       depth_values[[i]][,2]=rp[[i]]
#     }
#     
#     
#     
#     #KW
#     test_stats=lapply(depth_values,getKW_all_depths,N1=N1)
#     
#     crit=qchisq(.95,1)
#     
#     prop_rej=c(prop_rej,colMeans(do.call(rbind,lapply(test_stats,'>',crit))));       prop_rej
#     
#     
#     
#     ####PM
#     test_stats=lapply(depth_values,getKWPM_all_depths,groups=c(1,N1+1),r=0.8,k=2,N=N1+N1)
#     
#     crit=qchisq(.95,1)
#     
#     prop_rej=c(prop_rej,colMeans(do.call(rbind,lapply(test_stats,'>',crit))))
#     
#     names(prop_rej)=c(names(prop_rej)[1:10],
#                       paste(c("FMD","RPD","MBD","SPAT","KSPAT","FMD_PM","RPD_PM","MBD_PM","SPAT_PM","KSPAT_PM"),"_ND",sep=""))
#     
#     
#     
#     
#     fileName=paste0("Two_Sample_Depths_Cov_Kernels/RPD_MBD/Eigen Scenario/other_dist_N_k_eigen_",scen,"_samesum_N_",N1,"_ot_0_num_runs_50.Rda",sep="")
#     load(fileName)
#     pvalues=do.call(rbind, p_values)
#     prop_rej=c(prop_rej,apply(pvalues<0.05,2,mean))
#     prop_rejs=rbind(prop_rejs,prop_rej)
#     
#     
#   }
#   print(prop_rejs)
#   FN=paste0("Two_Sample_Depths_Cov_Kernels/RPD_MBD/Power Files/eigen_",dist,"_All_power_N_",N1,"_num_runs_50.Rda",sep="")
#   xtable(t(prop_rejs))
#   save(scens,prop_rejs,file=FN)
#   prop_rejs=NULL
# }
# 
# 
# 
# 
# 
