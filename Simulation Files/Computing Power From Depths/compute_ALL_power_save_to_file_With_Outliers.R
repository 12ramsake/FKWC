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
# 
# getKW<-function(depths,N1){
#   ranks=rank(depths,ties.method = "average")
#   rbar=mean(ranks)
#   return((2*N1-1)*(N1*(mean(ranks[1:N1])-rbar)^2+N1*(mean(ranks[(N1+1):length(ranks)])-rbar)^2)/sum((ranks-rbar)^2))
# }
# 
# #takes matrix of depths, and then gets ts for each
# getKW_all_depths<-function(depths,N1){
#   apply(depths,2,getKW,N1)
# }

setwd("Functional Data Covariance Files/Nov_2020_sim_results")



# a=NULL
# for(i in 1:200){
#   a=c(a,getKW(depth_values[[i]][,1],N1)>qchisq(.95,1))
# }
# 
# mean(a)
########################################################### Beta ###############################################################
###########################################################  J=2 ###############################################################

dist='N'
print("sim 1")
#changes in beta
prop_rejs=NULL
count=0
tmp=seq(0.5,1,length.out = 15)
for(N1 in c(50,100,250)){
  N2=N1
    for(out in 1:4){
      for(both in c(T,F)){
        for(betai in 1:15){
      
      
      #for N1=100 I ran LTR separately
      beta=tmp[betai]
      fileName=paste0("REVISION_RESULTS/LTR/cl_",cl,"_both_",both,"_out_",out,"_dist_",dist,"_LTR_depths_k_exp_alpha1_p05_p05_beta_p",betai,"_p5_N_",N1,"_num_runs_200.Rda",sep="")
      #REPLACE WITH REVISED DEPTHZ:
      load(fileName)
      depth_values1=depth_values
      
      fileName=paste0("REVISION_RESULTS/outlier/cl_",cl,"_both_",both,"_out_",out,"_dist_", dist,"_deriv_SIM_TD_k_exp_alpha1_p05_p05_beta_p", betai, "_p5_N_",N1,"_num_runs_200.Rda", sep = "")
      load(fileName)

      
      for(i in 1:num_runs){
        depth_values[[i]][,4]= depth_values1[[i]]
      }
      
      
      
      
      #KW
      test_stats=lapply(depth_values,getKWPM_all_depths,groups=c(1,N1+1),r=1,k=2,N=N1+N2)
      
      crit=qchisq(.95,1)
      
      prop_rej=colMeans(do.call(rbind,lapply(test_stats,'>',crit)));       prop_rej
      
      
      
      ####PM
      test_stats=lapply(depth_values,getKWPM_all_depths,groups=c(1,N1+1),r=0.8,k=2,N=N1+N2)
      
      crit=qchisq(.95,1)
      
      prop_rej=c(prop_rej,colMeans(do.call(rbind,lapply(test_stats,'>',crit))))
      
      # names(prop_rej)=c("FMD","RPD","MBD","LTR","SPAT","KSPAT","FMD_PM","RPD_PM","MBD_PM","LTR","SPAT_PM","KSPAT_PM")
      names(prop_rej)=c("FMD","RPD","MBD","LTR","RPD_Like","FMD_PM","RPD_PM","MBD_PM","LTR_PM","RPD_Like_PM")
      
      
      
      
      fileName=paste0("REVISION_RESULTS/outlier/cl_",cl,"_both_",both,"_out_",out,"_dist_", dist,"_noderiv_depths_k_exp_alpha1_p05_p05_beta_p",betai,"_p5_N_",N1,"_num_runs_200.Rda",sep="")
      load(fileName)
      depth_values1=depth_values
      
      
      
      
      #KW
      test_stats=lapply(depth_values,getKWPM_all_depths,groups=c(1,N1+1),r=1,k=2,N=N1+N1)
      
      crit=qchisq(.95,1)
      
      prop_rej2=colMeans(do.call(rbind,lapply(test_stats,'>',crit)))
      prop_rej=c(prop_rej, prop_rej2)
      prop_rej
      
      
      ####PM
      test_stats=lapply(depth_values,getKWPM_all_depths,groups=c(1,N1+1),r=0.8,k=2,N=N1+N1)
      
      crit=qchisq(.95,1)
      
      prop_rej=c(prop_rej,colMeans(do.call(rbind,lapply(test_stats,'>',crit))))
      
      # names(prop_rej)=c(names(prop_rej)[1:10],
      #                     paste(c("FMD","RPD","MBD","LTR","SPAT","KSPAT",
      #                             "FMD_PM","RPD_PM","MBD_PM","LTR_PM","SPAT_PM","KSPAT_PM"),"_ND",sep=""))
      names(prop_rej)=c(names(prop_rej)[1:10],
                        paste(c("FMD","RPD","MBD","LTR",
                                "FMD_PM","RPD_PM","MBD_PM","LTR_PM"),"_ND",sep=""))
      
      
      
      
      ##other tests, had to run guo seperately here
      
      fileName=paste0("REVISION_RESULTS/outlier/cl_",cl,"_both_",both,"_out_",out,"dist_",dist,"_other_k_exp_alpha1_p05_p05_beta_p",betai,"_p5_N_",N1,"_num_runs_200.Rda",sep="")      
      a=try({load(fileName)},T)
      errorsp=inherits(a, "try-error")
      if(errorsp){
        fileName=paste0("REVISION_RESULTS/outlier/cl_",cl,"both_",both,"_out_",out,"_dist_",dist,"_other_k_exp_alpha1_p05_p05_beta_p",betai,"_p5_N_",N1,"_num_runs_200.Rda",sep="")      
        load(fileName)
      }
      pvalues=do.call(rbind, p_values)
      
      
      prop_rej=c(prop_rej,apply(pvalues<0.05,2,mean))
      prop_rejs=rbind(prop_rejs,prop_rej)
      
      
    }
    bet=matrix(tmp,ncol=ncol(prop_rejs),nrow=nrow(prop_rejs))
    print(prop_rejs)
    FN=paste0("Power Files/cl_",cl,"_both_",both,"_out_",out,"dist_",dist,"_All_power_rev_k_exp_beta_N_",N1,"J_2_num_runs_200.Rda",sep="")
    save(bet,prop_rejs,file=FN)
    prop_rejs=NULL
  }
}
}

a=list.files('REVISION_RESULTS/outlier/')
a1=a[grepl("cl_",a)]
a2=a1[grepl("noderiv",a1)]
# a2F=a2[grepl("TRUE",a2)]
# #Should have 3*4*2*15*2
# a21=a2F[grepl("out_4",a2F)]
# 
# 
# 
# 
# file.remove(paste('REVISION_RESULTS/outlier/',a21[91:93],sep=""))
# 










##J=2
# ALPHA changes
########################################################### Alpha ###############################################################
###########################################################  J=2 ###############################################################
dist="N"
prop_rejs=NULL
count=0
tmp=seq(0.05,0.1,length.out=15)
for(N1 in c(50,100,250)){
  for(out in 1:4){
    for(both in c(T,F)){
    for(alphai in 1:15){
      
      alpha=tmp[alphai]
      
      

      fileName=paste0("REVISION_RESULTS/LTR/cl_",cl,"both_",both,"_out_",out,"_dist_",dist,"_LTR_depths_k_exp_alpha1_p05_",alphai,"_beta_p5_p5_N_",N1,"_num_runs_200.Rda",sep="")
      #REPLACE WITH REVISED DEPTHZ:
      load(fileName)
      depth_values1=depth_values
      
      fileName=paste0("REVISION_RESULTS/outlier/cl_",cl,"_both_",both,"_out_",out,"_dist_", dist,"_deriv_SIM_TD_k_exp_alpha_p05_",alphai,"_beta_p5_p5_N_",N1,"_num_runs_200.Rda",sep="")
      load(fileName)
      
      for(i in 1:num_runs){
        depth_values[[i]][,4]= depth_values1[[i]]
      }
      
      
      #KW
      test_stats=lapply(depth_values,getKWPM_all_depths,groups=c(1,N1+1),r=1,k=2,N=N1+N1)
      
      crit=qchisq(.95,1)
      
      prop_rej=colMeans(do.call(rbind,lapply(test_stats,'>',crit)))
      prop_rej
      
      
      ####PM
      test_stats=lapply(depth_values,getKWPM_all_depths,groups=c(1,N1+1),r=0.8,k=2,N=N1+N1)
      
      crit=qchisq(.95,1)
      
      prop_rej=c(prop_rej,colMeans(do.call(rbind,lapply(test_stats,'>',crit))))
      
      names(prop_rej)=c("FMD","RPD","MBD","LTR","RPD_Like","FMD_PM","RPD_PM","MBD_PM","LTR_PM","RPD_Like_PM")
      
      
      #NO DERIVATIVE
      
      fileName=paste0("REVISION_RESULTS/outlier/cl_",cl,"both_",both,"_out_",out,"_dist_",dist,"_noderiv_depths_k_exp_alpha1_p05_",alphai,"_beta_p5_p5_N_",N1,"_num_runs_200.Rda",sep="")
      load(fileName)
      depth_values1=depth_values
      
      
      
      
      #KW
      test_stats=lapply(depth_values,getKWPM_all_depths,groups=c(1,N1+1),r=1,k=2,N=N1+N1)
      
      crit=qchisq(.95,1)
      
      prop_rej2=colMeans(do.call(rbind,lapply(test_stats,'>',crit)))
      prop_rej=c(prop_rej, prop_rej2)
      prop_rej
      
      
      ####PM
      test_stats=lapply(depth_values,getKWPM_all_depths,groups=c(1,N1+1),r=0.8,k=2,N=N1+N1)
      
      crit=qchisq(.95,1)
      
      prop_rej=c(prop_rej,colMeans(do.call(rbind,lapply(test_stats,'>',crit))))
      
      names(prop_rej)=c(names(prop_rej)[1:10],
                        paste(c("FMD","RPD","MBD","LTR","FMD_PM","RPD_PM","MBD_PM","LTR_PM"),"_ND",sep=""))
      
      
      ##other tests, had to run guo seperately here
      fileName=paste0("REVISION_RESULTS/outlier/cl_",cl,"_both_",both,"_out_",out,"dist_",dist,"_other_k_exp_alpha1_p05_",alphai,"_beta_p5_p5_N_",N1,"_num_runs_200.Rda",sep="") 
      load(fileName)
      pvalues=do.call(rbind, p_values)
      
      
      prop_rej=c(prop_rej,apply(pvalues<0.05,2,mean))
      prop_rejs=rbind(prop_rejs,prop_rej)
      
      
    }
    alp=matrix(tmp,ncol=ncol(prop_rejs),nrow=nrow(prop_rejs))
    print(prop_rejs)
    FN=paste0("Power Files/cl_",cl,"_both_",both,"_out_",out,"dist_",dist,"_All_power_rev_k_exp_alpha_N_",N1,"J_2_num_runs_200.Rda",sep="")
    save(alp,prop_rejs,file=FN)
    prop_rejs=NULL
    }
  }
}








num_runs=200
library(xtable)
#Eigenvalue Decay
prop_rejs=NULL
scens=c("1","2","3","21","22","23")
for(N1 in c(50,100)){
  for(out in 1:4){
    for(both in c(T,F)){
    dist="N"
    for(scen in scens){
    
    # fileName=paste0("EIGENVALUE/allD_dist_N_deriv_SIM_TD_k_eigen_",scen,"_samesum_N_",N1,"_num_runs_200.Rda",sep="")
    # load(fileName)
    # depth_values1=depth_values
    
    fileName=paste0("EIGENVALUE/Outliers/cl_",cl,"_both_",both,"_out_",out,"_allD_dist_N_deriv_depths_k_eigen_",scen,"_samesum_N_",N1,"_num_runs_200.Rda",sep="")
    load(fileName)
    if(length(depth_values)<200)
      print(paste0("num_runs ",length(depth_values)))
    # nams=names(depth_values[[1]])
    # nams=c(nams,'RPD_Like')
    # for(i in 1:num_runs){
    #   depth_values[[i]][,4]=  depth_values1[[i]][,2]
    #   like=depth_values[[i]][,2]
    #   depth_values[[i]][,2]=  depth_values1[[i]][,1]
    #   depth_values[[i]]=cbind(depth_values[[i]],like)
    #   names(depth_values[[i]])=nams
    # }
    
    
    
    
    #KW
    test_stats=lapply(depth_values,getKWPM_all_depths,groups=c(1,N1+1),r=1,k=2,N=N1+N1)
    
    crit=qchisq(.95,1)
    
    prop_rej=colMeans(do.call(rbind,lapply(test_stats,'>',crit)));       prop_rej
    
    
    
    ####PM
    test_stats=lapply(depth_values,getKWPM_all_depths,groups=c(1,N1+1),r=0.8,k=2,N=N1+N1)
    
    crit=qchisq(.95,1)
    
    prop_rej=c(prop_rej,colMeans(do.call(rbind,lapply(test_stats,'>',crit))))
    
    names(prop_rej)=c("FMD","RPD","MBD","LTR","RPD_Like","FMD_PM","RPD_PM","MBD_PM","LTR_PM","RPD_Like_PM")
    
    
    fileName=paste0("EIGENVALUE/Outliers/cl_",cl,"_both_",both,"_out_",out,"_allD_dist_N_noderiv_depths_k_eigen_",scen,"_samesum_N_",N1,"_num_runs_200.Rda",sep="")
    load(fileName)
    if(length(depth_values)<200)
      print(paste0("num_runs ",length(depth_values)))
    
    
    
    #KW
    test_stats=lapply(depth_values,getKWPM_all_depths,groups=c(1,N1+1),r=1,k=2,N=N1+N1)
    
    crit=qchisq(.95,1)
    
    prop_rej=c(prop_rej,colMeans(do.call(rbind,lapply(test_stats,'>',crit))));       prop_rej
    
    
    
    ####PM
    test_stats=lapply(depth_values,getKWPM_all_depths,groups=c(1,N1+1),r=0.8,k=2,N=N1+N1)
    
    crit=qchisq(.95,1)
    
    prop_rej=c(prop_rej,colMeans(do.call(rbind,lapply(test_stats,'>',crit))))
    
    names(prop_rej)=c(names(prop_rej)[1:10],
                      paste(c("FMD","RPD","MBD","LTR","FMD_PM","RPD_PM","MBD_PM","LTR_PM"),"_ND",sep=""))
    
    
    
    
    fileName=paste0("EIGENVALUE/other/cl_",cl,"_both_",both,"_out_",out,"_other_dist_N_k_eigen_",scen,"_samesum_N_",N1,"_num_runs_200.Rda",sep="")
    a=try(load(fileName),T)
    errorsp=inherits(a, "try-error")
    if(errorsp){
      fileName=paste0("EIGENVALUE/other/cl_",cl,"_both_",both,"_out_",out,"_other_dist_N_depths_k_eigen_",scen,"_samesum_N_",N1,"_num_runs_200.Rda",sep="")
      load(fileName)
    }
    # pvalues=do.call(rbind, p_values)
    pvalues=t(do.call(rbind,p_values))
    if(dim(pvalues)[2]<200)
      print(paste0(fileName," num_runs ",dim(p_values)[2]))
    # pvalues=p_values
    
    
    
    
    prop_rej=c(prop_rej,apply(pvalues<0.05,1,mean))
    # names(prop_rej)[28:34]=rownames( p_values)
    prop_rejs=rbind(prop_rejs,prop_rej)
    
    
  }
  print(prop_rejs)
  FN=paste0("Power Files/cl_",cl,"_both_",both,"_out_",out,"eigen_",dist,"_All_power_rev_N_",N1,"_num_runs_200.Rda",sep="")
  xtable(t(prop_rejs))
  save(scens,prop_rejs,file=FN)
  prop_rejs=NULL
  }
  }
  }
beepr::beep()


############# Additional contamination levels




for(contamination_level in c(0.01,0.025,0.05)){
  
  cl=contamination_level*100
  
  dist='N'
  print("sim 1")
  #changes in beta
  prop_rejs=NULL
  count=0
  tmp=seq(0.5,1,length.out = 15)
  N1=100
  N2=N1
    for(out in 1:4){
      for(both in c(T,F)){
        for(betai in 1:15){
          
          

          
          fileName=paste0("REVISION_RESULTS/outlier/cl_",cl,"_both_",both,"_out_",out,"_dist_", dist,"_deriv_SIM_TD_k_exp_alpha1_p05_p05_beta_p", betai, "_p5_N_",N1,"_num_runs_200.Rda", sep = "")
          load(fileName)
          
          
  
          
          
          
          #KW
          test_stats=lapply(depth_values,getKWPM_all_depths,groups=c(1,N1+1),r=1,k=2,N=N1+N2)
          
          crit=qchisq(.95,1)
          
          prop_rej=colMeans(do.call(rbind,lapply(test_stats,'>',crit)));       prop_rej
          
          
          
          ####PM
          test_stats=lapply(depth_values,getKWPM_all_depths,groups=c(1,N1+1),r=0.8,k=2,N=N1+N2)
          
          crit=qchisq(.95,1)
          
          prop_rej=c(prop_rej,colMeans(do.call(rbind,lapply(test_stats,'>',crit))))
          
          # names(prop_rej)=c("FMD","RPD","MBD","LTR","SPAT","KSPAT","FMD_PM","RPD_PM","MBD_PM","LTR","SPAT_PM","KSPAT_PM")
          names(prop_rej)=c("FMD","RPD","MBD","LTR","RPD_Like","FMD_PM","RPD_PM","MBD_PM","LTR_PM","RPD_Like_PM")
          
          
          
   
          
          
          ##other tests, had to run guo seperately here
          
          fileName=paste0("REVISION_RESULTS/outlier/cl_",cl,"_both_",both,"_out_",out,"dist_",dist,"_other_k_exp_alpha1_p05_p05_beta_p",betai,"_p5_N_",N1,"_num_runs_200.Rda",sep="")      
          a=try({load(fileName)},T)
          errorsp=inherits(a, "try-error")
          if(errorsp){
            fileName=paste0("REVISION_RESULTS/outlier/cl_",cl,"both_",both,"_out_",out,"_dist_",dist,"_other_k_exp_alpha1_p05_p05_beta_p",betai,"_p5_N_",N1,"_num_runs_200.Rda",sep="")      
            load(fileName)
          }
          pvalues=do.call(rbind, p_values)
          
          
          prop_rej=c(prop_rej,apply(pvalues<0.05,2,mean))
          prop_rejs=rbind(prop_rejs,prop_rej)
          
          
        }
        bet=matrix(tmp,ncol=ncol(prop_rejs),nrow=nrow(prop_rejs))
        print(prop_rejs)
        FN=paste0("Power Files/cl_",cl,"_both_",both,"_out_",out,"dist_",dist,"_All_power_rev_k_exp_beta_N_",N1,"J_2_num_runs_200.Rda",sep="")
        save(bet,prop_rejs,file=FN)
        prop_rejs=NULL
      }
    }
  
  
  a=list.files('REVISION_RESULTS/outlier/')
  a1=a[grepl("cl_",a)]
  a2=a1[grepl("noderiv",a1)]
  # a2F=a2[grepl("TRUE",a2)]
  # #Should have 3*4*2*15*2
  # a21=a2F[grepl("out_4",a2F)]
  # 
  # 
  # 
  # 
  # file.remove(paste('REVISION_RESULTS/outlier/',a21[91:93],sep=""))
  # 
  

  ##J=2
  # ALPHA changes
  ########################################################### Alpha ###############################################################
  ###########################################################  J=2 ###############################################################
  dist="N"
  prop_rejs=NULL
  count=0
  tmp=seq(0.05,0.1,length.out=15)
  
    for(out in 1:4){
      for(both in c(T,F)){
        for(alphai in 1:15){
          
          alpha=tmp[alphai]
          
          
          

          
          fileName=paste0("REVISION_RESULTS/outlier/cl_",cl,"_both_",both,"_out_",out,"_dist_", dist,"_deriv_SIM_TD_k_exp_alpha_p05_",alphai,"_beta_p5_p5_N_",N1,"_num_runs_200.Rda",sep="")
          load(fileName)
          

          
          #KW
          test_stats=lapply(depth_values,getKWPM_all_depths,groups=c(1,N1+1),r=1,k=2,N=N1+N1)
          
          crit=qchisq(.95,1)
          
          prop_rej=colMeans(do.call(rbind,lapply(test_stats,'>',crit)))
          prop_rej
          
          
          ####PM
          test_stats=lapply(depth_values,getKWPM_all_depths,groups=c(1,N1+1),r=0.8,k=2,N=N1+N1)
          
          crit=qchisq(.95,1)
          
          prop_rej=c(prop_rej,colMeans(do.call(rbind,lapply(test_stats,'>',crit))))
          
          names(prop_rej)=c("FMD","RPD","MBD","LTR","RPD_Like","FMD_PM","RPD_PM","MBD_PM","LTR_PM","RPD_Like_PM")
          
          
         
          ##other tests, had to run guo seperately here
          fileName=paste0("REVISION_RESULTS/outlier/cl_",cl,"_both_",both,"_out_",out,"dist_",dist,"_other_k_exp_alpha1_p05_",alphai,"_beta_p5_p5_N_",N1,"_num_runs_200.Rda",sep="") 
          load(fileName)
          pvalues=do.call(rbind, p_values)
          
          
          prop_rej=c(prop_rej,apply(pvalues<0.05,2,mean))
          prop_rejs=rbind(prop_rejs,prop_rej)
          
          
        }
        alp=matrix(tmp,ncol=ncol(prop_rejs),nrow=nrow(prop_rejs))
        print(prop_rejs)
        FN=paste0("Power Files/cl_",cl,"_both_",both,"_out_",out,"dist_",dist,"_All_power_rev_k_exp_alpha_N_",N1,"J_2_num_runs_200.Rda",sep="")
        save(alp,prop_rejs,file=FN)
        prop_rejs=NULL
      }
    }
  
  
  
  
  

}



