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



#J=2 
# Beta Changes


print("sim 1")
#changes in beta
prop_rejs=NULL
count=0
tmp=seq(0.5,1,length.out = 15)
N1=100
  for(dist in c("N","t")){
    for(betai in 1:15){
      
      
      #for N1=100 I ran FM seperately
      beta=tmp[betai]
      
      #REPLACE WITH REVISED DEPTHZ:
      
      fileName=paste0("REVISION_RESULTS/missing/dist_",dist,"_missing_k_exp_alpha1_p05_p05_beta_p",betai,"_p5_N_",N1,"_num_runs_200.Rda",sep="")
      load(fileName)

      
      
      
      #KW
      # test_stats=lapply(depth_values,getKW_all_depths,N1=N1)
      test_stats=lapply(depth_values,getKWPM_all_depths,groups=c(1,N1+1),r=1,k=2,N=N1+N2)
      crit=qchisq(.95,1)
      pr=colMeans(do.call(rbind,lapply(test_stats,'>',crit))); pr
      prop_rej=colMeans(do.call(rbind,lapply(test_stats,'>',crit)));       prop_rej
      
      
      names(prop_rej)=c("FMD","RPD","MBD","LTR","RPD_Like")
      # names(prop_rej)=c("FMD","RPD","MBD","LTR")
      

      ##other tests, had to run guo seperately here
      
      fileName=paste0("REVISION_RESULTS/missing/missing_dist_",dist,"_other_k_exp_alpha1_p05_p05_beta_p",betai,"_p5_N_",N1,"_num_runs_200.Rda",sep="")      
      load(fileName)
      pvalues=do.call(rbind, p_values)
      
      
      prop_rej=c(prop_rej,apply(pvalues<0.05,2,mean))
      prop_rejs=rbind(prop_rejs,prop_rej)
      
      
    }
    bet=matrix(tmp,ncol=ncol(prop_rejs),nrow=nrow(prop_rejs))
    print(prop_rejs)
    FN=paste0("Power Files/missing_dist_",dist,"_All_power_rev_k_exp_beta_N_",N1,"J_2_num_runs_200.Rda",sep="")
    save(bet,prop_rejs,file=FN)
    prop_rejs=NULL
  }



##J=2
# ALPHA changes

prop_rejs=NULL
count=0
tmp=seq(0.05,0.1,length.out=15)
N1=100
  for(dist in c("N","t")){
    for(alphai in 1:15){
      N2=N1
      
      alpha=tmp[alphai]
      
      fileName=paste0("REVISION_RESULTS/missing/dist_",dist,"_missing_k_exp_alpha_p05_",alphai,"_beta_p5_p5_N_",N1,"_num_runs_200.Rda",sep="")
      load(fileName)
    
      
      #KW
      # test_stats=lapply(depth_values,getKW_all_depths,N1=N1)
      test_stats=lapply(depth_values,getKWPM_all_depths,groups=c(1,N1+1),r=1,k=2,N=N1+N2)
      
      crit=qchisq(.95,1)
      
      prop_rej=colMeans(do.call(rbind,lapply(test_stats,'>',crit)))
      prop_rej

      
      names(prop_rej)=c("FMD","RPD","MBD","LTR","RPD_Like")
      
      

      
      
      ##other tests, had to run guo seperately here
      fileName=paste0("REVISION_RESULTS/missing/missing_dist_",dist,"_other_k_exp_alpha1_p05_",alphai,"_beta_p5_p5_N_",N1,"_num_runs_200.Rda",sep="")      
      load(fileName)
      pvalues=do.call(rbind, p_values)
      
      
      prop_rej=c(prop_rej,apply(pvalues<0.05,2,mean))
      prop_rejs=rbind(prop_rejs,prop_rej)
      
      
    }
    alp=matrix(tmp,ncol=ncol(prop_rejs),nrow=nrow(prop_rejs))
    print(prop_rejs)
    FN=paste0("Power Files/missing_dist_",dist,"_All_power_rev_k_exp_alpha_N_",N1,"J_2_num_runs_200.Rda",sep="")
    save(alp,prop_rejs,file=FN)
    prop_rejs=NULL
  }







