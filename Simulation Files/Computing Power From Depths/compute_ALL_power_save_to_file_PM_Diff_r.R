
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
  groups=groups[-(k+1)]
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


get_all_results=function(r,N,k=2){
  test_stats=lapply(depth_values,getKWPM_all_depths,groups=groups,r=r,k=k,N=N)
  crit=qchisq(.95,k-1)
  return(colMeans(do.call(rbind,lapply(test_stats,'>',crit))))
}



setwd("Functional Data Covariance Files/Nov_2020_sim_results")

rs=c(0.65,0.7,0.75,0.8,0.85,0.9,0.95)
num_runs=200


######### ######### ######### #########  changes in beta ######### ######### ######### #########  c######### ######### ######### #########  c
tmp=seq(0.5,1,length.out = 15)
print("sim 1")
for(N in c(50,100,250)*2){
  for(ratio in c(0.2,0.3,0.4)){
    N1=floor(ratio*N)
    N2=N-N1
    groups=c(1,N1+1,N1+N2+1)
    k=2
    N=N1+N2
    for(dist in c("N","t","sn")){
      power=array(0,dim=c(5,length(rs),length(c(1,7,15))))
      rownames(power)=c("FMD","RPD","MBD","LTR","RPD_Like")
      for(betai in c(1,7,15)){
        
        
        #for N1=100 I ran FM separately
        beta=tmp[betai]
        
        if(betai!=1){
          fileName=paste0("REVISION_RESULTS/LTR/dist_",dist,"_deriv_LTR_k_exp_alpha1_p05_p05_beta_p",betai,"_p5_N1_",N1,"_N2_",N2,"_num_runs_200.Rda",sep="")
          load(fileName)
          depth_values1=depth_values
        }
        
        fileName=paste0("REVISION_RESULTS/DGS/dist_",dist,"_deriv_SIM_TD_k_exp_alpha1_p05_p05_beta_p",betai,"_p5_N1_",N1,"_N2_",N2,"_num_runs_200.Rda",sep="")
  
        load(fileName)
        
        
        if(betai!=1){
          for(i in 1:num_runs){
            depth_values[[i]][,4]= depth_values1[[i]]
          }
          
        }
        
        
        
        
        
        ####PM
        power[,,which(betai==c(1,7,15))]=sapply(rs,get_all_results,N)
        
        
        
        
        
        
        
      }
      bet=matrix(tmp[c(1,7,15)],ncol=dim(power)[1],nrow=dim(power)[3])
      print(power)
      FN=paste0("Power Files/dist_",dist,"_All_power_rs_beta_N1_",N1,"_N2_",N2,"J_2_num_runs_200.Rda",sep="")
      save(bet,power,file=FN)

    }
  }
}




######### ######### ######### #########  changes in alpha ######### ######### ######### #########  c######### ######### ######### #########  c
tmp=seq(0.05,0.1,length.out=15)
print("sim 1")
for(N in c(50,100,250)*2){
  for(ratio in c(0.2,0.3,0.4)){
    N1=floor(ratio*N)
    N2=N-N1
    groups=c(1,N1+1,N1+N2+1)
    k=2
    N=N1+N2
    for(dist in c("N","t","sn")){
      power=array(0,dim=c(5,length(rs),length(c(1,7,15))))
      rownames(power)=c("FMD","RPD","MBD","LTR","RPD_Like")
      for(alphai in c(1,7,15)){
        
        
        #for N1=100 I ran FM separately
        alpha=tmp[alphai]
        
        if(alphai!=1){
          fileName=paste0("REVISION_RESULTS/LTR/dist_",dist,"_deriv_LTR_k_exp_alpha_p05_",alphai,"_beta_p5_p5_N1_",N1,"_N2_",N2,"_num_runs_200.Rda",sep="")
          load(fileName)
          depth_values1=depth_values
        }
        
        fileName=paste0("REVISION_RESULTS/DGS/dist_",dist,"_deriv_SIM_TD_k_exp_alpha_p05_",alphai,"_beta_p5_p5_N1_",N1,"_N2_",N2,"_num_runs_200.Rda",sep="")
        load(fileName)
        
        if(alphai!=1){
          for(i in 1:num_runs){
            depth_values[[i]][,4]= depth_values1[[i]]
          }
          
        }
        
        
        
        
        
        ####PM
        power[,,which(alphai==c(1,7,15))]=sapply(rs,get_all_results,N)
        
        
        
        
        
        
        
      }
      alp=matrix(tmp[c(1,7,15)],ncol=dim(power)[1],nrow=dim(power)[3])
      print(power)
      FN=paste0("Power Files/dist_",dist,"_All_power_rs_alpha_N1_",N1,"_N2_",N2,"J_2_num_runs_200.Rda",sep="")
      save(alp,power,file=FN)
      
    }
  }
}









######### ######### ######### #########  changes in beta j is 3 ######### ######### ######### #########  c######### ######### ######### #########  




tmp=seq(0.5,1,length.out = 15)
print("sim 1")
for(N1 in c(50,100)){
  for(ratio in c(0.4,0.6,0.8)){
    N2=floor(N1*ratio)
    N3=N1
    groups=c(1,N1+1,N1+N2+1,N1+N2+N3+1)
    k=3
    N=N1+N2+N3
    for(dist in c("N","t","sn")){
      power=array(0,dim=c(5,length(rs),length(c(1,7,15))))
      rownames(power)=c("FMD","RPD","MBD","LTR","RPD_Like")
      for(betai in c(1,7,15)){
        
        
        #for N1=100 I ran FM separately
        beta=tmp[betai]
        
        if(betai!=1){
          fileName=paste0("REVISION_RESULTS/LTR/dist_",dist,"_deriv_LTR_k_exp_alpha1_p05_p05_p05_beta_p5_p",betai,"_p5_N1_",N1,"_N2_",N2,"_N3_",N3,"_num_runs_200.Rda",sep="")
          load(fileName)
          depth_values1=depth_values
        }
        
        fileName=paste0("REVISION_RESULTS/DGS/dist_",dist,"_deriv_depths_k_exp_alpha1_p05_p05_p05_beta_p5_p",betai,"_p5_N1_",N1,"_N2_",N2,"_N3_",N3,"_num_runs_200.Rda",sep="")
        
        load(fileName)
        
        
        if(betai!=1){
          for(i in 1:num_runs){
            depth_values[[i]][,4]= depth_values1[[i]]
          }
          
        }
        
        
        
        
        
        ####PM
        power[,,which(betai==c(1,7,15))]=sapply(rs,get_all_results,N,k=k)
        
        
        
        
        
        
        
      }
      bet=matrix(tmp[c(1,7,15)],ncol=dim(power)[1],nrow=dim(power)[3])
      print(power)
      FN=paste0("Power Files/dist_",dist,"_All_power_rs_beta_N1_",N1,"_N2_",N2,"_N3_",N3,"_J_3_num_runs_200.Rda",sep="")
      save(bet,power,file=FN)
      
    }
  }
}









num_runs=200
print("sim 2")
#changes in alpha
prop_rejs=NULL
count=0
tmp=seq(0.05,0.1,length.out=15)
for(N1 in c(50,100)){
  for(ratio in c(0.4,0.6,0.8)){
    N2=floor(N1*ratio)
    N3=N1
    groups=c(1,N1+1,N1+N2+1,N1+N2+N3+1)
    k=3
    N=N1+N2+N3
    for(dist in c("N","t","sn")){
      for(alphai in c(1,7,15)){
        alpha=tmp[alphai]
        
        
        if(alphai!=1){
          fileName=paste0("REVISION_RESULTS/LTR/dist_",dist,"_deriv_LTR_k_exp_alpha1_p05_",alphai,"_p05_beta_p5_p5_p5_N1_",N1,"_N2_",N2,"_N3_",N3,"_num_runs_200.Rda",sep="")
          load(fileName)
          depth_values1=depth_values
        }
        
        fileName=paste0("REVISION_RESULTS/DGS/dist_",dist,"_deriv_depths_k_exp_alpha1_p05_",alphai,"_p05_beta_p5_p5_p5_N1_",N1,"_N2_",N2,"_N3_",N3,"_num_runs_200.Rda",sep="")
        load(fileName)
        if(alphai!=1){
          for(i in 1:num_runs){
            depth_values[[i]][,4]= depth_values1[[i]]
          }
        }
        
        
        
        #KW
        test_stats=lapply(depth_values,getKW_all_depths,groups,k,N)
        
        crit=qchisq(.95,2)
        
        prop_rej=colMeans(do.call(rbind,lapply(test_stats,'>',crit)))
        prop_rej
        
        
        ####PM
        test_stats=lapply(depth_values,getKWPM_all_depths,groups,r=0.8,k,N)
        
        crit=qchisq(.95,2)
        
        prop_rej=c(prop_rej,colMeans(do.call(rbind,lapply(test_stats,'>',crit))))
        
        names(prop_rej)=c("FMD","RPD","MBD","LTR","RPD_Like","FMD_PM","RPD_PM","MBD_PM","LTR_PM","RPD_Like_PM")
        prop_rej
        
        #NO DERIVATIVE
        
        fileName=paste0("REVISION_RESULTS/DGS/dist_",dist,"_noderiv_depths_k_exp_alpha1_p05_",alphai,"_p05_beta_p5_p5_p5_N1_",N1,"_N2_",N2,"_N3_",N3,"_num_runs_200.Rda",sep="")
        load(fileName)
        depth_values1=depth_values
        
        
        
        
        #KW
        test_stats=lapply(depth_values,getKW_all_depths,groups,k,N)
        
        crit=qchisq(.95,2)
        
        prop_rej2=colMeans(do.call(rbind,lapply(test_stats,'>',crit)))
        prop_rej=c(prop_rej, prop_rej2)
        prop_rej
        
        
        ####PM
        test_stats=lapply(depth_values,getKWPM_all_depths,groups,r=0.8,k,N)
        
        crit=qchisq(.95,2)
        
        prop_rej=c(prop_rej,colMeans(do.call(rbind,lapply(test_stats,'>',crit))))
        
        names(prop_rej)=c(names(prop_rej)[1:10],
                          paste(c("FMD","RPD","MBD","LTR","FMD_PM","RPD_PM","MBD_PM","LTR_PM"),"_ND",sep=""))
        
        
        ##other tests, had to run guo seperately here
        fileName=paste0("REVISION_RESULTS/DGS/dist_",dist,"_other_k_exp_alpha1_p05_",alphai,"_p05_beta_p5_p5_p5_N1_",N1,"_N2_",N2,"_N3_",N3,"_num_runs_200.Rda",sep="")
        load(fileName)
        pvalues=do.call(rbind, p_values)
        
        
        prop_rej=c(prop_rej,apply(pvalues<0.05,2,mean))
        prop_rejs=rbind(prop_rejs,prop_rej)
        
        
      }
      alp=matrix(tmp[c(1,7,15)],ncol=ncol(prop_rejs),nrow=nrow(prop_rejs))
      print(prop_rejs)
      FN=paste0("Power Files/dist_",dist,"_All_power_rev_k_exp_alpha_N1_",N1,"_N2_",N2,"_N3_",N3,"J_3_num_runs_200.Rda",sep="")
      save(alp,prop_rejs,file=FN)
      prop_rejs=NULL
    }
  }
}




beepr::beep()






num_runs=200
library(xtable)
#Eigenvalue Decay
prop_rejs=NULL
scens=c("1","2","3","21","22","23")
for(N in c(50,100)*2){
  for(ratio in c(0.2,0.3,0.4)){
    N1=floor(ratio*N)
    N2=N-N1
    groups=c(1,N1+1,N1+N2+1)
    k=2
    N=N1+N2
    dist="N"
    for(scen in scens){
      
      # fileName=paste0("EIGENVALUE/allD_dist_N_deriv_SIM_TD_k_eigen_",scen,"_samesum_N1_",N1,"_N2_",N2,"_num_runs_200.Rda",sep="")
      # load(fileName)
      # depth_values1=depth_values
      
      fileName=paste0("EIGENVALUE/DGS/allD_dist_N_deriv_depths_k_eigen_",scen,"_samesum_N1_",N1,"_N2_",N2,"_num_runs_200.Rda",sep="")
      load(fileName)
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
      test_stats=lapply(depth_values,getKW_all_depths,groups,k,N)
      
      crit=qchisq(.95,1)
      
      prop_rej=colMeans(do.call(rbind,lapply(test_stats,'>',crit)));       prop_rej
      
      
      
      ####PM
      test_stats=lapply(depth_values,getKWPM_all_depths,groups=groups,r=0.8,k=2,N=N)
      
      crit=qchisq(.95,1)
      
      prop_rej=c(prop_rej,colMeans(do.call(rbind,lapply(test_stats,'>',crit))))
      
      names(prop_rej)=c("FMD","RPD","MBD","LTR","RPD_Like","FMD_PM","RPD_PM","MBD_PM","LTR_PM","RPD_Like_PM")
      
      
      fileName=paste0("EIGENVALUE/DGS/allD_dist_N_noderiv_depths_k_eigen_",scen,"_samesum_N1_",N1,"_N2_",N2,"_num_runs_200.Rda",sep="")
      load(fileName)
      
      
      
      
      #KW
      test_stats=lapply(depth_values,getKW_all_depths,groups,k,N)
      
      crit=qchisq(.95,1)
      
      prop_rej=c(prop_rej,colMeans(do.call(rbind,lapply(test_stats,'>',crit))));       prop_rej
      
      
      
      ####PM
      test_stats=lapply(depth_values,getKWPM_all_depths,groups=groups,r=0.8,k=2,N=N)
      
      crit=qchisq(.95,1)
      
      prop_rej=c(prop_rej,colMeans(do.call(rbind,lapply(test_stats,'>',crit))))
      
      names(prop_rej)=c(names(prop_rej)[1:10],
                        paste(c("FMD","RPD","MBD","LTR","FMD_PM","RPD_PM","MBD_PM","LTR_PM"),"_ND",sep=""))
      
      
      
      
      fileName=paste0("EIGENVALUE/other/other_dist_N_k_eigen_",scen,"_samesum__N1_",N1,"_N2_",N2,"_num_runs_200.Rda",sep="")
      load(fileName)
      # pvalues=do.call(rbind, p_values)
      pvalues=p_values
      
      prop_rej=c(prop_rej,apply(pvalues<0.05,1,mean))
      # names(prop_rej)[28:34]=rownames( p_values)
      prop_rejs=rbind(prop_rejs,prop_rej)
      
      
    }
    print(prop_rejs)
    FN=paste0("Power Files/eigen_",dist,"_All_power_rev_N1_",N1,"_N2_",N2,"_num_runs_200.Rda",sep="")
    xtable(t(prop_rejs))
    save(scens,prop_rejs,file=FN)
    prop_rejs=NULL
  }
}

beepr::beep()








# setwd("EIGENVALUE/other")
# files=list.files()
# for(file in files){
#   
#   file.rename(file, gsub('depths_','',file))
#   
# }







