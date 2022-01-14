setwd("Functional Data Covariance Files")
library(BBmisc)
##J=2

setwd("")


# files=list.files(path="Two_Sample_Depths_Cov_Kernels\\RPD_MBD\\Spatial No Deriv",pattern="N_250")
# files
# files=list.files(pattern="N_250")
# files
# files=list.files(pattern="spat2")
# files2=gsub("spat2","spat",files)

is.error=function(ob){class(ob)=="try-error"}



#SIM 1
tmp=seq(0.5,1,length.out = 15)
for(N1 in c(50,100,250)){
  for(dist in c("N","t","sn")){
    for(betai in 1:15){
      
      beta=tmp[betai]
      fileName=paste0("FKWCResults/RPDMBD/dist_",dist,"_deriv_depths_k_exp_alpha1_p05_p05_beta_p",betai,"_p5_N_",N1,"_num_runs_200.Rda",sep="")
      att1=try(load(fileName),T)

      fileName=paste0("FKWCResults/RPDMBD/dist_",dist,"_noderiv_depths_k_exp_alpha1_p05_p05_beta_p",betai,"_p5_N_",N1,"_num_runs_200.Rda",sep="")
      att7=try(load(fileName),T)
      
      #load spatial and k spatial
      if(betai<=6){
        fileName=paste0("FKWCResults/OSPAT_D/ospat_dist_",dist,"_deriv_depths_k_exp_alpha1_p05_p05_beta_p",betai,"_p5_N_",N1,"_num_runs_50.Rda",sep="")   
        # fileName=paste0("FKWCResults/ospat_dist_",dist,"_deriv_depths_k_exp_alpha1_p05_p05_beta_p",betai,"_p5_N_",N1,"_num_runs_50.Rda",sep="")   
        att2=try(load(fileName),T)
      
        fileName=paste0("FKWCResults/OSPAT_ND/o_spat_dist_",dist,"_noderiv_depths_k_exp_alpha1_p05_p05_beta_p",betai,"_p5_N_",N1,"_num_runs_50.Rda",sep="")   
        # fileName=paste0("FKWCResults/o_spat_dist_",dist,"_noderiv_depths_k_exp_alpha1_p05_p05_beta_p",betai,"_p5_N_",N1,"_num_runs_50.Rda",sep="")   
        att3=try(load(fileName),T)
      

        fileName=paste0("FKWCResults/KSPAT_D/kspat_dist_",dist,"_deriv_depths_k_exp_alpha1_p05_p05_beta_p",betai,"_p5_N_",N1,"_num_runs_50.Rda",sep="")   
        # fileName=paste0("FKWCResults/kspat_dist_",dist,"_deriv_depths_k_exp_alpha1_p05_p05_beta_p",betai,"_p5_N_",N1,"_num_runs_50.Rda",sep="")   
        
        att4=try(load(fileName),T)
         
        
        fileName=paste0("FKWCResults/KSPAT_ND/kspat_dist_",dist,"_noderiv_depths_k_exp_alpha1_p05_p05_beta_p",betai,"_p5_N_",N1,"_num_runs_50.Rda",sep="")   
        # fileName=paste0("FKWCResults/kspat_dist_",dist,"_noderiv_depths_k_exp_alpha1_p05_p05_beta_p",betai,"_p5_N_",N1,"_num_runs_50.Rda",sep="")   
        
        att5=try(load(fileName),T)
      }
      else{
        att2=att3=att4=att5="d"
      }
      #no derivative
    
      ##other tests, had to run guo seperately here

      fileName=paste0("OTHER_TEST_RESULTS/dist_",dist,"_other_k_exp_alpha1_p05_p05_beta_p",betai,"_p5_N_",N1,"_num_runs_200.Rda",sep="")      
      att6=try(load(fileName),T)
      

      if(any(unlist(lapply(list(att1,att2,att3,att4,att5,att6,att7),is.error)))){
        print("***********************************")
        print(paste0("beta ",beta,"dist ",dist," N ",N1))
        print(which(unlist(lapply(list(att1,att2,att3,att4,att5,att6,att7),is.error))))
        print("***********************************")
        print("***********************************")
      }
    }
  }
}




#SIM 2
tmp=seq(0.05,0.1,length.out=15)
for(N1 in c(50,100,250)){
  for(dist in c("N","t","sn")){
    for(alphai in 2:length(tmp)){
      alpha=tmp[alphai]
       N2=N1
       fileName=paste0("FKWCResults/RPDMBD/dist_",dist,"_deriv_depths_k_exp_alpha_p05_",alphai,"_beta_p5_p5_N_",N1,"_num_runs_200.Rda",sep="")
       att1=try(load(fileName),T)
       fileName=paste0("FKWCResults/RPDMBD/dist_",dist,"_noderiv_depths_k_exp_alpha1_p05_",alphai,"_beta_p5_p5_N_",N1,"_num_runs_200.Rda",sep="")
       att7=try(load(fileName),T)
       
       
       
       fileName=paste0("OTHER_TEST_RESULTS/dist_",dist,"_other_k_exp_alpha1_p05_",alphai,"_beta_p5_p5_N_",N1,"_num_runs_200.Rda",sep="")      
       att6=try(load(fileName),T)
      
      if(any(unlist(lapply(list(att1,att2,att3,att4,att5,att6,att7),is.error)))){
        print("**********************************")
        print(paste0("alpha ",alpha,"dist ",dist,"N ",N1))
        print(which(unlist(lapply(list(att1,att2,att3,att4,att5,att6,att7),is.error))))
        print("***********************************")
        print("***********************************")
      }
    }
  }
}



##J=3


#SIM 1
tmp=seq(0.5,1,length.out = 15)
for(N1 in c(50,100)){
  for(dist in c("N","t","sn")){
    for(beta in 1:length(tmp)){
      beta=tmp[betai]
      fileName=paste0("FKWCResults/RPDMBD/dist_",dist,"_deriv_depths_k_exp_alpha1_p05_p05_p05_beta_p5_p",betai,"_p5_N_",N1,"_num_runs_200.Rda",sep="")
      att1=try(load(fileName),T)
      
      fileName=paste0("FKWCResults/RPDMBD/dist_",dist,"_noderiv_depths_k_exp_alpha1_p05_p05_p05_beta_p5_p",betai,"_p5_N_",N1,"_num_runs_200.Rda",sep="")
      att7=try(load(fileName),T)
      
      #load spatial and k spatial
      if(betai<=6){
        fileName=paste0("FKWCResults/OSPAT_D/ospat_dist_",dist,"_deriv_depths_k_exp_alpha1_p05_p05_beta_p",betai,"_p5_N_",N1,"_num_runs_50.Rda",sep="")   
        # fileName=paste0("FKWCResults/ospat_dist_",dist,"_deriv_depths_k_exp_alpha1_p05_p05_beta_p",betai,"_p5_N_",N1,"_num_runs_50.Rda",sep="")   
        att2=try(load(fileName),T)
        
        fileName=paste0("FKWCResults/OSPAT_ND/o_spat_dist_",dist,"_noderiv_depths_k_exp_alpha1_p05_p05_beta_p",betai,"_p5_N_",N1,"_num_runs_50.Rda",sep="")   
        # fileName=paste0("FKWCResults/o_spat_dist_",dist,"_noderiv_depths_k_exp_alpha1_p05_p05_beta_p",betai,"_p5_N_",N1,"_num_runs_50.Rda",sep="")   
        att3=try(load(fileName),T)
        
        
        fileName=paste0("FKWCResults/KSPAT_D/kspat_dist_",dist,"_deriv_depths_k_exp_alpha1_p05_p05_beta_p",betai,"_p5_N_",N1,"_num_runs_50.Rda",sep="")   
        # fileName=paste0("FKWCResults/kspat_dist_",dist,"_deriv_depths_k_exp_alpha1_p05_p05_beta_p",betai,"_p5_N_",N1,"_num_runs_50.Rda",sep="")   
        
        att4=try(load(fileName),T)
        
        
        fileName=paste0("FKWCResults/KSPAT_ND/kspat_dist_",dist,"_noderiv_depths_k_exp_alpha1_p05_p05_beta_p",betai,"_p5_N_",N1,"_num_runs_50.Rda",sep="")   
        # fileName=paste0("FKWCResults/kspat_dist_",dist,"_noderiv_depths_k_exp_alpha1_p05_p05_beta_p",betai,"_p5_N_",N1,"_num_runs_50.Rda",sep="")   
        
        att5=try(load(fileName),T)
        att2=att3=att4=att5="d"
      }
      else{
        att2=att3=att4=att5="d"
      }
      #no derivative
      
      ##other tests, had to run guo seperately here
      fileName=paste0("OTHER_TEST_RESULTS/dist_",dist,"_other_k_exp_alpha1_p05_p05_p05_beta_p5_p",betai,"_p5_N_",N1,"_num_runs_200.Rda",sep="")
      att6=try(load(fileName),T)
    
      
      if(any(unlist(lapply(list(att1,att2,att3,att4,att5,att6,att7),is.error)))){
        print("***********************************")
        print(paste0("beta ",beta,"dist ",dist," N ",N1))
        print(which(unlist(lapply(list(att1,att2,att3,att4,att5,att6,att7),is.error))))
        print("***********************************")
        print("***********************************")
      }
    }
  }
}


#SIM 2
tmp=seq(0.05,0.1,length.out=15)
for(N1 in c(50,100)){
  for(dist in c("N","t","sn")){
    for(alpha in tmp){
      
      N2=N1
      alpha=tmp[alphai]
      N2=N1
      fileName=paste0("FKWCResults/RPDMBD/dist_",dist,"_deriv_depths_k_exp_alpha1_p05_",alphai,"_p05_beta_p5_p5_p5_N_",N1,"_num_runs_200.Rda",sep="")
      att1=try(load(fileName),T)
      fileName=paste0("FKWCResults/RPDMBD/dist_",dist,"_noderiv_depths_k_exp_alpha1_p05_",alphai,"_p05_beta_p5_p5_p5_N_",N1,"_num_runs_200.Rda",sep="")
      att7=try(load(fileName),T)
      
      
      if(alphai<=6){
        #load spatial and k spatial
        fileName=paste0("FKWCResults/OSPAT_D/ospat_dist_",dist,"_deriv_depths_k_exp_alpha_p05_",alphai,"_beta_p5_p5_N_",N1,"_num_runs_50.Rda",sep="")   
        att2=try(load(fileName),T)
        
        fileName=paste0("FKWCResults/OSPAT_ND/o_spat_dist_",dist,"_noderiv_depths_k_exp_alpha1_p05_",alphai,"_beta_p5_p5_N_",N1,"_num_runs_50.Rda",sep="")   
        att3=try(load(fileName),T)
        
        
        fileName=paste0("FKWCResults/KSPAT_D/kspat_dist_",dist,"_deriv_depths_k_exp_alpha_p05_",alphai,"_beta_p5_p5_N_",N1,"_num_runs_50.Rda",sep="")   
        att4=try(load(fileName),T)
        
        fileName=paste0("FKWCResults/KSPAT_ND/kspat_dist_",dist,"_noderiv_depths_k_exp_alpha_p05_",alphai,"_beta_p5_p5_N_",N1,"_num_runs_50.Rda",sep="")   
        att5=try(load(fileName),T)
        att2=att3=att4=att5="d"
      }
      else{
        att2=att3=att4=att5="d"
      }
      fileName=paste0("OTHER_TEST_RESULTS/dist_",dist,"_other_k_exp_alpha1_p05_",alphai,"_p05_beta_p5_p5_p5_N_",N1,"_num_runs_200.Rda",sep="")
           att6=try(load(fileName),T)
      
      if(any(unlist(lapply(list(att1,att2,att3,att4,att5,att6,att7),is.error)))){
        print("**********************************")
        print(paste0("alpha ",alpha,"dist ",dist,"N ",N1))
        print(which(unlist(lapply(list(att1,att2,att3,att4,att5,att6,att7),is.error))))
        print("***********************************")
        print("***********************************")
      }
    }
  }
}





