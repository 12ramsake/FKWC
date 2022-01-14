# k_xy<-function(x,y,sig){
#   
#   exp(-norm.fd(x-y)^2/(sig^2))
#   
# }

##takes fda obj
K_spatial_depth<-function(fdat,N){
  
  #each col is a fn
  # mat=eval.fd(grid,fdat)
  
  #each col is a fn
  # mat=eval.fd(grid,fdat)
  pairs=combn(N,2)
  
  tmp=function(inds){fdat[inds[1]]-fdat[inds[2]]}
  diffs=apply(pairs,2,tmp)
  
  norms=unlist(lapply(diffs,norm.fd))
  sig=quantile(norms,.15)
  
  dist_mat=matrix(1,nrow=N,ncol=N)
  for(i in 1:ncol(pairs)){
    dist_mat[pairs[1,i],pairs[2,i]]=exp(-norms[i]^2/(sig^2))
    dist_mat[pairs[2,i],pairs[1,i]]=dist_mat[pairs[1,i],pairs[2,i]]
  }
  inds=1:N
  # fn_k=Vectorize(function(i,j){k_xy(fdat[i],fdat[j],sig)},vectorize.args = c('i','j'))
  # dist_mat=outer(inds,inds,fn_k)
  

  #kernel entry
  get_sum_entry<-function(x,y,z){
    num=dist_mat[x,x]+dist_mat[y,z]-dist_mat[x,y]-dist_mat[x,z]
    denom1=dist_mat[x,x]+dist_mat[y,y]-2*dist_mat[x,y]
    denom2=dist_mat[x,x]+dist_mat[z,z]-2*dist_mat[x,z]
    return(num/sqrt(denom1*denom2))
  }
  
  
  get_depth<-function(ind){
   
    pairs_ind=which(mapply("&&",pairs[1,]!=ind,pairs[2,]!=ind))
    pairs_red=pairs[,pairs_ind]
    vals=rbind(rep(ind,dim(pairs_red)[2]),pairs_red)
    summ=sum(mapply(get_sum_entry,rep(ind,dim(pairs_red)[2]),pairs_red[1,],pairs_red[2,]))^.5
    return(1-(1/N)*summ)
    
  }
  
  depths=sapply(1:N,get_depth)
  
  return(depths) 
}




##takes fda obj
spatial_depth<-function(fdat,N){
  
  #each col is a fn
  # mat=eval.fd(grid,fdat)
  pairs=combn(N,2)
  
  tmp=function(inds){fdat[inds[1]]-fdat[inds[2]]}
  diffs=apply(pairs,2,tmp)
  
  norms=unlist(lapply(diffs,norm.fd))
  
  diffs2=lapply(1:length(diffs),function(x){diffs[[x]]*(1/norms[x])})
  
  
  get_dp=function(ind){
    a=which(pairs[1,]==ind)
    b=which(pairs[2,]==ind)
    # inds=c(a,b)
    
    if( length(a)>0)
      mean_d_1=Reduce("+",diffs2[a])*(1/N)
    else 
      mean_d_1=0
    if( length(b)>0)
      mean_d_2=-Reduce("+",diffs2[b])*(1/N)
    else 
      mean_d_2=0
    
    dp=1-norm.fd(mean_d_1+mean_d_2)
    
    return(dp)
  }
  
  depths=sapply(1:N,get_dp)
  return(depths) 
}


