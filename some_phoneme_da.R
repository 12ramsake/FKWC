library(fds)
library(fda)
library(fda.usc)
data(aa)
data(ao)
data(dcl)
data(iy)
data(sh)
# 
plot(aa)
plot(ao)
plot(dcl)
plot(iy)
plot(sh)
# matplot(aa$y,type='l',col=1:10)
# 
# basisobj <- create.bspline.basis(c(0,1),50)
# #  convert to a functional data object that interpolates the data.
# result <- smooth.basis(seq(0,1,l=150), aa$y, basisobj)
# yfd  <- result$fd
# 
# #  set up a functional parameter object with smoothing
# #  parameter 1e-6 and a penalty on the 3rd derivative.
# yfdPar <- fdPar(yfd, 2, 1e-6)
# yfd1 <- smooth.fd(yfd, yfdPar)
# plot(yfd1)
# 
# 
# # basisobj1 <- create.bspline.basis(c(0,1),50)
# #  convert to a functional data object that interpolates the data.
# result1 <- smooth.basis(seq(0,1,l=150), sh$y, basisobj)
# yfd11  <- result1$fd
# 
# #  set up a functional parameter object with smoothing
# #  parameter 1e-6 and a penalty on the 3rd derivative.
# yfdPar1 <- fdPar(yfd11, 2, 1e-5)
# yfd111<- smooth.fd(yfd11, yfdPar1)
# plot(yfd111)
# # reg=register.fd(yfd111)
# # plot(reg$regfd)
# # plot(reg$warpfd)
# # # plot(reg$Wfd)
# # fda::AmpPhaseDecomp(yfd111)
# 
# # fdaata_obj=fdata(,seq(0,1,l=200))
# # plot(  fdaata_obj)
# MBD_reg=depthTools::MBD(t(eval.fd(seq(0,1,l=100),yfd111)))
# which.max(MBD_reg$MBD)
# # med_reg=reg$regfd[92]+reg$regfd[65]
# # med_reg$coefs=med_reg$coefs/2
# # plot(reg$regfd[65])
# # 
# # plot(reg$regfd,col=c(rep(1,24),2,rep(1,2222)))
# # plot(reg$regfd[33],add=T,col="red",lwd=3)
# 
# # regged=reg$regfd
# 
# yfd111$coefs=sweep(yfd111$coefs,1,yfd111$coef[,103],"-")
# plot(yfd111)
# 
# 
# plot(var.fd(yfd111))
# 
# # plot(  fdaata_obj)
# MBD_reg2=depthTools::MBD(t(eval.fd(seq(0,1,l=100),yfd1)))
# which.max(MBD_reg2$MBD)
# 
# 
# yfd1$coefs=sweep(yfd1$coefs,1,yfd1$coef[,391],"-")
# plot(yfd1)
# plot(var.fd(yfd1))



pre_process=function(ys){
  
  
  basisobj <- create.bspline.basis(c(0,1),50)
  #  convert to a functional data object that interpolates the data.
  results <- smooth.basis(seq(0,1,l=150), ys, basisobj)
  yfds  <- results$fd
  
  #  set up a functional parameter object with smoothing
  #  parameter 1e-6 and a penalty on the 3rd derivative.
  yfdPars <- fdPar(yfds, 2, 1e-5)
  yfd1s <- smooth.fd(yfds, yfdPars)
  reg=register.fd(yfd1s)
  # MBD_reg=depthTools::MBD(t(eval.fd(seq(0,1,l=100),reg$regfd)))
  d_RPD=depth.RPD(fdata(t(eval.fd(seq(0,1,l=100),reg$regfd))))$dep
  medz=which.max( d_RPD)
  mdc=reg$regfd$coef[,medz]
  reg$regfd$coefs=sweep(reg$regfd$coefs,1,reg$regfd$coef[,medz],"-")
  return(list(curves=reg,med_coef= mdc))
}




aa_pp=pre_process(aa$y)
ao_pp=pre_process(ao$y)
dcl_pp=pre_process(dcl$y)
iy_pp=pre_process(iy$y)
sh_pp=pre_process(sh$y)

# 
# save(aa_pp,
#      ao_pp,
#      dcl_pp,
#      iy_pp,
#      sh_pp,file="phoneme_pp_rpd.Rda")
# 
# load("phoneme_pp.Rda")
load("phoneme_pp_rpd.Rda")
plot(aa_pp$curves$regfd,ylim=c(-6,6))
plot(ao_pp$curves$regfd,ylim=c(-6,6))
plot(dcl_pp$curves$regfd,ylim=c(-6,6))
plot(iy_pp$curves$regfd,ylim=c(-6,6))
plot(sh_pp$curves$regfd,ylim=c(-6,6))

library(RColorBrewer)
pink_grad=colorRampPalette(c("#470f32","#821e5c", "#d63696","#ffd1ee"))

colfunc <- colorRampPalette(c("yellow","red","magenta", "purple","navyblue"))
colfunc2 <- colorRampPalette(c("#c97c00","#990b0b","#821e5c","#d63696", "#5a18ab","#0e1575"))


av=seq(0,1,l=1000)

aa_ev=eval.fd(av, aa_pp$curves$regfd[1:200])
ao_ev=eval.fd(av, ao_pp$curves$regfd[1:200])
dcl_ev=eval.fd(av, dcl_pp$curves$regfd[1:200])
iy_ev=eval.fd(av, iy_pp$curves$regfd[1:200])
sh_ev=eval.fd(av, sh_pp$curves$regfd[1:200])

make_plot=function(vals,cols,main="aa"){
  matplot(av,vals,type='l', ylim=c(-8,8),col=cols,
          ylab="",xlab="",lty=1,bty="n",
          cex.lab=2,yaxt="n",lwd=2,main=main, cex.main=3)
  mtext(expression(italic("t")),1,1.5,cex=2)
  mtext("Log-Periodogram",2,1.5,cex=2)
  axis(2,line=-1,cex=1.25)
}

make_plot(aa_ev,pink_grad(200))
make_plot(ao_ev,brewer.pal(200, "Greens"))
make_plot(dcl_ev,brewer.pal(200, "Blues"))
make_plot(iy_ev,brewer.pal(200, "Purples"))
make_plot(sh_ev,brewer.pal(200, "Reds"))



Cairo::CairoPDF("aa.pdf",width=12,height=8)
make_plot(aa_ev,pink_grad(200),main='aa')
dev.off()
Cairo::CairoPDF("ao.pdf",width=12,height=8)
make_plot(ao_ev,brewer.pal(200, "Greens"),main='ao')
dev.off()
Cairo::CairoPDF("dcl.pdf",width=12,height=8)
make_plot(dcl_ev,rev(brewer.pal(200, "Blues")),main='dcl')
dev.off()
Cairo::CairoPDF("iy.pdf",width=12,height=8)
make_plot(iy_ev,brewer.pal(200, "Purples"),main='iy')
dev.off()
Cairo::CairoPDF("sh.pdf",width=12,height=8)
make_plot(sh_ev,brewer.pal(200, "Reds"),main='sh')
dev.off()



Cairo::CairoPDF("aa_St.pdf",width=16,height=8)

matplot(av,aa_ev,type='l', ylim=c(-8,8),col=rev(brewer.pal(200, "Blues")),
        ylab="",xlab="",lty=1,bty="n",
        cex.lab=2,yaxt="n",lwd=2,xaxt='n')

dev.off()

# legend(0,6,c("aa"),lty=c(1,2),col=c("#821e5c"),lwd=2,cex=2)





par(mfrow=c(3,2))
# deriv.fd(aa_pp$curves$regfd)
plot(deriv.fd(aa_pp$curves$regfd),ylim=c(-50,50))
plot(deriv.fd(ao_pp$curves$regfd),ylim=c(-50,50))
plot(deriv.fd(dcl_pp$curves$regfd),ylim=c(-50,50))
plot(deriv.fd(iy_pp$curves$regfd),ylim=c(-50,50))
plot(deriv.fd(sh_pp$curves$regfd),ylim=c(-50,50))

library(RColorBrewer)
pink_grad=colorRampPalette(c("#470f32","#821e5c", "#d63696","#ffd1ee"))

colfunc <- colorRampPalette(c("yellow","red","magenta", "purple","navyblue"))
colfunc2 <- colorRampPalette(c("#c97c00","#990b0b","#821e5c","#d63696", "#5a18ab","#0e1575"))

par(mfrow=c(1,1))
pdf("aa_vs_sh.pdf",width=12,height=8)
plot.fd(aa_pp$curves$regfd[1:200],
     ylim=c(-6,6),col=pink_grad(200),
     ylab="",xlab=""
     ,lty=1,bty="n",cex.lab=2,yaxt="n",lwd=2)
plot(sh_pp$curves$regfd[1:200],add=T,
     ylim=c(-6,6),col=brewer.pal(200, "Blues"),ylab="",xlab="",
     lty=2,bty="n",cex.lab=2,yaxt="n",lwd=2,main="Smoothed Data",cex.main=2)
mtext(expression(italic("t")),1,1.5,cex=2)
mtext("Log-Periodogram",2,1.5,cex=2)
axis(2,line=-1,cex=1.25)
legend(0,6,c("aa","sh"),lty=c(1,2),col=c("#821e5c",brewer.pal(9, "Blues")[9]),lwd=2,cex=2)
dev.off()


#plotting the curves for paper
un_centred_aa=aa_pp$curves$regfd
un_centred_aa$coefs=sweep(un_centred_aa$coefs ,1,aa_pp$med_coef,"+")
un_centred_sh=sh_pp$curves$regfd
un_centred_sh$coefs=sweep(un_centred_sh$coefs ,1,sh_pp$med_coef,"+")
un_centred_dcl=dcl_pp$curves$regfd
un_centred_dcl$coefs=sweep(un_centred_dcl$coefs ,1,dcl_pp$med_coef,"+")
plot(un_centred_aa[1:200],col=colfunc2(200),ylab="",xlab="",lty=1,bty="n",cex.lab=2,yaxt="n",lwd=2)
plot(un_centred_sh[1:200],col=colfunc2(200),ylab="",xlab="",lty=1,bty="n",cex.lab=2,yaxt="n",lwd=2)

un_centred_ao=ao_pp$curves$regfd
un_centred_ao$coefs=sweep(un_centred_ao$coefs ,1,ao_pp$med_coef,"+")
un_centred_iy=iy_pp$curves$regfd
un_centred_iy$coefs=sweep(un_centred_iy$coefs ,1,iy_pp$med_coef,"+")




un_centred_aa_ev=eval.fd(av, un_centred_aa[1:200])
un_centred_ao_ev=eval.fd(av, un_centred_ao[1:200])
un_centred_dcl_ev=eval.fd(av, un_centred_dcl[1:200])
un_centred_iy_ev=eval.fd(av, un_centred_iy[1:200])
un_centred_sh_ev=eval.fd(av, un_centred_sh[1:200])

make_plot=function(vals,cols,main="aa"){
  matplot(av,vals,type='l', ylim=c(-0,25),col=cols,
          ylab="",xlab="",lty=1,bty="n",
          cex.lab=2,yaxt="n",lwd=2,main=main, cex.main=3)
  mtext(expression(italic("t")),1,1.5,cex=2)
  mtext("Log-Periodogram",2,1.5,cex=2)
  axis(2,line=-1,cex=1.25)
}

make_plot(un_centred_aa_ev,pink_grad(200))
make_plot(un_centred_ao_ev,brewer.pal(200, "Greens"))
make_plot(un_centred_dcl_ev,brewer.pal(200, "Blues"))
make_plot(un_centred_iy_ev,brewer.pal(200, "Purples"))
make_plot(un_centred_sh_ev,brewer.pal(200, "Reds"))



Cairo::CairoPDF("un_centred_aa.pdf",width=12,height=8)
make_plot(un_centred_aa_ev,pink_grad(200),main='aa')
dev.off()
Cairo::CairoPDF("un_centred_ao.pdf",width=12,height=8)
make_plot(un_centred_ao_ev,brewer.pal(200, "Greens"),main='ao')
dev.off()
Cairo::CairoPDF("un_centred_dcl.pdf",width=12,height=8)
make_plot(un_centred_dcl_ev,rev(brewer.pal(200, "Blues")),main='dcl')
dev.off()
Cairo::CairoPDF("un_centred_iy.pdf",width=12,height=8)
make_plot(un_centred_iy_ev,brewer.pal(200, "Purples"),main='iy')
dev.off()
Cairo::CairoPDF("un_centred_sh.pdf",width=12,height=8)
make_plot(un_centred_sh_ev,brewer.pal(200, "Reds"),main='sh')
dev.off()







make_plot(aa_ev,pink_grad(200))
make_plot(ao_ev,brewer.pal(200, "Greens"))
make_plot(dcl_ev,brewer.pal(200, "Blues"))
make_plot(iy_ev,brewer.pal(200, "Purples"))
make_plot(sh_ev,brewer.pal(200, "Reds"))












pdf("aa_vs_dcl_uc.pdf",width=12,height=8)
plot(un_centred_dcl[1:400],col=pink_grad(400),ylab="",xlab="",lty=1,bty="n",cex.lab=2,yaxt="n",lwd=2,ylim=c(-5,25),main="Smoothed Data",cex.main=2)
plot(un_centred_aa[1:400],col=rev(brewer.pal(9, "Blues")),ylab="",xlab="",lty=1,bty="n",cex.lab=2,yaxt="n",lwd=2,add=T)
mtext(expression(italic("t")),1,1.5,cex=2)
mtext("Log-Periodogram",2,1.5,cex=2)
axis(2,line=-1,cex=1.25)
legend_image <- grDevices::as.raster(t(matrix(pink_grad(400), ncol=400,nrow=400)))
graphics::rasterImage(legend_image,.4,-5,.6,-4.5)
text(.62,-4.8,"dcl",cex=1.7)
legend_image2 <- grDevices::as.raster(t(matrix(rev(brewer.pal(9, "Blues")), ncol=9,nrow=9)))
graphics::rasterImage(legend_image2,.4,-5+1,.6,-4.5+1)
text(.62,-4.8+1,"aa",cex=1.7)
# legend(0,6,c("aa","sh"),lty=c(1,2),col=c("#821e5c",brewer.pal(9, "Blues")[9]),lwd=2,cex=2)
dev.off()



pdf("aa_vs_dcl.pdf",width=12,height=8)
plot(dcl_pp$curves$regfd[1:400],col=pink_grad(400),ylab="",xlab="",lty=1,bty="n",cex.lab=2,yaxt="n",lwd=2,ylim=c(-10,10),main="Centered and Smoothed Data",cex.main=2)
plot(aa_pp$curves$regfd[1:400],col=rev(brewer.pal(9, "Blues")),ylab="",xlab="",lty=2,bty="n",cex.lab=2,yaxt="n",lwd=2,add=T)
mtext(expression(italic("t")),1,1.5,cex=2)
mtext("Log-Periodogram",2,1.5,cex=2)
axis(2,line=-1,cex=1.25)
legend_image <- grDevices::as.raster(t(matrix(pink_grad(400), ncol=400,nrow=400)))
graphics::rasterImage(legend_image,.4,-10,.6,-9.5)
text(.62,-9.8,"dcl",cex=1.7)
legend_image2 <- grDevices::as.raster(t(matrix(rev(brewer.pal(9, "Blues")), ncol=9,nrow=9)))
graphics::rasterImage(legend_image2,.4,-10+1,.6,-9.5+1)
text(.62,-9.8+1,"aa",cex=1.7)
dev.off()






ker_aa=var.fd(aa_pp$curves$regfd)
plot(ker_aa)
ker_ao=var.fd(ao_pp$curves$regfd)
plot(ker_ao)
ker_dcl=var.fd(dcl_pp$curves$regfd)
plot(ker_dcl)
ker_iy=var.fd(iy_pp$curves$regfd)
plot(ker_iy)
ker_sh=var.fd(sh_pp$curves$regfd)
plot(ker_sh)

ker_aad=var.fd(deriv.fd(aa_pp$curves$regfd))
plot(ker_aad)
ker_aod=var.fd(deriv.fd(ao_pp$curves$regfd))
plot(ker_aod)
ker_dcld=var.fd(deriv.fd(dcl_pp$curves$regfd))
plot(ker_dcld)



ddcl=eval.fd(seq(0,1,l=100),deriv.fd(dcl_pp$curves$regfd))
ddcl[-c(1,2,10,11)]
matplot((ddcl[-c(1:3,97:100),]),type='l')
ddcl2=fdata(t(ddcl[-c(1:3,97:100),]),seq(0,1,l=100)[-c(1:3,97:100)])
plot(ddcl2)



par(mfrow=c(3,2))
plot(aa_pp$curves$warpfd)
plot(ao_pp$curves$warpfd)
plot(dcl_pp$curves$warpfd)
plot(iy_pp$curves$warpfd)
plot(sh_pp$curves$warpfd)




# fdata_obj1=fdata(t(eval.fd(seq(0,1,l=100),dcl_pp$curves$regfd)))
depts=depth.RPD(ddcl2)$dep                  
depts_ordering=rank(depts)


plot(fdata2fd(ddcl2)[depts_ordering>5])
ker_dcld=var.fd(fdata2fd(ddcl2)[depts_ordering>5])
plot(ker_dcld)


ker_iyd=var.fd(deriv.fd(iy_pp$curves$regfd))
plot(ker_iyd)
ker_shd=var.fd(deriv.fd(sh_pp$curves$regfd))
plot(ker_shd)

# comb=depthTools::MBD(rbind(t(eval.fd(seq(0,1,l=100),aa_pp$curves$regfd)),
#                            t(eval.fd(seq(0,1,l=100),ao_pp$curves$regfd)),
#                            t(eval.fd(seq(0,1,l=100),dcl_pp$curves$regfd)),
#                            t(eval.fd(seq(0,1,l=100),iy_pp$curves$regfd)),
#                            t(eval.fd(seq(0,1,l=100),sh_pp$curves$regfd))))
# kruskal.test(comb$ordering,g=rep(c(1,2,3,4,5),each=400))
# rank_frame=data.frame(cbind(comb$ordering,rep(c(1,2,3,4,5),each=400)))
# aggregate(rank_frame[,1],by=list(rank_frame$X2),mean)







fdata_obj=fdata(rbind(t(eval.fd(seq(0,1,l=100),aa_pp$curves$regfd)),
                      t(eval.fd(seq(0,1,l=100),ao_pp$curves$regfd)),
                      t(eval.fd(seq(0,1,l=100),dcl_pp$curves$regfd)),
                      t(eval.fd(seq(0,1,l=100),iy_pp$curves$regfd)),
                      t(eval.fd(seq(0,1,l=100),sh_pp$curves$regfd))),seq(0,1,l=100))



groups=list(fdata_obj[1:400],fdata_obj[401:800],fdata_obj[801:1200],fdata_obj[1201:1600],fdata_obj[1601:2000])

steel_test=function(groups,RPD=T){
  
  pvalues=matrix(1,length(groups),length(groups))
  
  for(i in 2:length(groups)){
    for(j in 1:(i-1)){
      ni=length(groups[[i]])
      nj=length(groups[[j]])
      if(RPD)
        depths=depth.RPD(c(groups[[i]],groups[[j]]))$dep
      else
        depths=norm.fdata(c(groups[[i]],groups[[j]]))
      pval=wilcox.test(depths[1:ni],depths[(ni+1):(ni+nj)],correct=F)$p.value
      pvalues[i,j]=pvalues[j,i]=   pval
    }
  }
  return(pvalues)
}



RPD_rnks=depth.RPD(fdata_obj)
RPD_rnks2=rank(RPD_rnks$dep)

kruskal.test(RPD_rnks2,g=rep(c(1,2,3,4,5),each=400))
rank_frameRPD=data.frame(cbind(RPD_rnks2,rep(c(1,2,3,4,5),each=400)))
aggregate(rank_frameRPD[,1],by=list(rank_frameRPD$V2),mean)
res_rpd=dunn.test::dunn.test(RPD_rnks2,g=rep(c(1,2,3,4,5),each=400),method="bonferroni")
sidak_pvalues_rpd=1-(1-res_rpd$P)^20
res_rpd=steel_test(groups)
res_rpd=1-(1-res_rpd)^22
# res_rpd[1,2]=res_rpd[2,1]=sidak_pvalues_rpd[1]
# res_rpd[1,3]=res_rpd[3,1]=sidak_pvalues_rpd[2]
# res_rpd[2,3]=res_rpd[3,2]=sidak_pvalues_rpd[3]
# res_rpd[1,4]=res_rpd[4,1]=sidak_pvalues_rpd[4]
# res_rpd[2,4]=res_rpd[4,2]=sidak_pvalues_rpd[5]
# res_rpd[3,4]=res_rpd[4,3]=sidak_pvalues_rpd[6]
# res_rpd[1,5]=res_rpd[5,1]=sidak_pvalues_rpd[7]
# res_rpd[2,5]=res_rpd[5,2]=sidak_pvalues_rpd[8]
# res_rpd[3,5]=res_rpd[5,3]=sidak_pvalues_rpd[9]
# res_rpd[4,5]=res_rpd[5,4]=sidak_pvalues_rpd[10]
xtable(res_rpd)


LTR_rnks=norm.fdata(fdata_obj)^2
kruskal.test(LTR_rnks,g=rep(c(1,2,3,4,5),each=400))
LTR_rnkss=data.frame(cbind(rank(LTR_rnks),rep(c(1,2,3,4,5),each=400)))
aggregate(LTR_rnkss[,1],by=list(LTR_rnkss$X2),mean)
resltr=dunn.test::dunn.test(LTR_rnks,g=rep(c(1,2,3,4,5),each=400),method="bonferroni")
LTR_rnkss$X2=LTR_rnkss$X2%>%as.factor()
res_ltr=steel_test(groups,F)


library(xtable)
res_ltr=1-(1-res_ltr)^22
# res_ltr=matrix(1,ncol=5,nrow=5)
# res_ltr[1,2]=res_ltr[2,1]=sidak_pvalues_ltr[1]
# res_ltr[1,3]=res_ltr[3,1]=sidak_pvalues_ltr[2]
# res_ltr[2,3]=res_ltr[3,2]=sidak_pvalues_ltr[3]
# res_ltr[1,4]=res_ltr[4,1]=sidak_pvalues_ltr[4]
# res_ltr[2,4]=res_ltr[4,2]=sidak_pvalues_ltr[5]
# res_ltr[3,4]=res_ltr[4,3]=sidak_pvalues_ltr[6]
# res_ltr[1,5]=res_ltr[5,1]=sidak_pvalues_ltr[7]
# res_ltr[2,5]=res_ltr[5,2]=sidak_pvalues_ltr[8]
# res_ltr[3,5]=res_ltr[5,3]=sidak_pvalues_ltr[9]
# res_ltr[4,5]=res_ltr[5,4]=sidak_pvalues_ltr[10]
xtable(res_ltr)
tab=data.frame(cbind(res_rpd,res_ltr))
names(tab)=rep(c('aa ',' ao ',' dcl ',' iy ',' sh ',' aa ',' ao ',' dcl ',' iy ',' sh'),1)
rownames(tab)=c('aa ',' ao ',' dcl ',' iy ',' sh ')
xtable(tab)


plot(center.fd(aa_pp$curves$warpfd))
plot(center.fd(ao_pp$curves$warpfd))
plot(center.fd(dcl_pp$curves$warpfd))
plot(center.fd(iy_pp$curves$warpfd))
plot(center.fd(sh_pp$curves$warpfd))


fdata_warp_obj=fdata(rbind(t(eval.fd(seq(0,1,l=100),center.fd(aa_pp$curves$warpfd))),
                      t(eval.fd(seq(0,1,l=100),center.fd(ao_pp$curves$warpfd))),
                      t(eval.fd(seq(0,1,l=100),center.fd(dcl_pp$curves$warpfd))),
                      t(eval.fd(seq(0,1,l=100),center.fd(iy_pp$curves$warpfd))),
                      t(eval.fd(seq(0,1,l=100),center.fd(sh_pp$curves$warpfd)))),seq(0,1,l=100))
plot(fdata_warp_obj)
depth_ltr=norm.fdata(fdata_warp_obj)^2
dunn.test::dunn.test(depth_ltr,g=rep(c(1,2,3,4,5),each=400),method="bonferroni")


par(mfrow=c(2,3))

#we wanna trim the data, and then divide by square root of k(t,t) , the trace of the covariance

eval_trace=function(kernel,grid=seq(0,1,l=100)){
  tr=mapply(function(x,y){eval.bifd(x,y,kernel)},grid,grid)
  plot(tr,type='l')
  return(tr)
}

ker_aa=var.fd(aa_pp$curves$regfd)
# plot(ker_aa,type='l')
eval_trace(ker_aa)
ker_ao=var.fd(ao_pp$curves$regfd)
# plot(ker_ao)
eval_trace(ker_ao)
ker_dcl=var.fd(dcl_pp$curves$regfd[-101])
# plot(ker_dcl)
eval_trace(ker_dcl)
ker_iy=var.fd(iy_pp$curves$regfd)
# plot(ker_iy)
eval_trace(ker_iy)
ker_sh=var.fd(sh_pp$curves$regfd)
# plot(ker_sh)
eval_trace(ker_sh)


tr_aa=eval_trace(ker_aa)%>%sqrt()
aa_pp_norm=aa_pp
aa_pp_norm$curves$regfd$coefs=sweep(aa_pp$curves$regfd$coefs,2,tr_aa,"/")
plot(aa_pp_norm$curves$regfd,ylim=c(-6,6))


tr_dcl=eval_trace(ker_dcl)%>%sqrt()
dcl_pp_norm=dcl_pp
dcl_pp_norm$curves$regfd$coefs=sweep(dcl_pp$curves$regfd$coefs,2,tr_dcl,"/")
plot(dcl_pp_norm$curves$regfd[-101],ylim=c(-6,6))
dcl_pp_norm$curves$regfd[-101]


tr_iy=eval_trace(ker_iy)%>%sqrt()
iy_pp_norm=iy_pp
iy_pp_norm$curves$regfd$coefs=sweep(iy_pp$curves$regfd$coefs,2,tr_iy,"/")
plot(iy_pp_norm$curves$regfd[-101],ylim=c(-6,6))
tr_ao=eval_trace(ker_ao)%>%sqrt()
ao_pp_norm=ao_pp
ao_pp_norm$curves$regfd$coefs=sweep(ao_pp$curves$regfd$coefs,2,tr_ao,"/")
plot(ao_pp_norm$curves$regfd[-101],ylim=c(-6,6))
tr_sh=eval_trace(ker_sh)%>%sqrt()

sh_pp_norm=sh_pp
sh_pp_norm$curves$regfd$coefs=sweep(sh_pp$curves$regfd$coefs,2,tr_sh,"/")
plot(sh_pp_norm$curves$regfd,ylim=c(-6,6))





fdata_obj_norm=fdata(rbind(t(eval.fd(seq(0,1,l=100),aa_pp_norm$curves$regfd)),
                      t(eval.fd(seq(0,1,l=100),ao_pp_norm$curves$regfd)),
                      t(eval.fd(seq(0,1,l=100),dcl_pp_norm$curves$regfd)),
                      t(eval.fd(seq(0,1,l=100),iy_pp_norm$curves$regfd)),
                      t(eval.fd(seq(0,1,l=100),sh_pp_norm$curves$regfd))))



plot(fdata_obj_norm)

deps_norm=depth.RPD(fdata_obj_norm)$dep
plot(deps_norm,col=rep(c(1,2,3,4,5),each=400))
kruskal.test(deps_norm,g=rep(c(1,2,3,4,5),each=400))
dunn.test::dunn.test(deps_norm,g=rep(c(1,2,3,4,5),each=400),method="bonferroni")


deps_norm_2=norm.fdata(fdata_obj_norm)^2
plot(deps_norm_2,col=rep(c(1,2,3,4,5),each=400),ylim=c(0,100))
colMeans(matrix(rank(deps_norm_2),ncol=5))

kruskal.test(rank(deps_norm_2),g=rep(c(1,2,3,4,5),each=400))
plot(deps_norm_2,col=rep(c(1,2,3,4,5),each=400))



plot(deriv.fd(dcl_pp_norm$curves$regfd[-101]),ylim=c(-50,50))
plot(deriv.fd(aa_pp_norm$curves$regfd),ylim=c(-50,50))


datamx=rbind(t(eval.fd(seq(0,1,l=100),aa_pp$curves$regfd)),
      t(eval.fd(seq(0,1,l=100),ao_pp$curves$regfd)),
      t(eval.fd(seq(0,1,l=100),dcl_pp$curves$regfd)),
      t(eval.fd(seq(0,1,l=100),iy_pp$curves$regfd)),
      t(eval.fd(seq(0,1,l=100),sh_pp$curves$regfd)))

sample_n=rep(400,5)
PV_rps <- rcpparma_ECFRPall(t(datamx),sample_n,1000,1,1)













