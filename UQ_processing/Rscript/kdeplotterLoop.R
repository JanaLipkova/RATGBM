source("kdepairs.default.R")
source("kdepairs.R")
par(cex.axis=2)


for (GenID in 0:22){

test<-read.table( sprintf("../SyntheticSmall/Samples/ConveretedSamples/SyntSmall_gen_%s.txt", GenID) )

png(filename =sprintf("../SyntheticSmall/Samples/LoopOut/SyntSmall_gen_%s.png", GenID),width=1200,height=1200,units="px",pointsize=12,res=150,type=c("quartz"))


kdepairs(test,contour=FALSE,n=500,labels=c(expression(D[w]),expression(rho),expression(T), expression(sigma) ))

#kdepairs(test,contour=FALSE,n=100,labels=c(expression(D[w]),expression(rho),expression(T), expression(ic[x]),expression(ic[y]),expression(ic[z]),expression(sigma), expression(b), expression(uc[T1gD] ),expression(uc[FLAIR]), expression(sigma[alpha])))
}

dev.off()
