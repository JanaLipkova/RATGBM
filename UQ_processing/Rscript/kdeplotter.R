source("kdepairs.default.R")
source("kdepairs.R")
par(cex.axis=2)

test<-read.table("../Patient20/Samples/Patient_20_6K_final.txt")

png(filename ="../Patient20/Samples/Patient_20_6K_points.png",width=1200,height=1200,units="px",pointsize=12,res=150,type=c("quartz"))


kdepairs(test,contour=FALSE,n=300,labels=c(expression(D[w]),expression(rho),expression(T), expression(sigma) ))

#kdepairs(test,contour=FALSE,n=100,labels=c(expression(D[w]),expression(rho),expression(T), expression(ic[x]),expression(ic[y]),expression(ic[z]),expression(sigma), expression(b), expression(uc[T1gD] ),expression(uc[FLAIR]), expression(sigma[alpha])))


dev.off()
