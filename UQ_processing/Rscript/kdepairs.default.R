kdepairs.default <- function(x, n=25, density=TRUE, contour=TRUE,labels, ...) {
    require(MASS)
    require(sp)
    ly=length(x)-1
    last=length(x)
    y <- data.frame(x[,1:ly])
    assign("loglh",data.frame(x[,last]) , envir = .GlobalEnv)
    
    #loglh<<-data.frame(x[,last])
    
    
    fun.lower <- function(x1, x2, ...) {
        if (is.factor(x1)) {
            x1 <- as.integer(x1)
        }
        if (is.factor(x2)) {
            x1 <- as.integer(x2)
        }
        OK <- length(unique(x1))>2 && length(unique(x2))>2
        if (!density && !contour)
        n <- 0
        
        if (n>0 && OK) {
            if (density || contour)
            
            #d <- kde2d(x1, x2, n=n)
            # d <- kde2d(x1, x2, h=c( 0.25 *(max(x1) - min(x1)), 0.25 *( max(x2) - min(x2) )  ), n = n,... )
            d <- kde2d(x1, x2, h=c( 0.2 *(max(x1) - min(x1)), 0.2 *( max(x2) - min(x2) )  ), n = n,... )

            if (density) {
                #image(d, col=terrain.colors(50), add=TRUE)
                #  jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
                vv<-bpy.colors(n = 100, cutoff.tails = 0.3, alpha = 0.5)
                
                #vv[1] <- "#FFFFFF00"
                image(d, col=vv, add=TRUE)
                # image(d, xlim = c(0.01*min(x1),10*max(x1)), ylim = c(0.01*min(x2), 100*max(x2)), col= vv, add=TRUE)
                
                #image(d,  xlim = c(min(x1), 2*max(x1)), add=TRUE)
            }
            if (contour)
            graphics:::contour(d,add=TRUE)
        } else points(x1, x2)
        # for Helium: plot the correct answer
        #		points(0.2556,0.141,pch=8,col="green",lwd=2, cex=1.5)
        
            # D-rho
            points(0.7565,0.0375,pch=8,col="green",lwd=2, cex=1.5)
            points(1.5829,0.0673,pch=8,col="red",lwd=2, cex=1.5)
            points(0.7565,0.0673,pch=8,col="black",lwd=2, cex=1.5)
            points(1.5829,0.0375,pch=8,col="yellow",lwd=2, cex=1.5)
            
            # rho-Tend
            points(0.0375,151.3850,pch=8,col="green",lwd=2, cex=1.5)
            points(0.0673,275.2527,pch=8,col="red",lwd=2, cex=1.5)
            points(0.0375,275.2527,pch=8,col="black",lwd=2, cex=1.5)
            points(0.0673,151.3850,pch=8,col="yellow",lwd=2, cex=1.5)
            
            # Tend - sigma
            points(151.3850,0.1945,pch=8,col="green",lwd=2, cex=1.5)
            points(275.2527,0.2365,pch=8,col="red",lwd=2, cex=1.5)
            points(275.2527,0.1945,pch=8,col="black",lwd=2, cex=1.5)
            points(151.3850,0.2365,pch=8,col="yellow",lwd=2, cex=1.5)


    }
    
    
    
    fun.upper <- function(x1, x2, ...) {
        if (is.factor(x1)) {
            x1 <- as.integer(x1)
        }
        if (is.factor(x2)) {
            x1 <- as.integer(x2)
        }
        vv<-bpy.colors(n = 100, cutoff.tails = 0.3, alpha = 0.5)
        #vv<-cm.colors(10)
        #This adds a column of color values
        # based on the y values
        length(loglh)
        loglh
        datc <- vv[as.numeric(cut(as.numeric(loglh[,1]),breaks = 100))]
        points(x1,x2,col=datc)
        #        points(x1,x2, col="lightgrey")
        lines(lowess(x1,x2), col="darkgreen", lty=1)
        COR <- cor(x1, x2)
        text(mean(range(x1,na.rm=TRUE)), mean(range(x2,na.rm=TRUE)),
        round(COR, 3), cex=1+abs(COR))
        
        
        # D-rho
        points(0.0375,0.7565,pch=8,col="green",lwd=2, cex=1.5)
        points(0.0673,1.5829,pch=8,col="red",lwd=2, cex=1.5)
        points(0.0375,1.5829,pch=8,col="black",lwd=2, cex=1.5)
        points(0.0673,0.7565,pch=8,col="yellow",lwd=2, cex=1.5)
        
        # rho-Tend
        points(151.3850,0.0375,pch=8,col="green",lwd=2, cex=1.5)
        points(275.2527,0.0673,pch=8,col="red",lwd=2, cex=1.5)
        points(275.2527,0.0375,pch=8,col="black",lwd=2, cex=1.5)
        points(151.3850,0.0673,pch=8,col="yellow",lwd=2, cex=1.5)
        
        # Tend - sigma
        points(0.1945,151.3850,pch=8,col="green",lwd=2, cex=1.5)
        points(0.2365,275.2527,pch=8,col="red",lwd=2, cex=1.5)
        points(0.1945,275.2527,pch=8,col="black",lwd=2, cex=1.5)
        points(0.2365,151.3850,pch=8,col="yellow",lwd=2, cex=1.5)


    }
    
    
    panel.crap = function(x,...)
    {
        pu <- par("usr")
        
        #d <- density(x,...)
        d <- density(x, bw = 0.06* (max(x) - min(x)), adjust = 1, kernel = "gaussian", weights = NULL, window = kernel, width, give.Rkern = FALSE, ...)

        
        #      par("usr" = c(pu[1:2], 0, max(d$y)*1.9))
        #        lines(d)
        #        polygon(d, col="green", border="black")
        
        par("usr" = c(pu[1:2], 0, 1))
        h<-boxplot(x, at=0.5,varwidth=TRUE, boxwex=0.2, horizontal=TRUE, outline=FALSE,plot=FALSE,col = NULL)
        bxp(h,at=0.5,varwidth=TRUE,boxwex=0.2,horizontal=TRUE,outline=FALSE,add=TRUE,axes=FALSE,col = NULL)
        
        par("usr" = c(pu[1:2], 0, max(d$y)*1.9))
        lines(d)
        polygon(d, col="green", border="black")
        
        
        # .... tmp....
        #points(0.0312,0.141,pch=8,col="black",lwd=2, cex=1.5)
        #points(0.0833,0.141,pch=8,col="black",lwd=2, cex=1.5)
        
        # ... end tmp...
        
        # h<-boxplot(x, at=0.5,varwidth=TRUE, boxwex=0.3, horizontal=TRUE, outline=FALSE,plot=FALSE)
        # bxp(h,at=0.5,varwidth=TRUE,boxwex=0.3,horizontal=TRUE,outline=FALSE,add=TRUE,axes=FALSE)
        
        #		sp = 5.0
        #		mp = 2.0
        #		sn = 1.0
        #		d = 1.0
        #		s <- sqrt(1.0/(1.0/sp^2 + 1.0/sn^2))
        #		m <- (mp/sp^2 + d/sn^2)/(1.0/sp^2 + 1.0/sn^2)
        #		print(m)
        #		print(s^2)
        #		l <- 100*(m-4*s)
        #		r <- 100*(m+4*s)
        #		xa <- array(l:r)/100
        #		da <- (dnorm(xa, mean=m, sd=s))
        #		lines(xa, da, col = "red",type="l")
        
        #		eps <- 0.9
        #		n <- 100
        #		l <- -10*n
        #		r <- 10*n
        #		xa <- array(l:r)/n
        #		epsa <- array(eps, c(1,n))
        #		ya <- pnorm(epsa-xa) - pnorm(-(epsa+xa)) + pnorm(10*(epsa-xa)) - pnorm(-10*(epsa+xa))
    }
    
    panel.hist <- function(x, ...) {
        usr <- par("usr"); on.exit(par(usr))
        par(usr = c(usr[1:2], 0, 1.5) )
        h <- hist(x, breaks = 30,plot=FALSE)
        breaks <- h$breaks; nB <- length(breaks)
        y <- h$density; #y <- y/max(y)
        rect(breaks[-nB], 0, breaks[-1], y, col="green", ...)
        legend("center", "(x,y)", pch=1, title="center")
        box()
        h<-boxplot(x, at=0.5,varwidth=TRUE, boxwex=0.3, horizontal=TRUE, outline=FALSE,plot=FALSE)
        bxp(h,at=0.5,varwidth=TRUE,boxwex=0.3,horizontal=TRUE,outline=FALSE,add=TRUE,axes=FALSE)
        
        #		xa <- sort(x)
        sp = 10.0
        mp = 0.0
        sn = 1.0
        s <- sqrt(1.0/(1.0/sp^2 + 1.0/sn^2))
        m <- (mp/sp^2 + 0.0/sn^2)/(1.0/sp^2 + 1.0/sn^2)
        print(m)
        print(s^2)
        l <- 100*(m-4*s)
        r <- 100*(m+4*s)
        xa <- array(l:r)/100
        da <- (dnorm(xa, mean=m, sd=s))
        lines(xa, da, col = "red",type="l")
    }
    pairs.default(y,labels=labels,lower.panel=fun.lower, upper.panel=fun.upper, diag.panel=panel.crap)
    invisible(NULL)
}
