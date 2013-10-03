plot.teigen <- function(x, xmarg=1, ymarg=2, res=200, levels=c(seq(.01,1,by=0.01), 0.001), what=c("contour","uncertainty"), ...){
	teigen <- x
	classcolours <- rainbow(length(unique(teigen$class)))
	numscr <- length(what)
	if(ncol(teigen$x)>1){
		if(numscr==2){
			par(mfrow=c(2,1))
		}
		if("contour"%in%what){
			plot(teigen$x[,c(xmarg,ymarg)], col=classcolours[teigen$clas], pch=20, main="Marginal Contour Plot", ...)
      lims <- par()$usr
			xseq <- seq(lims[1], lims[2], length.out=res)
			yseq <- seq(lims[3], lims[4], length.out=res)
			seqmat <- matrix(NA, res^2, 2)
			seqmat[,1] <- rep(xseq, each=res)
			seqmat[,2] <- rep(yseq, res)
			val <- matrix(NA,res,res)
			if(!teigen$info$univar){
				sigmainv <- array(NA, dim=c(2,2,teigen$G))
				for(g in 1:teigen$G){
					sigmainv[,,g] <- solve(teigen$par$sigma[c(xmarg,ymarg),c(xmarg,ymarg),g])
				}
			}
			delt <- deltaup(seqmat, matrix(teigen$par$mean[,c(xmarg,ymarg)],nrow=teigen$G,ncol=2), teigen$par$sigma[c(xmarg,ymarg),c(xmarg,ymarg),],sigmainv, teigen$G, res^2, teigen$info$univar)
			dens <- rowSums(exp(tft(seqmat,teigen$G,colSums(teigen$fuzzy)/nrow(teigen$x),teigen$par$df,2,matrix(teigen$par$mean[,c(xmarg,ymarg)],nrow=teigen$G,ncol=2),sigmainv,res^2,array(teigen$par$sigma[c(xmarg,ymarg),c(xmarg,ymarg),], dim=c(2,2,teigen$G)),teigen$info$univar,delt,teigen$info$gauss)))
			val <- matrix(dens, res, res, byrow=TRUE)
			contour(x=xseq, y=yseq, z=val, add=TRUE, levels=levels, col=rgb(0.5,0.5,0.5,alpha=0.7))
		}
		if("uncertainty"%in%what){
			plot(teigen$x[,c(xmarg,ymarg)], col=classcolours[teigen$clas],cex=(2*(1-apply(teigen$fuzzy,1,max))), pch=20, main="Uncertainty Plot", ...)
		}
	}
	else{
		G <- ncol(teigen$fuzzy)
		pigs <- colSums(teigen$fuzzy)/nrow(teigen$x)
		dunivt <- function(xdum,df,sig,mean,pig){ 
			exp(log(pig)+lgamma((df+1)/2)-(1/2)*log(sig)-((1/2)*(log(pi)+log(df))+lgamma(df/2)+((df+1)/2)*(log(1+ mahalanobis(matrix(xdum,ncol=1), mean, 1/sig, inverted=TRUE)/df))))
		}
		bigdens <- NA
		for(g in 1:G){
			bigdens[g] <- dunivt(teigen$par$mean[g,],teigen$par$df[g],teigen$par$sigma[,,g],teigen$par$mean[g,],pigs[g])
		}
		plot(density(teigen$x), ylim=c(0,max(bigdens)+0.04*max(bigdens)), main="Univariate Density Plot")
		#univt <- function(xdum){ log(pig[g])+lgamma((teigen$par$df[g]+1)/2)-(1/2)*log(teigen$par$sigma[,,g])-((p/2)*(log(pi)+log(teigen$par$df[g]))+lgamma(teigen$par$df[g]/2)+((teigen$par$df[g]+p)/2)*(log(1+ mahalanobis(matrix(xdum,nrow(teigen$fuzzy),1), teigen$par$mug[g,], 1/teigen$par$sigma[,,g], inverted=TRUE)/teigen$par$df[g])))}
		for(g in 1:G){
			curve(dunivt(x,teigen$par$df[g],teigen$par$sigma[,,g],teigen$par$mean[g,],pigs[g]),add=TRUE, col=classcolours[g])
		}
	}
	par(mfrow=c(1,1))
}
