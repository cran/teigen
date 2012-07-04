twinit <-
function(x,n,G,mug,sigmainv,vg,p,sg,zmat,univar,sigma){
	w <- matrix(0,n,G)
	delta <- matrix(0,n,G)
	if(!univar){
		for(g in 1:G){
			delta[,g] <- mahalanobis(x, mug[g,], sigmainv[,,g], inverted=TRUE)
			w[,g] <- (vg[g]+p)/(vg[g]+delta[,g])
		}
	}
	else{
		for(g in 1:G){
			delta[,g] <- mahalanobis(x, mug[g,], 1/sigma[,,g],inverted=TRUE)
			w[,g] <- (vg[g]+p)/(vg[g]+delta[,g])
		}
	}
	w
}

