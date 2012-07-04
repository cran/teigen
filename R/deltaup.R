deltaup <-
function(x,mug,sigma,sigmainv,G,n,univar){
	delta <- matrix(0,n,G)
	if(univar){
		for(g in 1:G){
			delta[,g] <- mahalanobis(x, mug[g,],1/sigma[,,g], inverted=TRUE)
		}
	}
	else{
		for(g in 1:G){
			delta[,g] <- mahalanobis(x, mug[g,],sigmainv[,,g], inverted=TRUE)
		}
	}
	delta
}

