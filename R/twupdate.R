twupdate <-
function(x,n,G,mug,sigmainv,dhfgs78,p,univar,sigma,delta){
	w <- matrix(0,n,G)
	if(!univar){
		for(g in 1:G){
			w[,g] <- (dhfgs78[g]+p)/(dhfgs78[g]+delta[,g])
		}
	}
	else{
		for(g in 1:G){
			w[,g] <- (dhfgs78[g]+p)/(dhfgs78[g]+delta[,g])
		}
	}
	w
}

