tuniformz <-
function(n,G,clas,kno,known){
	zmat <- matrix(0, n, G)
	if(clas>0){
		for(i in 1:G){
				zmat[as.numeric(known)==i,i] <- 1
		}
		zmat[kno==0,] <- 1/G
	}
	zmat
}

