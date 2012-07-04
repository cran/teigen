tdiscrandz <-
function(n,G){
	zmat <- matrix(0, n, G)
	dum <- sample(1:G,n,replace=TRUE)
	for(i in 1:G){
		zmat[dum==i,i] <- 1
	}
	zmat
}

