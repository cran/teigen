tgivenz <-
function(n,G,known,initg,testindex,clas){
	zmat <- matrix(0,n,G)
	if(clas>0){
		matchtab <- table(known[testindex],initg[testindex])
		matchit <- NA
		for(i in 1:G){
			matchit[i] <- which.max(matchtab[,i])
			
		}
		if(length(unique(matchit))==length(matchit)){
			initnew <- initg
			for(i in 1:G){
				initnew[initg==i] <- matchit[i]
			}
			initg <- initnew
		}
	}
	for(i in 1:G){
		zmat[initg==i, i]<-1
	}
	zmat
}

