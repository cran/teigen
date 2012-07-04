tkmeansz <-
function(x,n,G,known,kno,testindex,clas){
	kclus <- kmeans(x,G)$cluster
	zmat <- matrix(0,n,G)
	if(clas>0){
		matchtab <- table(known[testindex],kclus[testindex])
		matchit <- NA
		for(i in 1:G){
			matchit[i] <- which.max(matchtab[,i])

		}
		if(length(unique(matchit))==length(matchit)){
			knew <- kclus
			for(i in 1:G){
				knew[kclus==i] <- matchit[i]
			}
			kclus <- knew
		}
	}
	for(i in 1:G){
		zmat[kclus==i, i]<-1
	}
	zmat
}

