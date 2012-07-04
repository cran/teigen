tmugupdate <-
function(G,zmat,w,x,p,univar){
	mug <- matrix(0,G,p)
	if(univar){
		for(g in 1:G){
			mug[g,1] <- colSums(zmat[,g]*w[,g]*x)/sum(zmat[,g]*w[,g])
		}
	}
	else{
		for(g in 1:G){
			mug[g,] <- colSums(zmat[,g]*w[,g]*x)/sum(zmat[,g]*w[,g])
		}
	}
	mug
}

