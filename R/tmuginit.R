tmuginit <-
function(G,p,x,zmat,ng,univar){
	mug <- matrix(0,G,p)
		for(g in 1:G){
			mug[g,] <- colSums(zmat[,g]*x)/ng[g]
		}
	mug
}

