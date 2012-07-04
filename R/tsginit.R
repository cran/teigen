tsginit <-
function(p,G,x,mug,zmat,n,ng){
	sg <- array(0, dim=c(p, p, G))
	for(g in 1:G){
		sg[,,g] <- cov.wt(x,wt=zmat[,g],center=mug[g,],method="ML")$cov
	}
	sg 
}

