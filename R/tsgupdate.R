tsgupdate <-
function(p,G,n,x,mug,zmat,w,ng,mod,pig,submod13){
	sg <- array(0, dim=c(p, p, G))
	sgc <- matrix(0,p,p)
	for(g in 1:G){
		sg[,,g] <- cov.wt(x,wt=zmat[,g]*w[,g],center=mug[g,],method="ML")$cov*(sum(zmat[,g]*w[,g])/ng[g])
	}
	if(any(submod13==c("CCC","CIC","CII")) | mod=="univCC"|mod=="univCU"){
		for(g in 1:G){
			sgc <- sgc + pig[g]*sg[,,g]
		}
			sg[,,] <- sgc
	}
	sg	
}

