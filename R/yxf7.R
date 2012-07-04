yxf7 <-
function(mod,dhfgs78,ng,zmat,w,G,p,n,x,mug,sigmainv){
	if(substring(mod,4,4)=="U"|mod=="univUU"|mod=="univCU"){
		dfoldg <- dhfgs78
		for(g in 1:G){
			dhfgs78[g] <- uniroot( function(v) 1 - digamma(v/2) + log(v/2) + 
					(1/ng[g]) * sum(zmat[,g]*(log(w[,g])-w[,g])) +
					digamma((dfoldg[g]+p)/2) - log((dfoldg[g]+p)/2), 
					lower=0.0000000001, upper=10000)$root
			if(dhfgs78[g]>200){
				dhfgs78[g]<-200
			}
			if(dhfgs78[g]<2){
				dhfgs78[g]<-2
			}
		}
	}
	if(substring(mod,4,4)=="C"|mod=="univUC"|mod=="univCC"){
		dfoldg <- dhfgs78[1]
		dfsamenewg <- uniroot( function(v) 1 - digamma(v/2) + log(v/2) + 
				(1/n) * sum(zmat * (log(w)-w)) +
				digamma((dfoldg+p)/2) - log((dfoldg+p)/2), 
				lower=0.00000000001, upper=10000)$root
		if(dfsamenewg>200){
			dfsamenewg<-200
		}
		if(dfsamenewg<2){
			dfsamenewg<-2
		}
		dhfgs78 <- c(rep(dfsamenewg,G))
	}
	dhfgs78
}

