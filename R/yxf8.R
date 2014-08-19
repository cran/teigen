yxf8 <- function(mod,dhfgs78,ng,zmat,w,G,p,n,x,mug,sigmainv){
	if(substring(mod,4,4)=="U"|mod=="univUU"|mod=="univCU"){
		dfoldg <- dhfgs78
		for(g in 1:G){
      constn <- 1 + (1/ng[g]) * sum(zmat[,g]*(log(w[,g])-w[,g])) + digamma((dfoldg[g]+p)/2) - log((dfoldg[g]+p)/2)
# 			dhfgs78[g] <- uniroot( function(v) log(v/2) - digamma(v/2) + constn , 
# 					lower=0.0001, upper=1000, tol=0.00001)$root
      constn <- -constn
      #print(constn)
      #print(paste("fix", dhfgs78[g] - ((-exp(constn)+2*(exp(constn))*(exp(digamma(dfoldg[g]/2))-(dfoldg[g]/2-1/2)))/(1-exp(constn)))))
      #print(paste("old", dhfgs78[g] - (-exp(constn)/(1-exp(constn)))))
      #print(c((2-sqrt(4+4*4*exp(constn)))/2, (2+sqrt(4+4*4*exp(constn)))/2))
      #print((-exp(constn))/(1-exp(constn)))
      dhfgs78[g] <- (-exp(constn)+2*(exp(constn))*(exp(digamma(dfoldg[g]/2))-(dfoldg[g]/2-1/2)))/(1-exp(constn))
#     dhfgs78[g] <- (-exp(constn))/(1-exp(constn))
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
    constn <- 1 + (1/n) * sum(zmat * (log(w)-w)) + digamma((dfoldg+p)/2) - log((dfoldg+p)/2)
    constn <- -constn
#		dfsamenewg <- uniroot( function(v) log(v/2) - digamma(v/2) + constn, 
#				lower=0.0001, upper=1000, tol=0.01)$root
    dfsamenewg <- (-exp(constn)+2*(exp(constn))*(exp(digamma(dfoldg/2))-(dfoldg/2-1/2)))/(1-exp(constn))
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

