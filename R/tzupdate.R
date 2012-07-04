tzupdate <-
function(x,G,pig,dhfgs78,p,mug,sigmainv,n,sigma,clas,kno,known,unkno,univar,delta){
	num <- matrix(0,n,G)
	if(univar){
		log_num <- matrix(0,n,G)
		for(g in 1:G){
			log_num[,g]<-log(pig[g])+lgamma((dhfgs78[g]+1)/2)-(1/2)*log(sigma[,,g])-
					((p/2)*(log(pi)+log(dhfgs78[g]))+lgamma(dhfgs78[g]/2)+((dhfgs78[g]+p)/2)*
					 (log(1+ delta[,g]/dhfgs78[g])))
		}
		num <- exp(log_num)
		zmat <- num/rowSums(num)
	}
	else{
		log_num <- matrix(0,n,G)
		for(g in 1:G){
			log_num[,g]<-log(pig[g])+lgamma((dhfgs78[g]+p)/2)-(1/2)*log(det(sigma[,,g]))-((p/2)*(log(pi)+log(dhfgs78[g]))
					+lgamma(dhfgs78[g]/2)+((dhfgs78[g]+p)/2)*(log(1+ delta[,g]/dhfgs78[g])))
		}	
		num <- exp(log_num)
		zmat <- num/rowSums(num)
		
	}
	if(clas>0){
		zmat <- unkno*zmat
		for(i in 1:n){
			if(kno[i]==1){
				zmat[i, known[i]] <- 1
			}
		}
	}
	zmat
}

