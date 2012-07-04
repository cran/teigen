tsigmaup <-
function(p,G,sg,lambdag,dg,ag,mod,univar,submod13){
	sigma <- array(0, dim=c(p, p, G))
	if(any(submod13==c("CCC","UUU")) | univar){
		sigma <- sg
	}
	else{
		if(substring(mod,2,2)=="I"){
			for(g in 1:G){
				sigma[,,g] <- lambdag[g] * ag[,,g]
			}
		}
		else{		
			for(g in 1:G){
				sigma[,,g] <- lambdag[g] * (dg[,,g] %*% ag[,,g] %*% t(dg[,,g]))
			}
		}
	}
	sigma
}

