tfminup <-
function(mod,G,sg,dg,ag,p,ng,lambdag,submod13){
	fmin <- 0
	if(submod13=="UCC"){
		for(g in 1:G){
			fmin <- fmin + sum(diag((ng[g]*sg[,,g])%*%solve(dg[,,g]%*%ag[,,g]%*%t(dg[,,g]))))/lambdag[g] + p*ng[g]*log(lambdag[g])
		} 
	}
	if(submod13=="UUC"){
		for(g in 1:G){
			fmin <- fmin + sum(diag((ng[g]*sg[,,g])%*%solve(dg[,,g]%*%ag[,,g]%*%t(dg[,,g]))))/lambdag[g] + p*ng[g]*log(lambdag[g])
		} 
	}
	if(submod13=="CCU"){
		for(g in 1:G){
			fmin <- fmin + sum(diag(dg[,,g] %*% solve(ag[,,g]) %*% t(dg[,,g]) %*% (ng[g]*sg[,,g])))/lambdag[g] + p*ng[g]*log(lambdag[g])
		} 
	}
	fmin
}

