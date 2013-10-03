tlambdagupdate <-
function(G,mod,sg,dg,ag,p,n,ng,submod13){
	lambdag <- rep.int(0,G)
	if(submod13=="CUC"){
		omegag <- matrix(0,p,p)
		for(g in 1:G){
			diag(omegag) <- diag(omegag) + eigen(ng[g]*sg[,,g])$values
		}
		lambdag <- rep((det(omegag)^(1/p))/n,G)
	}
	if(submod13=="CUU"){
		thing <- 0
		for(g in 1:G){
			thing <- thing + (det(ng[g]*sg[,,g])^(1/p))/n
		}
		lambdag <- rep(thing,G)
	}
	if(submod13=="CIC"){
		for(g in 1:G){
			lambdag[g] <- (prod(diag(n*sg[,,g]))^(1/p))/n
		}
	}
	if(submod13=="CIU"){
		thing <- 0
		for(g in 1:G){
			thing <- thing + (prod(diag(ng[g]*sg[,,g]))^(1/p))/n
		}
		lambdag <- rep(thing,G)
	}
	if(submod13=="UIU"){
		for(g in 1:G){
			lambdag[g] <- (prod(diag(ng[g]*sg[,,g]))^(1/p))/ng[g]
		}
	}
	if(any(submod13==c("CII","UII"))){
		for(g in 1:G){
			lambdag[g] <- sum(diag(sg[,,g]))/p
		}
	}
	#IP
	if(submod13=="UCC"){
		for(g in 1:G){
			lambdag[g] <- sum(diag(sg[,,g] %*% solve(dg[,,g] %*% ag[,,g] %*% t(dg[,,g]))))/p
		}
	}
	if(submod13=="UUC"){
		for(g in 1:G){
			lambdag[g] <- sum(diag(sg[,,g] %*% dg[,,g] %*% diag(1/diag(ag[,,g])) %*% t(dg[,,g])))/p
		}
	}
	if(submod13=="UIC"){
		for(g in 1:G){
			lambdag[g] <- sum(diag(sg[,,g] %*% diag(1/diag(ag[,,g]))))/p
		}
	}
	if(submod13=="CCU"){
		for(g in 1:G){
			lambdag <- lambdag + sum(diag(ng[g]*sg[,,g] %*% dg[,,g] %*% diag(1/diag(ag[,,g])) %*% t(dg[,,g])))/(n*p)
		}
	}
	if(submod13=="UCU"){
		for(g in 1:G){
			lambdag[g] <- sum(diag(sg[,,g] %*% dg[,,g] %*% diag(1/diag(ag[,,g])) %*% t(dg[,,g])))/p
		}
	}
	lambdag
}

