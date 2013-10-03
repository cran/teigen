tagupdate <-
function(p,G,mod,sg,lambdag,ng,n,dg,submod13){
	ag <- array(0, dim=c(p, p, G))
	negdet <- FALSE
	if(submod13=="CUC"){
		omegag <- matrix(0,p,p)
		for(g in 1:G){
			diag(omegag) <- diag(omegag) + eigen(ng[g]*sg[,,g])$values
		}
		ag[,,] <- omegag/(det(omegag)^(1/p))
	}
	if(submod13=="CUU"){
		for(g in 1:G){
			diag(ag[,,g]) <- eigen((ng[g]*sg[,,g])/(det(ng[g]*sg[,,g])^(1/p)))$values
		}
	}
	if(submod13=="CIC"){
		for(g in 1:G){
			diag(ag[,,g]) <- diag(n*sg[,,g])/(prod(diag(n*sg[,,g]))^(1/p))
		}
	}
	if(any(submod13==c("CIU","UIU"))){
		for(g in 1:G){
			diag(ag[,,g]) <- diag(ng[g]*sg[,,g])/(prod(diag(ng[g]*sg[,,g]))^(1/p))
		}
	}
	if(substring(mod,3,3)=="I"){
			ag[,,] <- diag(p)
	}
	if(submod13=="CCU"){
		for(g in 1:G){
			dwd <- diag(t(dg[,,g]) %*% (ng[g]*sg[,,g]) %*% dg[,,g])
			diag(ag[,,g]) <- dwd/(prod(dwd)^(1/p))
		}
	}
	if(submod13=="UCU"){
		for(g in 1:G){
			dwd <- diag(t(dg[,,g]) %*% (ng[g]*sg[,,g]) %*% dg[,,g])
			diag(ag[,,g]) <- dwd/(prod(dwd)^(1/p))
		}
	}
	#IP
	if(any(submod13==c("UCC","UUC"))){
		dum <- matrix(0,p,p)
		for(g in 1:G){
			dum <- dum + ng[g]*sg[,,g]/lambdag[g]
		}
		detdum <- det(dum)
		if(is.finite(detdum)){
			if(detdum<=0){
				negdet <- TRUE
			}
			else{
				dum2 <- dum/(detdum^(1/p))
				dum3 <- eigen(dum2)$values
				for(g in 1:G){
					diag(ag[,,g]) <- dum3
				}
			}
		}
		else{
			negdet <- TRUE
		}
	}
	if(submod13=="UIC"){
		dum <- matrix(0,p,p)
		for(g in 1:G){
			dum <- dum + ng[g]*sg[,,g]/lambdag[g]
		}
		dum2 <- diag(dum)/(prod(diag(dum))^(1/p))
		for(g in 1:G){
			diag(ag[,,g]) <- dum2
		}
	}
	if(negdet){
		ag <- FALSE
	}
	ag
}

