tdgupdate <-
function(p,G,mod,sg,lambdag,ng,ag,submod13){
	dg <- array(0, dim=c(p, p, G))
	negdet <- FALSE
	if(any(submod13==c("CUC","UUC"))){
		for(g in 1:G){
			dg[,,g] <- eigen(ng[g]*sg[,,g])$vectors
		}
	}
	else{
		if(submod13=="CUU"){
			for(g in 1:G){
				detngsg <- det(ng[g]*sg[,,g])
				if(detngsg<=0){
					negdet <- TRUE
					break
				}
				dg[,,g] <- eigen((ng[g]*sg[,,g])/(detngsg^(1/p)))$vectors
			}
		}
		else{
			if(submod13=="UCC"){
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
						dum3 <- eigen(dum2)$vectors
						dg[,,] <- dum3
					}
				}
				else{
					negdet <- TRUE
				}
			}
			else{
				dg[,,] <- diag(p)
			}
		}
	}
	if(negdet){
		dg <- FALSE
	}
	dg
}

