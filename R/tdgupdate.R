tdgupdate <-
function(p,G,mod,sg,lambdag,ng,ag,submod13, dg, flipswitch=NULL){
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
				if(any(submod13==c("UCU","CCU"))){
					if(flipswitch){
					#MM1
						dum <- matrix(0,p,p)
						ainv <- dum
						for(g in 1:G){
							wg <- ng[g]*sg[,,g]
							wk <- eigen(wg, symmetric=TRUE, only.values=TRUE)$values[1]
							diag(ainv) <- diag(1/(lambdag[g]*ag[,,g]))
							dum <- dum + ainv %*% t(dg[,,g]) %*% wg - wk*ainv %*% t(dg[,,g])
						}
						#print(dum)
						if(all(is.finite(dum))){
							svdum <- svd(dum)
							dg[,,] <- svdum$v %*% t(svdum$u)
						}
						else{
							negdet <- TRUE
						}
					}
					else{
					#MM2
						dum <- matrix(0,p,p)
						ainv <- dum
						for(g in 1:G){
							wg <- ng[g]*sg[,,g]
							diag(ainv) <- diag(1/(lambdag[g]*ag[,,g]))
							#dum <- dum + wg %*% dg[,,g] %*% ainv  - (lambdag[g]*ag[p,p,g])*wg %*% dg[,,g]
							dum <- dum + wg %*% dg[,,g] %*% ainv  - (max(ainv))*wg %*% dg[,,g]
						}
						svdum <- svd(dum)
						dg[,,] <- svdum$v %*% t(svdum$u)
					
					}
				}
				else{
					dg[,,] <- diag(p)
				}
			}
		}
	}
	if(negdet){
		dg <- FALSE
	}
	dg
}

