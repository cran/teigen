taginit <-
function(p,G,sg,mod,sgc){
	ag <- array(0, dim=c(p, p, G))
	if(substring(mod,3,3)=="U"){
		for(g in 1:G){	
			diag(ag[,,g]) <- (eigen(sg[,,g])$values)/(det(sg[,,g])^(1/p))
		}
	}
	if(substring(mod,3,3)=="C"){
		for(g in 1:G){	
			diag(ag[,,g]) <- (eigen(sgc)$values)/(det(sgc)^(1/p))
		}
	}
	ag
}

