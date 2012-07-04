tlambdaginit <-
function(p,G,sg,mod){
	lambdag <- rep.int(0,G)
	if(substring(mod,1,1)=="U"){
		for(g in 1:G){
			lambdag[g] <- det(sg[,,g])^(1/p)
		}
	}
	if(substring(mod,1,1)=="C"){
		thing <- 0
		for(g in 1:G){
			thing <- thing + det(sg[,,g])^(1/p)
		}
		lambdag <- rep(thing,G)
	}
	lambdag
}

