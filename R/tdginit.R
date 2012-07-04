tdginit <-
function(p,G,sg,mod,sgc){
	dg <- array(0, dim=c(p, p, G))
	if(substring(mod,3,3)=="U"){
		for(g in 1:G){	
			dg[,,g] <- eigen(sg[,,g])$vectors
		}
	}
	if(substring(mod,3,3)=="C"){
			dg[,,] <- eigen(sgc)$vectors
	}
	dg
}

