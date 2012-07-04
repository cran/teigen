tsigmainvup <-
function(p,G,sigma){
	sigmainv <- array(0, dim=c(p, p, G))
	for(g in 1:G){
		sigmainv[,,g] <- solve(sigma[,,g])
	}
	sigmainv
}

