tcontrandz <-
function(n,G){
	zmat <- matrix(runif(n*G,0,1), n, G)
	zmat <- zmat/rowSums(zmat)
	zmat
}

