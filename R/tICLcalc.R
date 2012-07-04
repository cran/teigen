tICLcalc <-
function(conv,n,zmat,bic,modnum,G){
	if(conv==0){
		ICL <- -Inf
	}
	if(conv==1){
		ent <- apply(zmat,1,max)
		ICL <- bic[modnum,G] + 2*sum(log(ent))
	}
	ICL
}

