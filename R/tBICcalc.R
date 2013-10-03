tBICcalc <-
function(conv,G,p,mod,logl,n,gauss,univar,submod13){
	if(conv==0){
		BIC <- -Inf
	}
	if(conv==1){
		freepar <- G-1 + G*p
		if(univar){
			if(mod=="univUU"){
				freepar <- freepar + G + G
			}
			if(mod=="univUC"){
				freepar <- freepar + G + 1
			}
			if(mod=="univCU"){
				freepar <- freepar + 1 + G
			}
			if(mod=="univCC"){
				freepar <- freepar + 1 + 1
			}
		}
		if(submod13=="UUU"){
			freepar <- freepar + G*p*(p+1)/2
		}
		if(submod13=="CCC"){
			freepar <- freepar + p*(p+1)/2
		}
		if(submod13=="CUC"){
			freepar <- freepar + G*p*(p+1)/2 - (G-1)*p
		}
		if(submod13=="CUU"){
			freepar <- freepar + G*p*(p+1)/2 - (G-1)
		}
		if(submod13=="CIC"){
			freepar <- freepar + p
		}
		if(submod13=="CIU"){
			freepar <- freepar + G*p - (G-1)
		}
		if(submod13=="UIU"){
			freepar <- freepar + G*p
		}
		if(submod13=="UII"){
			freepar <- freepar + G
		}
		if(submod13=="CII"){
			freepar <- freepar + 1
		}
		#ITERATIVE PROCEDURE MODELS
		if(submod13=="UCC"){
			freepar <- freepar + p*(p+1)/2 + G-1
		}
		if(submod13=="UUC"){
			freepar <- freepar + G*p*(p+1)/2 - (G - 1)*(p-1)
		}
		if(submod13=="CCU"){
			freepar <- freepar + p*(p+1)/2 + (G - 1)*(p-1)
		}
		if(submod13=="UCU"){
			freepar <- freepar + p*(p+1)/2 + (G - 1)*p
		}
		if(submod13=="UIC"){
			freepar <- freepar + p + G-1
		}
		if(substring(mod,4,4)=="U"){
			freepar <- freepar + G
		}
		if(substring(mod,4,4)=="C"){
			freepar <- freepar + 1
		}
		if(gauss){freepar <- freepar-1}
		BIC <- 2 * max(logl) - freepar*log(n)
	}
	BIC
}

