summary.teigen <- function(object, ...){
	x <- object
	    cat("---------- Summary for teigen ----------", "\n\n")
	if(x$G==x$iclresults$G & x$modelname==x$iclresults$modelname){
#		cat("     BIC and ICL select same model", "\n")
		cat("         ----   RESULTS   ----    ", "\n")
		cat("         Loglik:     ", x$logl, "\n")
		cat("         BIC:        ", x$bic, "\n")
		cat("         ICL:        ", x$iclresults$icl, "\n")
		cat("         Model:      ", x$modelname, "\n")
		cat("         # Groups:   ", x$G, "\n")
	}
	else{
#		cat("      BIC and ICL disagree...", "\n\n")
		
		cat("          ---- BIC RESULTS ----    ", "\n")
		cat("          Loglik:     ", x$logl, "\n")
		cat("          BIC:        ", x$bic, "\n")
		cat("          Model:      ", x$modelname, "\n")
		cat("          # Groups:   ", x$G, "\n\n")
		
		cat("          ---- ICL RESULTS ----    ", "\n")
		cat("          Loglik:     ", x$iclresults$logl, "\n")
		cat("          ICL:        ", x$iclresults$icl, "\n")
		cat("          Model:      ", x$iclresults$modelname, "\n")
		cat("          # Groups:   ", x$iclresults$G, "\n")
	}
}
