teigen.parallel <- function(x, Gs=1:9, numcores=NULL, models="all", init="kmeans", scale=TRUE, dfstart=50, clas=0, known=NULL, training=NULL, gauss=FALSE, dfupdate=TRUE, eps=c(0.001,0.1), anneal=NULL, maxit=c(20,1000)){
	require("parallel")
	if(is.null(numcores)){
		numcores <- detectCores()
	}
	teigenModels <- list()
	teigenModels[["altnames"]] <- c("VVVV", "VVVE","EEVV", "EEVE","EVVV", "EVVE","EEEV","EEEE",
	                                "EVIV", "EVIE","EEIV", "EEIE","VIIV","VIIE","EIIV","EIIE",
	                                "VVIV","VVIE","VEEV", "VEEE","VEVV", "VEVE","VEIV", "VEIE",
	                                "VVEV","VVEE","EVEV","EVEE")
	teigenModels[["altunivariate"]] <- c("univVV", "univVE", "univEV", "univEE")
	teigenModels[["multivariate"]] <- c("UUUU", "UUUC","CUCU", "CUCC","CUUU", "CUUC","CCCU","CCCC",
	                                    "CIUU", "CIUC","CICU", "CICC","UIIU","UIIC","CIIU","CIIC",
	                                    "UIUU","UIUC","UCCU", "UCCC","UUCU", "UUCC","UICU", "UICC",
	                                    "UCUU","UCUC","CCUU","CCUC")
	teigenModels[["univariate"]] <- c("univUU", "univUC", "univCU", "univCC")
	teigenModels[["dfconstrained"]] <- c("UUUC","CUCC","CUUC","CCCC",
	                                     "CIUC", "CICC","UIIC","CIIC",
	                                     "UIUC","UCCC","UUCC", "UICC",
	                                     "UCUC","CCUC")
	teigenModels[["altdfconstrained"]] <- c("VVVE", "EEVE", "EVVE","EEEE",
	                                        "EVIE","EEIE","VIIE","EIIE",
	                                        "VVIE", "VEEE", "VEVE", "VEIE","VVEE","EVEE")
	teigenModels[["dfunconstrained"]] <- c("UUUU","CUCU","CUUU","CCCU",
	                                       "CIUU", "CICU","UIIU","CIIU",
	                                       "UIUU","UCCU","UUCU", "UICU",
	                                       "UCUU","CCUU")
	teigenModels[["altdfunconstrained"]] <- c("VVVV", "EEVV", "EVVV","EEEV",
	                                          "EVIV","EEIV","VIIV","EIIV",
	                                          "VVIV", "VEEV", "VEVV", "VEIV","VVEV","EVEV")
	if(length(models)==1){
		if(models=="dfunconstrained"){
			models <- teigenModels$dfunconstrained
		}
		else{
			if(models=="all"){
				if(ncol(x)==1){
					models <- teigenModels$univariate
				}
				else{
					models <- teigenModels$multivariate
				}
			}
			else{
				if(models=="gaussian"){
					models <- teigenModels$dfconstrained
				}
				else{
					if(models=="mclust"){
						models <- c("UUUC","CUCC","CCCC","CIUC","CICC","UIIC","CIIC","UIUC","UUCC","UICC")
					}
					else{
						if(models=="dfconstrained"){
							models <- teigenModels$dfconstrained
						}
						else{
							if(models=="univariate"){
								models <- teigenModels$univariate
							}
              else{
                if(models=="altall"){
                  if(ncol(x)==1){
                    models <- teigenModels$altunivariate
                  }
                  else{
                    models <- teigenModels$altnames
                  }
                }
              }
						}
					}
				}
			}
		}
	}
	modrep <- rep(models,length(Gs))
	backwards <- sort(Gs, decreasing=TRUE)
	grep <- rep(Gs, each=length(models))
	if(length(grep[grep==1])>0){
		mod1 <- NA
		cuont <- 1
		CCCCgroup <- c("CCCC", "CCCU", "CUCC", "CUCU", "CUUC","CUUU","UCCC","UCCU","UUCU","UUCC", "UUUC","UUUU")
    CCCCgroup <- c(CCCCgroup,teigenModels$altnames[teigenModels$multivariate%in%CCCCgroup])
		cccdum <- models[models %in% CCCCgroup]
		if(length(cccdum)>0){
			mod1[cuont] <- cccdum[1]
			cuont <- cuont+1
		}
		CICCgroup <- c("CICC","CICU","UICC","UICU","CIUC","CIUU","UIUC","UIUU")
		CICCgroup <- c(CICCgroup,teigenModels$altnames[teigenModels$multivariate%in%CICCgroup])
		cicdum <- models[models %in% CICCgroup]
		if(length(cicdum)>0){
			mod1[cuont] <- cicdum[1]
			cuont <- cuont+1
		}
		CIICgroup <- c("CIIC", "CIIU", "UIIC","UIIU")
		CIICgroup <- c(CIICgroup,teigenModels$altnames[teigenModels$multivariate%in%CIICgroup])
		ciidum <- models[models %in% CIICgroup]
		if(length(ciidum)>0){
			mod1[cuont] <- ciidum[1]
			cuont <- cuont+1
		}
		univgroup <- c("univCC","univCU","univUC","univUU","univVV", "univVE", "univEV", "univEE")
		unidum <- models[models %in% univgroup]
		if(length(unidum)>0){
			mod1[cuont] <- unidum[1]
			cuont <- cuont+1
		}
		modrep <- c(mod1, modrep[!grep==1])
		grep <- c(rep(1,length(mod1)), grep[!grep==1])
	}
	runvec <- 1:length(modrep)
	clus <- makeCluster(numcores)
	clusterEvalQ(clus, library(teigen))
	clusterExport(clus, ls(environment()), envir=environment())
  testparallel <- try(runlist <- clusterApplyLB(clus, runvec, function(g) teigen(x, grep[g], models=modrep[g], verbose=FALSE, init=init, scale=scale, dfstart=dfstart, clas=clas, known=known, training=training, gauss=gauss, dfupdate=dfupdate, eps=eps, anneal=anneal, maxit=maxit )), silent=TRUE)
  stopCluster(clus)
	if(class(testparallel)=="try-error"){
		stop(testparallel)
	}
	gmap <-  unlist(lapply(runlist, function(x) x$G))
	bicmap <-  unlist(lapply(runlist, function(x) x$bic))
	iclmap <-  unlist(lapply(runlist, function(x) x$iclresults$icl))
	modmap <-  unlist(lapply(runlist, function(x) x$modelname))
#	modmapicl <- unlist(lapply(runlist, function(x) x$iclresults$modelname))
	#dim(nodemap) <- c(3, length(models)*length(Gs))
#	nodemap <- t(nodemap)
	bestbic <- runlist[[which.max(bicmap)]]
	besticl <- runlist[[which.max(iclmap)]]
	bigbic <- matrix(-Inf, nrow=length(models), ncol=length(Gs))
	colnames(bigbic) <- paste("G=", Gs, sep="")
	rownames(bigbic) <- models
	bigicl <- bigbic
	for(i in 1:length(gmap)){
		bigbic[models==modmap[i],Gs==gmap[i]] <- bicmap[i]
		bigicl[models==modmap[i],Gs==gmap[i]] <- iclmap[i]
	}
	store <- list()
	store <- bestbic
#	store[["runlist"]] <- runlist
	store[["iclresults"]] <- besticl$iclresults
	store[["allbic"]] <- bigbic
	store$iclresults[["allicl"]] <- bigicl 
	store
}
