modelgen <- function(){
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
  teigenModels
}
