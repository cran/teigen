estimateTime <- function(stage, start.time, totmod){
  if(stage=="init"){
    curwidth <- getOption("width")
    if(curwidth>=80){
      cat("Time taken:", format("???", width=6, justify="right"), "  |  (Approx) remaining:", format("???", width=6, justify="right"), "  |  " ,format("???", width=4, justify="right"), "% complete", "\r", sep="")
    }
    else{
      if(curwidth>=50){
        cat("(Approx) remaining:", format("???", width=6, justify="right"), "  |  " ,format("???", width=4, justify="right"), "% complete", "\r", sep="")
      }
      else{
        if(curwidth>=16){
          cat(format("0", width=4, justify="right"), "% complete", "\r", sep="")
        }
      }
    }
    flush.console()
  }
  else{
    modelcount <- stage+1
    modsleft <- totmod-modelcount
    timerun <- difftime(Sys.time(),start.time, units="secs")
    timeremain <- (timerun/modelcount)*modsleft
    if(timeremain>60){
      if(timeremain>3600){
        if(timeremain>86400){
          if(timeremain>604800){
            units(timeremain) <- "weeks"
          }
          else{
            units(timeremain) <- "days"
          }
        }
        else{
          units(timeremain) <- "hours"
        }
      }
      else{
        units(timeremain) <- "mins"
      }
    }
    curwidth <- getOption("width")
    if(curwidth>=80){
      cat("Time taken:", format(round(timerun,1), width=6, justify="right"), "  |  (Approx) remaining:", format(round(timeremain,1), width=6, justify="right"), "  |  " ,format(round((1-modsleft/totmod)*100), width=4, justify="right"), "% complete", "\r", sep="")
    }
    else{
      if(curwidth>=50){
        cat("Approx. remaining:", format(round(timeremain,1), width=6, justify="right"), "  |  " ,format(round((1-modsleft/totmod)*100), width=4, justify="right"), "% complete", "\r", sep="")
      }
      else{
        if(curwidth>=16){
          cat(format(round((1-modsleft/totmod)*100), width=4, justify="right"), "% complete", "\r", sep="")
        }
      }
    }
    flush.console()
    modelcount
  }
}