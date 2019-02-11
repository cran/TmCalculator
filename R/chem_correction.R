chem_correction <-
function(Tm, DMSO=0, fmd=0, DMSOfactor=0.75,fmdfactor=0.65,fmdmethod="concentration",ptGC=NULL){
  if(!DMSOfactor %in% c(0.5,0.6,0.65,0.675,0.75)){
    stop("Only 0.5,0.6,0.65,0.675 and 0.75 are allowed for DMSOfactor")
  }
  if(!fmdfactor %in% c(0.6,0.65,0.72)){
    stop("Only 0.6,0.65 and 0.72 are allowed for fmdfactor")
  }
  if(is.null(DMSO)==FALSE){
    Tm <- Tm-DMSOfactor*DMSO
  }
  if(is.null(fmd)==FALSE){
    if(fmdmethod == "concentration"){
      Tm <- Tm-fmdfactor*fmd
    }else if(fmdmethod == "molar"){
      if(is.null(ptGC)==TRUE){
        stop("GC content in percent must be given for fmdmethod = molar")
      }
      Tm <- Tm+(0.453*(ptGC/100)-2.88)*fmd
    }else{
      stop("Only concentration and molar are allowed for fmdmethod")
    }
  }
  return(Tm)
}
