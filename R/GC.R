GC <-
function (ntseq, ambiguous = FALSE, totalnt = FALSE){
  if (length(ntseq) == 1 && is.na(ntseq)){
    return(NA)
  }
  if(class(ntseq)=="character"){
    ntseq <- toupper(ntseq)
    if(length(ntseq)==1 && nchar(ntseq)>1){
      vecSeq <- s2c(ntseq)
    }else if(length(ntseq) > 1){
      vecSeq <- ntseq
    }
  }else{
    stop("sequence is not characters")
  }
  nSeq <- length(vecSeq)
  if(!all(vecSeq %in% c("A","B","C","D","G","H","I","K","M","N","R","S","T","V","W","Y"))){
    warning("None Nucleic Acid Base in input Sequence")
  }
  nc <- sum(vecSeq %in% "C")
  ng <- sum(vecSeq %in% "G")
  na <- sum(vecSeq %in% "A")
  nt <- sum(vecSeq %in% "T")
  
  if(ambiguous == FALSE){
    ngc <- ng+nc
    nat <- na+nt
  }else{
    ngc <- ng+nc+sum(vecSeq %in% "S")
    nat <- na+nt+sum(vecSeq %in% "W")
    # for other ambiguous nucleatide acid base
    if (na+nc != 0) {     #M
      nm <- sum(vecSeq %in% "M")
      ngc <- ngc+nm*nc/(na+nc)
      nat <- nat+nm*na/(na+nc)
    }
    if (ng+nt != 0) {     #K
      nk <- sum(vecSeq %in% "K")
      ngc <- ngc+nk*ng/(ng+nt)
      nat <- nat+nk*nt/(ng+nt)
    }
    if (ng+na != 0) {    #R
      nr <- sum(vecSeq %in% "R")
      ngc <- ngc+nr*ng/(ng+na)
      nat <- nat+nr*na/(ng+na)
    }
    if (nc+nt != 0) {    #Y
      ny <- sum(vecSeq %in% "Y")
      ngc <- ngc+ny*nc/(nc+nt)
      nat <- nat+ny*nt/(nc+nt)
    }
    if (na+nc+ng != 0) {    #V
      nv <- sum(vecSeq %in% "V")
      ngc <- ngc+nv*(nc+ng)/(na+nc+ng)
      nat <- nat+nv*na/(na+nc+ng)
    }
    if (na+nc+nt != 0) {    #H
      nh <- sum(vecSeq %in% "H")
      ngc <- ngc+nh*nc/(na+nc+nt)
      nat <- nat+nh*(na+nt)/(na+nc+nt)
    }
    if (na+ng+nt != 0) {    #D
      nd <- sum(vecSeq %in% "D")
      ngc <- ngc+nd*ng/(na+ng+nt)
      nat <- nat+nd*(na+nt)/(na+ng+nt)
    }
    if (nc+ng+nt != 0) {    #B
      nb <- sum(vecSeq %in% "B")
      ngc <- ngc+nb*(nc+ng)/(nc+ng+nt)
      nat <- nat+nb*nt/(nc+ng+nt)
    }
  }
  if(totalnt){
    cat("argument totalnt is deprecated\n")
    ptGC <- 100*(ngc)/nSeq
    return(ptGC)
  }else{
    if (ngc+nat == 0) {
      ptGC <- NA
    }else {
      ptGC <- 100*ngc/(ngc+nat)
    }
  }
  return(ptGC)
}
