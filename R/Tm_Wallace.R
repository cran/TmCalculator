Tm_Wallace <-
function(ntseq,ambiguous=FALSE){
  mySeq <- check_filter(ntseq,method='Tm_Wallace')
  nSeq <- length(mySeq)
  nGC <- nSeq*GC(mySeq,ambiguous=ambiguous)/100
  nAT <- nSeq-nGC
  Tm <- 4*nGC+2*nAT
  return(Tm)
}
