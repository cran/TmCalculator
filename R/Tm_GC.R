Tm_GC <-
function(ntseq,ambiguous=FALSE,userset=NULL,variant="Primer3Plus",Na=50,
                  K=0,Tris=0, Mg=0, dNTPs=0, saltcorr=0, mismatch=TRUE){
  if(saltcorr == 5){
    stop('salt-correction method 5 not applicable to Tm_GC')
  }
  mySeq <- check_filter(ntseq,method='Tm_GC')
  nSeq <- length(mySeq)
  ptGC <- GC(mySeq,ambiguous=ambiguous)
  
  varTab <- matrix(c(69.3,0.41,650,1,0,81.5,0.41,675,1,0,81.5,0.41,675,1,2,81.5,0.41,500,1,3,78.0,0.7,500,1,3,
                     67.0,0.8,500,1,3,81.5,0.41,600,1,2,77.1,0.41,528,1,4),nrow=8,byrow=TRUE)
  rownames(varTab) <- c("Chester1993","QuikChange","Schildkraut1965","Wetmur1991_MELTING","Wetmur1991_RNA","Wetmur1991_RNA/DNA","Primer3Plus","vonAhsen2001")
  colnames(varTab) <- c("A","B","C","D","saltcorr")

  if(is.null(userset)){
    if(!variant %in% rownames(varTab)){
      stop("only Chester1993, QuikChange, Schildkraut1965, Wetmur1991_MELTING, Wetmur1991_RNA, Wetmur1991_RNA/DNA, Primer3Plus and vonAhsen2001 are allowed in variant")
    }else{
      gcCoef <- varTab[variant,]
    }
  }else{
    gcCoef <- as.numeric(userset)
  }
  saltcorr <- gcCoef[5]
  Tm = gcCoef[1]+gcCoef[2]*ptGC-gcCoef[3]/nSeq
  if(saltcorr !=0 ){
    tmp1 <- salt_correction(Na=Na,K=K,Tris=Tris,Mg=Mg,dNTPs=dNTPs,method=saltcorr,ntseq=mySeq)
    Tm <- Tm+tmp1
  }
  if(mismatch == TRUE){
    Tm <- Tm-gcCoef[4]*(sum(mySeq %in% 'X')*100/nSeq)
  }
  return(as.numeric(Tm))
}
