salt_correction <-
function(Na=0, K=0, Tris=0, Mg=0, dNTPs=0, method=1, ntseq=NULL,ambiguous=FALSE){
  if (method %in% c(5,6,7) && is.null(ntseq)==TRUE){
    stop("sequence is needed to calculate GC content and sequence length when method is 5,6,7")
  }
  if(method > 7){
    stop("allowed method is 1-7")
  }
  if(Na < 0 | K < 0 | Tris < 0 | Mg < 0 | dNTPs < 0){
    stop("all parameters 'Na','K','Tris','Mg','dNTP' should not be less than 0")
  }
  if (is.null(ntseq)==FALSE){
    corr = 0
    if(length(ntseq)==1){
      mySeq <- s2c(ntseq)
      mySeq <- toupper(mySeq)
    }else{
      mySeq <- toupper(ntseq)
    }
    
    nSeq <- length(mySeq)
    ptGC <- GC(mySeq,ambiguous=ambiguous)
    
    if(!method %in% c(1:7)){
      return(corr)
    }else{
      Mon = Na+K+Tris/2
      mg = Mg/1e+3
      
      if (sum(c(K,Mg,Tris,dNTPs)) > 0 && method != 7 && dNTPs < Mg){
        Mon <- Mon+120*sqrt(Mg-dNTPs)
      }
      mon <- Mon/1e+3
      if (method %in% c(1:7) && mon == 0){
        stop("total ion concentration of zero is not allowed in this method")
      }
      if (method == 1){
        corr <- 16.6*log10(mon)
      }else if (method == 2){
        corr <- 16.6*log10((mon)/(1.0+0.7*(mon)))
      }else if (method == 3){
        corr <- 12.5*log10(mon)
      }else if (method == 4){
        corr <- 11.7*log10(mon)
      }else if (method == 5){
        corr <- 0.368*(nSeq-1)*log(mon)
      }else if (method == 6){
        corr <- (4.29*ptGC/100-3.95)*1e-5*log(mon)+9.40e-6*log(mon) ^ 2
      }else if(method == 7){
        m7 <- c(3.92, -0.911, 6.26, 1.42, -48.2, 52.5, 8.31)
        dntps <- dNTPs*1e-3
        ka = 3e4
        mg <- (sqrt((ka*dntps-ka*mg+1)**2+4*ka*mg)-(ka*dntps-ka*mg+1))/(2*ka)
        R <- if (Mon > 0) sqrt(mg)/mon
        if (R < 0.22){
          corr <- (4.29*ptGC/100-3.95)*1e-5*log(mon)+9.40e-6*log(mon)**2
        }else if (R >= 0.22 && R < 6.0){
          m7[1] <- 3.92*(0.843-0.352*sqrt(mon)*log(mon))
          m7[4] <- 1.42*(1.279-4.03e-3*log(mon)-8.03e-3*log(mon)**2)
          m7[7] <- 8.31*(0.486-0.258*log(mon)+5.25e-3*log(mon)**3)
          corr <- (m7[1]+m7[2]*log(mg)+(ptGC/100)*(m7[3]+m7[4]*log(mg))+(1/(2.0*(nSeq-1))) *(m7[5]+m7[6]*log(mg)+m7[7]*log(mg)**2))*1e-5
        }
      }
    }
    return(corr)
  }
}
