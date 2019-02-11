check_filter <-
function(ntseq,method){
  mySeq <- s2c(ntseq)
  mySeq <- toupper(mySeq)
  if (method == 'Tm_Wallace'){
    baseset <- c("A","B","C","D","G","H","I","K","M","N","R","S","T","V","W","Y")
  }else if (method == 'Tm_GC'){
    baseset <- c("A","B","C","D","G","H","I","K","M","N","R","S","T","V","W","X","Y")
  }else if(method == 'Tm_NN'){
    baseset = c('A','C','G','T','I')
  }else{
    stop("Only methods 'Tm_Wallace' or 'Tm_GC' or 'Tm_NN' is allowed")
  }
  finalSeq <- NULL
  #i='A'
  for(i in mySeq){
    if(i %in% baseset){
      finalSeq <- append(finalSeq,i)
    }
  }
  return(finalSeq)
}
