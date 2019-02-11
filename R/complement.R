complement <-
function(ntseq,reverse=FALSE){
  complement_table <- matrix(c('A','T','B','V','C','G','D','H','G','C','H','D','M','K','N','N','R','Y','S',
                               'S','T','A','U','A','V','B','W','W','X','X','Y','R','a','t','b','v','c','g',
                               'd','h','g','c','h','d','m','k','n','n','r','y','s','s','t','a','u','a','v',
                               'b','w','w','x','x','y','r'),ncol = 2, byrow = TRUE)
  new_seq <- NULL
  lineseq <- s2c(ntseq)
  for (i in lineseq){
    if(i %in% complement_table){
      n <- which(complement_table[,1] %in% i)
      new_seq <- append(new_seq,complement_table[n,2])
    }
  }
  if(reverse==TRUE){
    new_seq <- rev(new_seq)
  }
  rc_seq <- c2s(new_seq)
  return(rc_seq)
}
