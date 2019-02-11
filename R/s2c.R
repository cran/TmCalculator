s2c <-
function(strings){
  vecChar <- vector()
  for (i in nchar(strings):1){
    vecChar <- append(unlist(strsplit(strings,''))[i],vecChar)
  }
  return(vecChar)
}
