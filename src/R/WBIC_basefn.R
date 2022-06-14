WBIC <- function(llk){
  wbic <- - mean(rowSums(llk))
  return(wbic)
}