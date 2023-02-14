asm_SpectraPadder <- function(a,max_mz){ 

  a_padded = numeric(max_mz);
  
  for(i in 1:dim(a)[1]){
    j = as.numeric(a[i,1]);
    a_padded[j] = as.numeric(a[i,2]);
  }
  
  return(a_padded);
}
