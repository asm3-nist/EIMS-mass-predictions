asm_TopPeakSelector <- function(x,MF1){ 

  a = strsplit(x[[1]][1],";")[[1]]
  b = length(a)
  mz = numeric(b);
  ab = numeric(b)
  
  for(i in 1:b){
    c = strsplit(a[i]," ")[[1]];
    k = which(c > 0)
    d = length(k);
    if(d>1){
      if (as.numeric(strsplit(a[i]," ")[[1]][k[1]])==0) next
      mz[i] = as.numeric(strsplit(a[i]," ")[[1]][k[1]]);
      ab[i] = as.numeric(strsplit(a[i]," ")[[1]][k[2]]);
    }
  }
  
  intermediateTable = as.data.table(cbind(mz,ab));
  intermediateTable = intermediateTable[order(-ab)];
  return(intermediateTable[ab>=MF1])
}



