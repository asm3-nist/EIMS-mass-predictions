maxPeak = function(x){ 
  
  a = strsplit(x[[1]][1],";")[[1]]
  b = length(a)
  mz = numeric(b);
  ab = numeric(b);
  
  for(i in 1:b){
    c = strsplit(a[i]," ")[[1]];
    k = which(c > 0)
      d = length(k);
      if(d>1){
        mz[i] = as.numeric(strsplit(a[i]," ")[[1]][k[1]]);
        ab[i] = as.numeric(strsplit(a[i]," ")[[1]][k[2]]);
      }
  }
  
  iter = numeric(0) 
  for(i in 1:length(mz)){
    if (is.na(mz[i]) == TRUE) { 
      iter = c(iter,i);
    }
  }
  if(length(iter)>0){
    mz = mz[-iter]
    ab = ab[-iter]
  }
  
  iter = numeric(0) 
  for(i in 1:length(mz)){
    if (mz[i] == 0) { 
      iter = c(iter,i);
    }
  }
  if (length(iter)>0){
    mz = mz[-iter];
    ab = ab[-iter];
  }
  
  peaks = cbind(mz,ab)

  lowest_ab = 3;
  
  numPeaks = dim(peaks)[1];
  
  # Check 1. Insufficient number of peaks
  if (length(numPeaks) == 0) { return(peaks[1]); }
  if (numPeaks <=1){return(peaks[1]); }
  
  return(mz[numPeaks]);

}