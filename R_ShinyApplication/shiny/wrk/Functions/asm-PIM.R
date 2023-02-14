PIM = function(x,ab_M){ 
  
  a = strsplit(x[[1]][1],";")[[1]]
  b = length(a)
  mz = numeric(b);
  ab = numeric(b);

  cl.br.StopCrit = 0.5;
  
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
  
  basePeak = max(as.numeric(ab));

  lowest_ab = (ab_M/100)*basePeak;
  
  numPeaks = dim(peaks)[1];
  
  if (length(numPeaks) == 0) { return(peaks[1]); }
  
  if (numPeaks <=1){return(peaks[1]); }
  
  for(i in numPeaks:2){

    m0 = peaks[i,1];
    a0 = peaks[i,2];
    
    m1 = peaks[(i-1),1];
    md1 = m0 - m1;
    a1 = peaks[(i-1),2];
    
    if((i-2) < 1){
      md2 = 0;
    } else{
      m2 = peaks[(i-2),1];
      md2 = m0 - m2;  
      a2 = peaks[(i-2),2];
    }
    
    if (a0 < lowest_ab) next
    
    
    if (md1 > 2){ return(m0)} # no peaks within 2 Da
    
    
    if (md1 == 2){ # possible Cl Br peaks
      
      if (a1 < as.integer(cl.br.StopCrit*a0)){ return(m0)}
      
    } else if (md1 == 1){ # check for c13 isotope peaks
      
      isocalc = (a1*0.011*m1) / 14;
      
        if(as.integer(3*isocalc) > a0 | (a0 - as.integer(isocalc)) < lowest_ab){
          if(i==2){ 
            return(m1);
          }
          next;
        }
      
        if(md2 > 2){
          return(m0);
        } else {
            if (i-2 >= 1){
              if (a2 < as.integer(cl.br.StopCrit*a0)){
              return(m0);
              }
            }
        }
    } # end else if md1==1
  } # end for loop
  
  return(peaks[numPeaks,1]);
}