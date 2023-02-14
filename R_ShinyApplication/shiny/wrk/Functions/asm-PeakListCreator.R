asmPeakListCreator <- function(x){ 

  y = length(x);
  z = "";
  for(j in 1:y){
    z = paste(z,x[j],sep="")
  }
  
    
  return(z);
}
