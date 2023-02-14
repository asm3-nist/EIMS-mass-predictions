rm(list=ls())

directory = "SWGDRUG3.4"
setwd(directory)
datafile = paste0(directory,".MSP")
data = readLines(datafile)

iNames = grep("Name: ",data)

for(i in 1:length(iNames)){
  if(i == length(iNames)){
    end = length(data);
  } else {
    end = iNames[i+1]-1
  }
  
  tempData = data[iNames[i]:end]
  id = strsplit(tempData[grep("ID: ",tempData)],": ")[[1]][2]
  filename = paste0(id,".MSP");
  j = 0;
  np = as.numeric(grep("Num peaks:",tempData));
  
  sink(filename)
  for(j in 1:length(tempData)){
    if(j > np){
      a = strsplit(tempData[j]," ")[[1]];
      if (length(a)>0){
        b = which(a != "")
        cat(a[b])
        cat(" \n") 
      }
    } else if(j == np){
      a = strsplit(tempData[j],": ")[[1]][2];
      cat(paste0("Num Peaks: ", a, "\n"))
    } else {
      cat(paste0(tempData[j],"\n"));  
    }
    
  }
  sink()
  
}

setwd("..")