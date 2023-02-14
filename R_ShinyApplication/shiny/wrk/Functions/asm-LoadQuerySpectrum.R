asm_load_QuerySpectrum <- function(InpDirName,values){
  
  data = readLines(InpDirName)
  values$qData = data;
  
  iName = grep("Name: ",data);
  if(length(iName[[1]])==0){
    qName = "Unknown";
  } else {
    qName = strsplit(data[iName],": ")[[1]][2];  
  }
  
  values$qname = qName
  
  sdata = data[-grep(":",data)]
  
  a = length(sdata);
  c = sdata[a];
  d = strsplit(c," ");
  
  while(length(d[[1]])==0){
    a = a - 1;
    c = sdata[a];
    d = strsplit(c," ");
  }
  
  sdata = sdata[1:a];
  if(length(grep(";",sdata[1]))==0)
  {
    tdata = cat(sdata[1],";",sep="");
    for(i in 2:a){
      tdata = paste(tdata," ",sdata[i],";",sep="")
    }
    sdata = tdata
  }
  
  q = asmPeakListCreator(sdata)
  values$qspectrum = q
  
  mq = maxPeak(q);
  
  dq = asm_TopPeakSelector(q,0)
  q = asm_SpectraPadder(dq,mq);
  
  values$qspectrum_padded = q;
  values$qspectrum_dt = dq;
  
  
  
  iMW = grep("MW: ",data);
  if(length(iMW)==0){
    values$mwinspectrum <- 0; 
  } else {
    values$mwinspectrum <- 1;
    values$rMW = as.numeric(strsplit(data[iMW],": ")[[1]][2])
  }
  

  
  return(1);
  
}