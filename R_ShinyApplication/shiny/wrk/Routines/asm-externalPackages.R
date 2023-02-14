if("devtools" %in% rownames(installed.packages()) == FALSE){
  install.packages("devtools")
  library(devtools)
} else { 
  library(devtools)
}

if("data.table" %in% rownames(installed.packages()) == FALSE) {
  install.packages("data.table")
  library(data.table)
} else {
  library(data.table)
}

if("DT" %in% row.names(installed.packages())==FALSE){
  devtools::install_github('rstudio/DT');
  library(DT);
} else {
  library(DT);
}

if("plotrix" %in% row.names(installed.packages())==FALSE){
  install.packages("plotrix");
  library(plotrix);
} else {
  library(plotrix);
}

if("scales" %in% row.names(installed.packages())==FALSE){
  install.packages("scales",dependencies=TRUE);
  library(scales)
} else {
  library(scales)
}

if("car" %in% row.names(installed.packages())==FALSE){
  install.packages("car",dependencies=TRUE);
  library(car)
} else {
  library(car)
}


if("shiny" %in% row.names(installed.packages())==FALSE){
  install.packages("shiny",dependencies=TRUE);
  require("shiny");
} else {
  require(shiny)
}
