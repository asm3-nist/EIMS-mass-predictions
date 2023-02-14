# Obtaining molecular mass from EI MS
#
# Arun S. Moorthy, arun.moorthy@nist.gov
# 2020
# =======================================================

rm(list=ls())

# ==============================================================================
# Arun's Personal Functions and Routines
# ==============================================================================
source("Routines/asm-externalPackages.R")
source("Routines/asm-colours.R")

system("cmd.exe",input="dir SearchLibraries /s /b /o:n /ad > Data/libDirectories.txt");
LibDirectories = readLines("Data\\libDirectories.txt");
nLibraries = length(LibDirectories);
libNames = array("a",dim=c(nLibraries));

for(i in 1:nLibraries){
  b = strsplit(LibDirectories[i],"\\\\");
  c = length(b[[1]]);
  libNames[i] = b[[1]][c];
}

LossProbs = as.data.table(read.csv("Data\\MainLosses_prob.csv")[,2:3]);
colnames(LossProbs)[1] = "Loss"
tau_P = 0.06;
NegLoss = LossProbs[MainLosses<=tau_P]
beta = 5;

source("Functions/asm-LoadQuerySpectrum.R")
source("Functions/asm-PeakListCreator.R")
source("Functions/asm-maxPeakSelector.R")
source("Functions/asm-TopPeaksSelector.R")
source("Functions/asm-SpectraPadder.R")
source("Functions/asm-PIM.R")


stime <- system.time({
  SWG33Library  <- readRDS("LibraryCorrections/SWGDRUG33.RDS")
})[3]
ReadInDTLibraryTime_SWG33 = stime

nSpecSWG33 = dim(SWG33Library)[1];

stime <- system.time({
  SWG34Library  <- readRDS("LibraryCorrections/SWGDRUG34.RDS")
})[3]
ReadInDTLibraryTime_SWG34 = stime

nSpecSWG34 = dim(SWG34Library)[1];

# Define UI for data upload app ----
## =============================================================================
ui <- fluidPage(
 
  theme="bootstrap.css",
  
  # App title ----
  titlePanel("Determining Molecular Mass from EI Mass Spectra"),
  
  # Horizontal line ----
  tags$hr(),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      
      tabsetPanel(
        
        tabPanel("Main",
                 br(),
      
          # Input: Select a file ----
            fileInput("file", "Choose Query Spectrum (MSP File)",
                    multiple = FALSE,
                    accept = c(".MSP")),
      
          # Input: Select a library ----
            selectInput("lib","Select Search Library",libNames),
      
          # Input: Expected Mass Range ---
            sliderInput("aThetas", label = "Expected Mass Range", min = 50, 
                    max = 700, value = c(150, 250)),
      
          # Input: Prediction Methods
            checkboxGroupInput("Methods", label = "Methods", 
                            choices = list("Peak Interpretation Method (PIM)" = 1, 
                                           "Simple Search Hitlist Method (SS-HM)" = 2,
                                           "Iterative Hybrid Search Hitlist Method (iHS-HM)" = 3),
                            selected = c(1,2,3)),
     
      
          # Horizontal line ----
            tags$hr(),
      
          # Compute Mass Estimate ----
            actionButton("compute", "Predict Molecular Mass"), 
      
          # Horizontal line ----
            tags$hr(),
      
      
            div(p(style="border: 1px solid gray"))
      
          ),
        
          tabPanel("Advanced Options",
             br(),
             
             tabsetPanel(
               tabPanel("PIM",
                        br(),
                        div(p("The molecular-ion (M) abundance threshold is the minimum abundance (as percentage of base peak) acceptable for a predicted M. The default value was chosen by experience."), style = "color:black;font-size:8pt"),
                        # Input: M threshold (0.3)
                        sliderInput("ab_M", label = "M abundance threshold", min = 0.1, 
                                    max = 100, value = 0.3,step=0.1,ticks=FALSE),
                        
                        br(),
                        div(p("The noise tolerance is the maximum abundance (as percentage of base peak) acceptable for an illogical peak. The default value was chosen by experience."), style = "color:black;font-size:8pt"),
                        # Input: a threshold (0.1)
                        sliderInput("ab_a", label = "abundance threshold", min = 0.0, 
                                    max = 100, value = 0.0, step=0.1,ticks=FALSE),
                        
                        br(),
                        div(p("Classification threshold values are used to decide if a given molecular mass prediction is likely to be correct. The default value was chosen because it optimized reciever operator characteristics for a test set of drug spectra."), style = "color:black;font-size:8pt"),
                        # Input: PIM classification threshold (0.07)
                        sliderInput("tau_g", label = "classification threshold", min = 0.01, 
                                    max = 1, value = 0.01,ticks=FALSE)
                        ),
               
               tabPanel("SS-HM",
                        br(),
                        div(p("The Hitlist depth value is the number of reference library spectra consider for identifying potential molecular masses of the query. The default value was chosen by experience."), style = "color:black;font-size:8pt"),
                        # Input: M threshold (0.3)
                        sliderInput("hl_D", label = "Hitlist Depth", min = 2, 
                                    max = 50, value = 25,step=1,ticks=FALSE),
                        
                        br(),
                        div(p("The match factor significance value is the difference between a pair of match factors that would be be deemed signficant. The default value was chosen by experience."), style = "color:black;font-size:8pt"),
                        # Input: M threshold (0.3)
                        sliderInput("B_SS", label = "SS-exponent", min = 60, 
                                    max = 100, value = 75,step=1,ticks=FALSE),
                        
                        br(),
                        div(p("Classification threshold values are used to decide if a given molecular mass prediction is likely to be correct. The default value was chosen because it optimized reciever operator characteristics for a test set of drug spectra."), style = "color:black;font-size:8pt"),
                        # Input: SS-HM classification threshold (0.39)
                        sliderInput("tau_s", label = "classification threshold", min = 0.01, 
                                    max = 1, value = 0.65,ticks=FALSE)
                        ),
               
               tabPanel("iHS-HM",
                        br(),
                        div(p("The minimum expected match factor is lowest match factor for which two molecules may be considered isospectral. The default value was chosen by experience."), style = "color:black;font-size:8pt"),
                        # Input: M threshold (0.3)
                        sliderInput("mEMF", label = "minimum expected match factor", min = 500, 
                                    max = 800, value = 700,step=1,ticks=FALSE),
                        
                        
                        br(),
                        div(p("Classification threshold values are used to decide if a given molecular mass prediction is likely to be correct. The default value was chosen because it optimized reciever operator characteristics for a test set of drug spectra."), style = "color:black;font-size:8pt"),
                        # Input: iHS-HM classification threshold (0.05)
                        sliderInput("tau_t", label = "classification threshold", min = 0.01, 
                                    max = 1, value = 0.01,ticks=FALSE)
                        )
               
             )
             
             
             

             
          )
      )
    ),
    
    
    
    
# Main panel for displaying outputs ----
    mainPanel(
      
      plotOutput("QuerySpectrum"),
      
      tabsetPanel(
        
        tabPanel("Results Summary",
     
                  div(p(style="border: 5px solid gray")),
                  htmlOutput("rMW",style="color:black;font-size:12pt"),
                  htmlOutput("gamma",style="color:black;font-size:12pt"),
                  htmlOutput("sigma",style="color:black;font-size:12pt"),
                  htmlOutput("theta",style="color:black;font-size:12pt"),
                  br()
                  

                 
        ),
        
        tabPanel("PIM",
                 
                 div(p(style="border: 5px solid gray")),
                 htmlOutput("PIM",style="color:black;font-size:12pt"),
                 br()
                 

        ),
        
        tabPanel("SS-HM",
                 
                 div(p(style="border: 5px solid gray")),
                 htmlOutput("SSHM",style="color:black;font-size:12pt"),
                 br()
                 

        ),
        
        tabPanel("iHS-HM",
                 
                 div(p(style="border: 5px solid gray")),
                 plotOutput("OmegaX"),
                 #htmlOutput("HSHM",style="color:black;font-size:12pt"),
                 br()
                 
                 # div(p(style="border: 5px solid gray")),
                 # 
                 # 
                 # plotOutput("OmegaX"),
                 # div(p(textOutput("MW")),style = "color:black;font-size:14pt"),
                 # div(p(textOutput("Confidence1")),style= "color:black;font-size:14pt")
        )
        
      ),

      div(p(style="border: 5px solid gray")),
      div(p("DISCLAIMER: The Molecular Mass Predictor (MMP) is a prototype implementation designed for demonstration purposes. The authors of the MMP cannot
      guarantee the accuracy of results generated using the MMP, and cannot validate
      claims of others using this software."),
      style = "color:black;font-size:12pt")

    )
    
  )
)


## =============================================================================









# Define server logic to read selected file ----
## =============================================================================
server <- function(input, output,session) {
  
  
    values <- reactiveValues()

      values$mwinspectrum <- 0;
      values$rMW <- 0;
  
      # PIM
      values$isgamma <- 0;
      values$gamma <- 0; 
      values$Igamma <- 0;
      values$Iratio <- 0;
      values$elosses <- NULL;
      values$elosses2 <- 0;
      
      # SS-HM
      values$issigma <- 0;
      values$sigma <- 0;
      values$Isigma <- 0;

      values$SigmaFullResultMatrix <- NULL;
    

      # HS-HM
      values$istheta <- 0;
      values$theta <- 0;
      values$Itheta <- 0;
      values$AllThetas <- NULL;
      values$AllOmegas <- NULL;
    
      values$is_qspectrum = 0;
    
      values$qData <- NULL
      values$qspectrum <- NULL
      values$qspectrum_dt <- NULL;
      values$qspectrum_padded <- seq(100,1,-1)
      values$qname <- "Random";
      
      values$elosses <- NULL
    
   
    
    
     
    trigger1 <- eventReactive(input$file$datapath,ignoreInit = F, {

      InpDirName = input$file$datapath
      a = asm_load_QuerySpectrum(InpDirName,values);
      
      values$is_qspectrum = a;
      values$isgamma <- 0;
      values$issigma <- 0;
      values$issigma_cor <-0;
      values$istheta <- 0;
      
      return(0);
      
    })
    
    output$QuerySpectrum <- renderPlot({
      
      a = trigger1();
      
      q = values$qspectrum_padded
      name = paste("Query Spectrum: ",values$qname);
      plot(q,type="h",xlab="mass",ylab="abundance");
      title(name,adj=0);
      
    })
    
    output$rMW <- renderText({
      trigger1()
      if(values$mwinspectrum==1){
        message = paste("MW from meta data: ",values$rMW, "Da");
      } else {
        message = paste("No MW in meta data.")
      }
      return(message)
    })

    
        
    
    
    trigger2 <- eventReactive(input$compute,{
      
      if(values$is_qspectrum == 0){
        return(0)
      } else {
      
        whichPreds = input$Methods;
      
      
          if (1 %in% whichPreds){  # COMPUTE PIM Prediction
            values$isgamma = 1;
            
            a = values$qspectrum_dt;  # local copy of the spectrum as data.table
            #print(a)
            
            dq = values$qspectrum;    # local copy of spectrum as vector
            #cat(dq)
            
            ab_M = input$ab_M;        # minimum abundance of acceptable molecular ion
            
            g = PIM(dq,ab_M);         # compute interpretation based prediction
            # cat(g)
             
            pAb = a[mz==g,ab];        # abundance of potential molecular ion
            bpAb = max(a[,ab]);       # abundance of base peak
            # cat(pAb)
            
            lambda = g - as.numeric(a[,mz])   # estimated neutral losses given g
            b = cbind(a,lambda)               # data table with mass, ab, neutral loss
            # print(b)
            lambda = lambda[lambda<43];
            #print(lambda)
            
            
            U = NegLoss[Loss%in%lambda,Loss]; # set of "unexpected" losses
            # print(length(U))

            if (length(U)>0){
              abU = b[lambda%in%U,ab]           # abundance of unexpected losses
              # abU = abU[abU>=(bpAb*input$ab_a)/100];
              # cat(abU)
              abU = abU-(bpAb*input$ab_a/100)
              abU[abU<0] = 0
              # cat(abU)

              U2 = length(abU[abU>0]);
              
              x = sum(abU)/(pAb)                # ratio of unexpected losses to ab of M
              
              I = 2 - 2/(1+exp(-beta*x));       # gamma classifier
              
            } else { # no unexpected losses
              I = 1; # gamma classifier
              x = 0;
              U2 = 0;
            }
            
            
            values$gamma = g; 
            values$Igamma = I;
            values$Iratio = x;
            values$elosses = U;
            values$elosses2 = U2;
            
          } else {
            values$isgamma = 0;
          }
        
        
        
      
          if (2 %in% whichPreds){  # COMPUTE SS-HM Prediction
            values$issigma = 1;
            
            a = values$qspectrum_dt;  # local copy of the spectrum as data.table
            dq = values$qspectrum;    # local copy of spectrum as vector
            ab_M = input$ab_M;        # minimum abundance of acceptable molecular ion
            g1 = PIM(dq,ab_M);         # compute interpretation based prediction
            
            # this is where sigma must be computed
            data = values$qData;
            a = length(data);
            c = data[a];
            d = strsplit(c," ");

              while(length(d[[1]])==0){
                a = a -1;
                c = data[a];
                d = strsplit(c," ");
              }

              e = length(d[[1]])-1;
              f = d[[1]][e];

              sink("Data\\Query_SS.MSP")
                cat(data[1:a],sep="\n")
              sink();
              
             g = input$lib;
             h = which(libNames == g);
             LibDirName = LibDirectories[h];

                if(g == "mainlib"){
                  libType = " /MAIN ";
                } else if (g == "replib"){
                  libType = " /REPL ";
                } else {
                  libType = " /LIB ";
                }
                        cat(g)
                        #NEW in MMP2020
                        if(g=="SWGDRUG3.3"){
                          RefLibrary = SWG33Library;
                        } else if (g == "SWGDRUG3.4"){
                          RefLibrary = SWG34Library;
                        }
              

                y = "Fast";
                if(y == "Default"){
                  presearch = " Sdva^\n";
                } else if (y == "Fast"){
                  presearch = " Sf^\n";
                } else {
                  presearch = " Ssva^\n";
                }
            
            
            
              BATFILE_NAME = "Data\\asm-SS-HM-BAT.bat"
              ProgName = "MSPepSearch\\2017_05_15_MSPepSearch\\x64\\MSPepSearch64.exe";
              Output = "Data\\HitList_SS.txt";
              nHits = input$hl_D;  # the number of hits stored (from SS_ESTMW.C)
              B_SS = input$B_SS;

                sink(BATFILE_NAME)
                cat(paste(ProgName,presearch,sep=""))
                cat(paste(libType,LibDirName,"^\n", sep=""))
                cat(paste(" /HITS ", nHits ,"^\n",sep=""))
                cat(paste(" /INP Data\\Query_SS.MSP ^\n",sep=""));
                cat(" /OutNumMP^\n")
                cat(" /OutMW^\n")
                cat(" /PROGRESS^\n")
                cat(" /LibInMem^\n")
                cat(paste(" /OUTTAB ", Output,sep=""))
                sink()

                withProgress(message = "completing SS-HM", value = 0.3, {
                  system2(BATFILE_NAME);

                })
            
                # generate prediction
                fileName = Output
                cat(fileName); cat("\n")
                Results = readLines(fileName);
                Results_w = Results[5:(length(Results)-1)];
                
                ListMass = c();
                ListCorrection = c();
                ListMF = c();
                
                ku = 1+nHits;
                
                ku = min(ku,length(Results_w));
                
                PartResult = Results_w[1:ku]
                
                for(j in 1:ku){
                  
                  rank = as.numeric(strsplit(PartResult[j],"\t")[[1]][2]);
                  
                  if(is.na(rank)) break
                  
                  if(rank != j) break 
                  
                  HitName = strsplit(PartResult[j],"\t")[[1]][9] 
                  k = which(RefLibrary[,Name]==HitName)
                  hgamma = PIM(RefLibrary[k,PeakList],ab_M)
                  hrMW = as.numeric(RefLibrary[k,MW])
                  HitCorrection = hrMW - hgamma;
                  
                  ListCorrection = c(ListCorrection,HitCorrection)
                  
                  HitMass = as.numeric(strsplit(PartResult[j],"\t")[[1]][6]);
                  ListMass = c(ListMass,HitMass)
                  
                  HitMF = as.numeric(strsplit(PartResult[j],"\t")[[1]][7]);
                  ListMF = c(ListMF,HitMF)
                
                }
                
                ListCorrection = as.numeric(ListCorrection) + g1;
                sigmas = unique(ListCorrection);
                probSigmas = numeric(length(sigmas));
                
                
                FullResultMatrix_rows = length(sigmas);
                sigma_freq = numeric(length(sigmas));
                FullResultMatrix_cols = max(table(ListCorrection));
                FullResultMatrix = array(0,dim=c(FullResultMatrix_rows,FullResultMatrix_cols));
                
                for(j in 1:length(sigmas)){
                  k = which(ListCorrection==sigmas[j]);
                  num = 0;
                  sigma_freq[j] = length(k);
                  
                  for(l in 1:length(k)){
                    num = num + 2^((ListMF[k[l]] - ListMF[1])/B_SS)
                    #num = num + 2^((ListMF[k[l]] - 999)/B_SS)
                    FullResultMatrix[j,l] = ListMF[k[l]]
                  }
                  probSigmas[j] = num;
                }
                
                den=0;
                for(j in 1:length(ListMF)){
                  den = den + 2^((ListMF[j] - ListMF[1])/B_SS)
                  #den = den + 2^((ListMF[j] - 999)/B_SS)
                }
                
                probSigmas = round((probSigmas/den),2);
                
                sigma_star = sigmas[which(probSigmas==max(probSigmas))[1]];
                prob_sigma_star= max(probSigmas);
                
            
            complete_table = cbind(sigmas,sigma_freq,probSigmas,FullResultMatrix);
                
                
            values$SigmaFullResultMatrix = complete_table;        
            values$sigma = sigma_star; 
            values$Isigma = prob_sigma_star;
            
            
          } else {
            values$issigma = 0;
          }
      
        
        
        
          if (3 %in% whichPreds){
            values$istheta = 1;
            
            # this is where theta must be computed
              data = values$qData;

              a = length(data);
              c = data[a];
              d = strsplit(c," ");

              while(length(d[[1]])==0){
                a = a -1;
                c = data[a];
                d = strsplit(c," ");
              }

              e = length(d[[1]])-1;
              f = d[[1]][e];

              lMass = input$aThetas[1];
              uMass = input$aThetas[2];
              NMEs = seq(lMass,uMass,1);

              m = grep("MW",data);

              sink("Data\\Query_HS.MSP")
              if(length(m) == 0){
                for(j in 1:length(NMEs)){
                  cat(paste(data[1],"\n",sep=""))
                  cat(paste("MW: ",NMEs[j],"\n",sep=""))
                  for(i in 2:a){
                      cat(paste(data[i],"\n",sep=""));
                  }
                  cat("\n\n")
                }
              } else {
                for(j in 1:length(NMEs)){
                  for(i in 1:a){
                    if(length(grep("MW: ",data[i])) == 1){
                      cat(paste("MW: ",NMEs[j],"\n",sep=""))
                    } else {
                      cat(paste(data[i],"\n",sep=""));
                    }
                  }
                  cat("\n\n")
                }
              }
              sink();

              g = input$lib;
              h = which(libNames == g);
              LibDirName = LibDirectories[h];

              if(g == "mainlib"){
                libType = " /MAIN ";
              } else if (g == "replib"){
                libType = " /REPL ";
              } else {
                libType = " /LIB ";
              }


              x= input$mEMF;

              y = "Fast";
              if(y == "Default"){
                presearch = " Hdva^\n";
              } else if (y == "Fast"){
                presearch = " Hfva^\n";
              } else {
                presearch = " Hsva^\n";
              }


              BATFILE_NAME = "Data\\asm-HS-BAT.bat"
              ProgName = "MSPepSearch\\2017_05_15_MSPepSearch\\x64\\MSPepSearch64.exe";
              Output = "Data\\HitList_HS.txt";
              nHits = 25;  # the number of hits stored (from SS_ESTMW.C)

              sink(BATFILE_NAME)
              cat(paste(ProgName,presearch,sep=""))
              cat(paste(libType,LibDirName,"^\n", sep=""))
              cat(paste(" /HITS ", nHits ,"^\n",sep=""))
              cat(paste(" /INP Data\\Query_HS.MSP ^\n",sep=""));
              cat(" /OutNumMP^\n")
              cat(" /OutMW^\n")
              cat(" /OutDeltaMW^\n")
              cat(" /PROGRESS^\n")
              cat(" /LibInMem^\n")
              cat(paste(" /OUTTAB ", Output,sep=""))
              sink()

              withProgress(message = "completing iHS-HM", value = 0.3, {
                system2(BATFILE_NAME);

              })
              
              Results = readLines(Output);
              Results_w = Results[5:(length(Results)-1)];
              
              withProgress(message = 'obtaining iHS-HM prediction', value=0,{
              
              kl = 1
              
              Omegas = numeric(length(NMEs));
              thetas = numeric(length(NMEs));
              HitsPerSearch = numeric(length(NMEs))
              mEMF = input$mEMF;
                
                for(aj in 1:length(NMEs)){
                  
                  ku = kl+24;
                  
                  ku = min(ku,length(Results_w));
                  if (ku>=length(Results_w)) break
                  
                  PartResult = Results_w[kl:ku]
                  
                  N1 = 0;  # number of potential cognates
                  N2 = 0;  # number of potential isomers
                  
                  ListN1 = c();  # list of hMF - sMF for hits with DeltaMass != 0
                  ListN2 = c();  # list of hmf for hits with DeltaMass == 0
                  
                  for(j in 1:25){
                    
                    name = strsplit(PartResult[j],"\t")[[1]][1]
                    rank = as.numeric(strsplit(PartResult[j],"\t")[[1]][2]);
                    hMF = as.numeric(strsplit(PartResult[j],"\t")[[1]][8]);
                    sMF = as.numeric(strsplit(PartResult[j],"\t")[[1]][13]);
                    
                    if(is.na(rank)) break
                    if(rank != j) break
                    
                    DeltaMass = as.numeric(strsplit(PartResult[j],"\t")[[1]][6])
                    if(hMF >= mEMF){
                      if (DeltaMass == 0){
                        N2 = N2 + 1;
                        diff = ((hMF-mEMF)/(1000-mEMF))*hMF
                        ListN2 = c(ListN2,diff);
                      } else {
                        N1 = N1 + 1;
                        diff = ((hMF-mEMF)/(1000-mEMF))*(hMF-sMF);
                        ListN1 = c(ListN1,diff);
                      }
                    }      
                    
                    
                    kl = kl+1;
                  }
                  
                  
                  if (j>1){
                    HitsPerSearch[aj] = j;
                    
                    if (N1 == 0 && N2 == 0){
                      Omega = 0
                    } else if (N1!=0 && N2 == 0){
                      Omega = 1/N1*sum(ListN1)
                    } else if (N1==0 && N2 == 0){
                      Omega = 1/N2*sum(ListN2)
                    } else {
                      Omega = 1/(N1+N2)*sum(ListN1) + 1/(N1+N2)*sum(ListN2)
                    }
                    
                    Omegas[aj] = Omega
                    thetas[aj] = as.numeric(strsplit(PartResult[j-1],"\t")[[1]][7]) - 
                      as.numeric(strsplit(PartResult[j-1],"\t")[[1]][6])
                    
                  }
                }
                
              
              })
              
                O1 = -sort(-Omegas)[1];
                O2 = -sort(-Omegas)[2];
                
                epsilon = 1e-8;
                maxTheta = thetas[order(-Omegas)][1];
                probTheta = (O1-O2)/999;  
              
            
            values$theta = maxTheta; 
            values$Itheta = probTheta; 
                
            values$AllThetas = thetas;
            values$AllOmegas = Omegas;
            
          } else {
            values$istheta=0;
          }
        
          
          if (4 %in% whichPreds){
            
            a = values$qspectrum_dt;  # local copy of the spectrum as data.table
            #print(a)
            
            dq = values$qspectrum;    # local copy of spectrum as vector
            #cat(dq)
            
            ab_M = input$ab_M;        # minimum abundance of acceptable molecular ion
            
            g = PIM(dq,ab_M);         # compute interpretation based prediction
            
            
            
            
            
            
            
          } else {
            values$isgammabar = 0;
          }
        
     
        return(1)
      
      }
      
    })
    
  
    output$gamma <- renderText({
      a = trigger2()
      
      if(a==0){
        return(" ")
      } else {
        if(values$isgamma==1){
          if(values$Igamma >= input$tau_g){
            message = paste("PIM Prediction: ","<font color=\"#008000\"><b>",values$gamma, "Da","</b></font>");
          } else {
            message = paste("PIM Prediction: ","<font color=\"#FF0000\"><b>",values$gamma, "Da","</b></font>");
          }
          return(message) 
        } else {
          return(" ")
        }
      }
      
    })
    
    output$PIM <- renderText({
      a = trigger2();
      
      if(a==0){
        return(" ")
      } else {
        if(values$isgamma==1){
          
          message1 = "Unexpected Estimated Losses (abundance > 0):"
          if(length(values$elosses)>0){
            message2 = "";
            for(i in 1:length(values$elosses)){
              message2 = paste(message2,values$elosses[i],sep=" ");
            }  
          } else {
            message2 = " ";
          }
          
          messageA = paste(message1," ",message2)
          
          
          
          message1 = "Number of \"significant\" Estimated Losses used in computation:"
          message2 = values$elosses2;
          
          messageB = paste(message1," ",message2)
          
          
          
          message1 = "Question peak abundance to molecular-ion abundance ratio:"
          message2 = round(values$Iratio,2);
          
          messageC = paste(message1," ", message2)
          
          
          
          message1 = "Gamma Classifier Index:"
          message2 = round(values$Igamma,2);
          
          messageD = paste(message1," ", message2)
          
          
          
          message = paste(messageA,"<br><br>", messageB, "<br><br>", messageC, "<br><br>", messageD)
          
          return(message) 
        } else {
          return(" ")
        }
      }
      
    })
    
    
    output$sigma <- renderText({
      a = trigger2()
      
      if (a==0){
        return(" ")
      } else {
        if(values$issigma==1){
          if(values$Isigma >= input$tau_s){
            message = paste("SS-HM Prediction: ","<font color=\"#008000\"><b>",values$sigma, "Da","</b></font>");
          } else {
            message = paste("SS-HM Prediction: ","<font color=\"#FF0000\"><b>",values$sigma, "Da","</b></font>");
          } 
          return(message) 
        } else {
          return(" ")
        }
        
      }
      
    })
    
    output$SSHM <- renderText({
      a = trigger2();
      
      if(a==0){
        return(" ")
      } else {
        if(values$issigma==1){
          
          A = values$SigmaFullResultMatrix;
          
          n = dim(A)[1];
          m = dim(A)[2];
          message = NULL;
          
        
          message = "<table style=\"width:100%\"> <tr> <th> MW </th>  <th> n </th>  <th> I </th>  <th> sMF <th>";
          
          for(i in 5:m){
            message = paste(message,"<th> * </th");
          }
          message = paste(message,"</tr>");

          #message = paste(message, "<br>")
          
          for(i in 1:n){
            message = paste(message, "<tr>")
            for(j in 1:m){
              
                message = paste(message," <td>  ", A[i,j], "  </td>"); 
                
            }
            message = paste(message, "</tr>")
          }
          
          message = paste(message, "</table>")
          
          
          return(message) 
        } else {
          return(" ")
        }
      }
      
    })
    
    
    
    output$theta <- renderText({
      a = trigger2()
      
      if(a==0){
        return(" ")
      } else {
        if(values$istheta==1){
          if(values$Itheta >= input$tau_t){
            message = paste("HS-HM Prediction: ","<font color=\"#008000\"><b>",values$theta, "Da","</b></font>");
          } else {
            message = paste("HS-HM Prediction: ","<font color=\"#FF0000\"><b>",values$theta, "Da","</b></font>");
          }
          return(message)
        } else {
          return(" ")
        }
        
      }
      
    })
    
   
    
    output$OmegaX <- renderPlot({
      
      a = trigger2();
      
      if(a==0){
        return(" ")
      } else {
        if(values$istheta==1){
          
          Omegas = values$AllOmegas
          Thetas = values$AllThetas
          
          Title = paste("Index = ",round(values$Itheta,2))
       
          plot(Thetas,Omegas,pch=19,cex=1.5, ylab = "OmegaX",xlab="Assumed Molecular Mass",main=Title);
          Mstar = Thetas[which.max(Omegas)];
          Mstary = max(Omegas);
          
          abline(v=Mstar,lty=2,col="blue")
          points(Mstar,Mstary,pch=18,cex=2,col="red")
          abline(v=Thetas,lty=3,lwd=0.5,col="grey")
         }
        
      }

    })
    
    session$onSessionEnded(function() {
        stopApp()
    })
    
    
    
}

# Create Shiny app ----
shinyApp(ui, server)
