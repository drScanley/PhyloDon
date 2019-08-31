  #these are a series of functions to plot the variation of Gamma / treelength / Branch lengths / Phylogenetic distance
  #as species are added from the first descriptions to present day
  
  
  #Command list
  #NewTaxNo = number of new taxa per year, or block of years
  #GammaPlot = plot the gamma values of rolled back tree. Also includes a simulated gamma for each time slice using the tree of same OTU number with a random set of taxa removed.
  #LTTTTVideo (Lineages Through Time Through Time) = create a series of sequential images of LTT plots for a range of years
  #TreeGrow =create a series of sequential images of trees for a range of years
  #PDPlot = plot the Phylogenetic distance of taxa described within a specific time-window, PD = total length of distance matrix: PDcorr= PD divided by number of taxa included in the tree: PDcorr2= PD divided by the total pairwise distance of the total tree.
  #TLPlot =plot the Treelength values of rolled back tree. Also includes a simulated treelength for each time slice using the tree of same OTU numebr with a random set of taxa removed.
  #BLplot = plot the average and max min lengths of new Branchlengths between sucessive time slices)
  #GammaPredict = predict the range of gama values of "final tree" run forward from various starting points throughout the taxonomic history 
  #GammaRandomPlot =this adds X random tips to a tree n, where X is the differece in taxa between n and n+1
  #GammaTaxRandomPlot = This function randomizes the taxonomic positions, and runs multiple Gammaplot itterations, comparing that to the itterated randomly removed specimens.
  #GammaModify = Modify branchlength distribution to match a particular value of Gamma
  
  ##################################################
  
  #NewTaxNo = number of new taxa per year, or block of years
  
  ##################################################
  
  #if ntax<nsims then nsims=ntax if not nsims=nsims
  NewTaxNo<-function(phy,data,steps=5)
  {treedata(phy,data,warnings = FALSE)->td
  td$data->data
  as.numeric(data)->data
  max(data,na.rm = TRUE)-min(data,na.rm = TRUE)->yearbreaks
  hist(data,breaks = (yearbreaks/steps),plot = FALSE)->HistOut
  cbind(HistOut$mids+0.5, HistOut$counts)->NewTaxDate
  colnames(NewTaxDate)<-c("year","count")
  NewTaxDate[NewTaxDate == 0] <- NA
  max(HistOut$counts)*1.2->MaxTax
  plot(NewTaxDate,ylim = c(0,MaxTax))
  return(NewTaxDate)}
  
  ##################################################
  
  #GammaPlot = plot the gamma values of rolled back tree. Also includes a simulated gamma for each time slice using the tree of same OTU numebr with a random set of taxa removed. 
  #nsims= number of simulated removals per timeslice, Steps = Number of years per timeslice, Slope window = number of years from present to measure corrected gamma slope, CI Confidence interval range of simulations to be plotted
  
  ##################################################
  
  GammaPlot<-function(phy,data,nsims=10,TITLE=NA,Steps=1,slope_window=30,CI = 95){
    require(phytools)
    require(geiger)
    GTT<-data.frame()
    simGAM<-data.frame()
    #Trim datasets to just include taxa that occur on the tree and the matrix
    treedata(phy,data,warnings = FALSE)->td
    td$phy->phy
    #apparently non-bifurcating trees fail to produce the simulated gamma, so constrain the tree to be bifurcating 
    phy<-multi2di(phy)
    #Identify the  date of the third described species
    sort(td$data)[3]->firstdate
    #Identify the number of taxa in the tree
    length(phy$tip.label)->fullTax
    #Identify the number of years that taxa were described
    nrow(td$data)->fulldate
    #Identify the last date a species was described
    sort(td$data)[fulldate]->lastdate
    #Make both first and last date numbers
    as.numeric(firstdate)->firstdate
    as.numeric(lastdate)->lastdate
    seq(from = firstdate, to = lastdate, by = Steps)->Sequence
    CI/100->CI1
    CI1-1->CI2
    CI2/2->CI3
    1+CI3->CIMax
    0-CI3->CIMin
    #estimate Gamma values for every date
    for (i in seq(from = firstdate, to = lastdate, by = Steps)){
      Sys.sleep(0.1)
      keep<-c(0:i)
      subset(data, data$year %in% keep)->dataR
      treedata(phy,dataR,warnings = FALSE)->td2
      ltt(td2$phy,plot = F)->LTTdata
      LTTdata$gamma->gamma
      GTT<-data.frame(rbind(GTT,gamma))
      length(td2$phy$tip.label)->NTAX
      missingTAX<-fullTax-NTAX
      simgam<-data.frame()
      message(paste("Calculating Gamma for year", i, sep=" "))
      #Estimate the gamma for trees with random taxa removed at the same rate as the original tree
      repeat 
      {
        tr2<-drop.tip(phy,sample(phy$tip.label)[1:missingTAX])
        ltt(tr2, plot = F)->simLTT
        as.numeric(simLTT$gamma)->simGamma
        simgam<-data.frame(rbind(simgam,simGamma))
        if(nrow(simgam)>nsims) break;} 
      colMeans(simgam)->simGamma
      quantile(simgam[,1],CIMin,na.rm = TRUE)->minGamma
      as.numeric(minGamma)->minGamma
      quantile(simgam[,1],CIMax,na.rm = TRUE)->maxGamma
      as.numeric(maxGamma)->maxGamma
      gamma-simGamma->DiffGamma
      as.numeric(DiffGamma)->DiffGamma
      cbind(simGamma,minGamma,maxGamma,DiffGamma)->MinMaxMeanDiff
      simGAM<-rbind(simGAM,MinMaxMeanDiff)}
    #  
    Dates<-(seq(from = firstdate, to = lastdate, by = Steps)+Steps)
    cbind(Dates,GTT,simGAM)->OUTPUT
    colnames(OUTPUT)<-c("Date","Gamma","SimGamma","min","max","diff")
    Min<-min(c(OUTPUT$Gamma,OUTPUT$min,OUTPUT$max))
    Max<-max(c(OUTPUT$Gamma,OUTPUT$min,OUTPUT$max))
    par(mfrow=c(2,1))
    plot(OUTPUT$Date,OUTPUT$SimGamma,xlab = "Year",ylab = "Gamma value", type="l",lty = 2, main = TITLE, ylim=c(Min,Max))
    lines(OUTPUT$Date,OUTPUT$Gamma)
    lines(OUTPUT$Date,OUTPUT$min,lty=3)
    lines(OUTPUT$Date,OUTPUT$max,lty=3)
    legend("topright",legend=c("Gamma plot", "Simulated gamma", "95% confidence interval"), lty=1:3, cex=0.6, bty ="n")
    lastdate-slope_window->prevXyrs
    keep3<-c(prevXyrs:lastdate)
    subset(OUTPUT, OUTPUT$Date %in% keep3)->OutRangeX
    lm(OutRangeX[,2]~OutRangeX[,1])->lmOutRangeX
    lmOutRangeX$coefficients->yrXSlope
    lm(OutRangeX[,3]~OutRangeX[,1])->lmSimOutRangeX
    lmSimOutRangeX$coefficients->SimyrXSlope
    yrXSlope[2]-SimyrXSlope[2]->CorSlope
    max(OUTPUT$Gamma-OUTPUT$min)->GGMax
    min(OUTPUT$Gamma-OUTPUT$max)->GGMin
    
    OutList<-list('output'=OUTPUT,'Gamma slope'=yrXSlope,'Sim slope'=SimyrXSlope,"corrected slope"=CorSlope)
    plot(OUTPUT$Date,(OUTPUT$Gamma-OUTPUT$min),lty=2,xlab = "Year",ylab = "observed-expected Gamma", type="l", ylim=c(GGMin,GGMax))
    lines(OUTPUT$Date,(OUTPUT$Gamma-OUTPUT$max),lty=2)
    abline(h=0,lty=2)
    return(OutList)}
  
  ##################################################
  
  #LTTTTVideo 
  #Lineages Through Time Through Time
  #create a series of sequential images of LTT plots for a range of years
  #Please note that I will try to update this function to produce .avi/.mov/.gif files as soon as I work out how to do so in the animation package!
  
  ##################################################
  
  
  LTTTTVideo<-function(phy,data, Steps=1){
  require(phytools)
  require(geiger)
  GTT<-data.frame()
  simGAM<-data.frame()
  treedata(phy,data,warnings = FALSE)->td
  sort(td$data)[3]->firstdate
  length(phy$tip.label)->fullTax
  sort(td$data)[fullTax]->lastdate
  as.numeric(firstdate)->firstdate
  as.numeric(lastdate)->lastdate
  ltt(phy)->GammOut
  max(log(GammOut$ltt))->GammOut
  
  for (i in seq(from = firstdate, to = lastdate, by = Steps)){
    keep<-c(1757:i)
    message(paste("Recovering tree at year", i, sep=" "))
    subset(data, data$year %in% keep)->dataR
    treedata(phy,dataR,warnings = FALSE)->td2
    filename<-paste("LTT",i,".jpg", sep = "")
    jpeg(file = filename,width = 1500,height = 1200)
    ltt(td2$phy,ylim=c(0,GammOut),main=i)
    dev.off()
  }}
  
  
  ##################################################
  #TreeGrow 
  #create a series of sequential images of trees for a range of years
  
  ##################################################
  
  
  TreeGrow<-function(phy,data,Steps=1){
  require(phytools)
  require(geiger)
  GTT<-data.frame()
  simGAM<-data.frame()
  treedata(phy,data,warnings = FALSE)->td
  sort(td$data)[3]->firstdate
  td$phy->phy
  length(phy$tip.label)->fullTax
  sort(td$data)[fullTax]->lastdate
  as.numeric(firstdate)->firstdate
  as.numeric(lastdate)->lastdate
    
  
  for (i in seq(from = firstdate, to = lastdate, by = Steps)){
    keep<-c(0:i)
    subset(data, data$year %in% keep)->dataR
    treedata(phy,dataR,warnings = FALSE)->td2
    filename<-paste("Tree",i,".jpg", sep = "")
    jpeg(file = filename,width = 1500,height = 1200)
    plot.phylo(td2$phy,type="fan",show.tip.label = FALSE,main=i, cex.main=4)
    dev.off()
  }}
  
  #Phylogenetic distance of taxa described within a specific time-window. Multiple metrics of PD measured
  #PD= sum of pairwise distances (branchlengths) for taxa described in a time window. if less than two taxa were described PD=NA
  #PDcorr= PD divided by number of taxa included in the tree. This acts as a measure of "average pairwise distance per taxon per timeslice"
  #PDcorr2= PD divided by the total pairwise distance of the total tree. This acts as a measure of "proportion of total PD for each timeslice"
  #PDbrln= Total branchlength of subtree
  
  
  ##################################################
  
  #PDPlot 
  #plot the Phylogenetic distance of taxa described within a specific time-window, PD = total length of distance matrix: PDcorr= PD divided by number of taxa included in the tree: PDcorr2= PD divided by the total pairwise distance of the total tree.
  
  ##################################################
  
  
  PDPlot<-function(phy,data,TITLE=NA,Steps=1){
  require(phytools)
  require(geiger)
  PDTT<-data.frame()
  treedata(phy,data,warnings = FALSE)->td
  td$phy->phy
  data.frame(td$data)->data
  length(phy$tip.label)->fullTax
  sort(td$data)[1]->firstdate
  nrow(td$data)->fulldate
  sort(td$data)[fulldate]->lastdate
  as.numeric(firstdate)->firstdate
  as.numeric(lastdate)->lastdate
  cophenetic.phylo(phy)->DistMatrix
  DistMatrix->TotalPD
  length(phy$tip.label)->fullTax
  simulatedPD<-data.frame()

  for (i in seq(from = firstdate, to = lastdate, by = Steps)){
    j<-i+Steps
    keep2<-c(i:j)
    subset(data, data[,1] %in% keep2)->dataR2
    nrow(dataR2)->NrowTest
    PD<-c(0,0,0,0)
    if(NrowTest > 1){
      Sys.sleep(0.1)
      message(paste("Calculating Phylogenetic distance for year", i, sep=" "))
      treedata(phy,dataR2,warnings = FALSE)->td2.2
      td2.2$phy$tip.label->TaxaKeep
      DistMatrix[TaxaKeep,TaxaKeep]->DistMatrixR
      (sum(DistMatrixR)/2)->PD
      PDcorr<- PD/(length(TaxaKeep))
      PDcorr2<-PD/length(TotalPD)
      PDbrln<-sum(td2.2$phy$edge.length)
      PD<-c(PD,PDcorr,PDcorr2,PDbrln)
   
    }
    PDTT<-data.frame(rbind(PDTT,PD))
  }
  
  Dates<-(seq(from = firstdate, to = lastdate, by = Steps)+(Steps/2))
  
  cbind(Dates,PDTT)->OUTPUT
  
  colnames(OUTPUT)<-c("Date","PD","PDcorr","PDcorr2", "PDbrln")
  OUTPUT[OUTPUT == 0] <- NA
  
  
  plot(OUTPUT$Date,OUTPUT$PDbrln,xlab = "Year",ylab = "Phylogenetic Distance",main = TITLE)
  
  
  return(OUTPUT)
  }
  
  ##################################################
  
  #TLPlot 
  #plot the Treelength values of rolled back tree. Also includes a simulated treelength for each time slice using the tree of same OTU numebr with a random set of taxa removed.
  
  ##################################################
  
  
  TLPlot<-function(phy,data,nsims=10,TITLE=NA,Steps = 5){
  require(phytools)
  require(geiger)
  GTT<-data.frame()
  simGAM<-data.frame()
  treedata(phy,data,warnings = FALSE)->td
  sort(td$data)[3]->firstdate
  nrow(td$data)->fulldate
  sort(td$data)[fulldate]->lastdate
  length(phy$tip.label)->fullTax
    as.numeric(firstdate)->firstdate
    as.numeric(lastdate)->lastdate
  
  
    pb <- txtProgressBar(min = firstdate, max = lastdate, style = 3)
    for (i in seq(from = firstdate, to = lastdate, by = Steps)){
      Sys.sleep(0.1)
      message(paste("Calculating treelength for year", i, sep=" "))
      
    keep<-c(0:i)
    subset(data, data$year %in% keep)->dataR
    treedata(phy,dataR,warnings = FALSE)->td2
    sum(td2$phy$edge.length)->gamma
    GTT<-data.frame(rbind(GTT,gamma))
    length(td2$phy$tip.label)->NTAX
    missingTAX<-fullTax-NTAX
    simgam<-data.frame()
    repeat 
    {
      tr2<-drop.tip(phy,sample(phy$tip.label)[1:missingTAX])
      sum(tr2$edge.length)->simGamma
      simgam<-data.frame(rbind(simgam,simGamma))
      if(nrow(simgam)>nsims) break;} 
    colMeans(simgam)->simGamma
    quantile(simgam[,1],0.025)->minGamma
    as.numeric(minGamma)->minGamma
    quantile(simgam[,1],0.975)->maxGamma
    as.numeric(maxGamma)->maxGamma
    simGamma-gamma->DiffGamma
    as.numeric(DiffGamma)->DiffGamma
    cbind(simGamma,minGamma,maxGamma,DiffGamma)->MinMaxMeanDiff
    simGAM<-rbind(simGAM,MinMaxMeanDiff)}
  Dates<-(seq(from = firstdate, to = lastdate, by = Steps)+Steps/2)
  cbind(Dates,GTT,simGAM)->OUTPUT
  colnames(OUTPUT)<-c("Date","Treelength","SimTreelength","min","max","diff")
  Min<-min(c(OUTPUT$Treelength,OUTPUT$min,OUTPUT$max))
  Max<-max(c(OUTPUT$Treelength,OUTPUT$min,OUTPUT$max))
  plot(OUTPUT$Date,OUTPUT$SimTreelength,xlab = "Year",ylab = "Tree length", type="l",lty = 2, main = TITLE, ylim=c(Min,Max))
  lines(OUTPUT$Date,OUTPUT$Treelength)
  lines(OUTPUT$Date,OUTPUT$min,lty=3)
  lines(OUTPUT$Date,OUTPUT$max,lty=3)
  legend("topleft",legend=c("sum tree length", "Simulated tree length", "95% confidence interval"), lty=1:3, cex=0.6, bty ="n")
  return(OUTPUT)
  }
  
  ##################################################
  
  #BLplot = plot the difference between Branchlength n and branchlength n+1 (treelength difference between time slices)
  
  ##################################################
  
  
  BLplot<-function(phy,data,nsims=10,TITLE=NA,graph="Delta",Color="black",Steps=1){
  TLDelta<-data.frame()
  treedata(phy,data,warnings = FALSE)->td
  data.frame(td$data)->data
  sort(td$data)[3]->firstdate
  length(phy$tip.label)->fullTax
  nrow(td$data)->fulldate
  sort(td$data)[fulldate]->lastdate
  as.numeric(firstdate)->firstdate
  as.numeric(lastdate)->lastdate

  
  for (i in seq(from = firstdate, to = lastdate, by = Steps)){
    Sys.sleep(0.1)
    message(paste("Calculating new branchlengths for year", i, sep=" "))
    keep<-c(0:i)
    i+Steps->Steps2
    keep1<-c(0:Steps2)
    subset(data, data[,1] %in% keep)->dataR
    subset(data, data[,1] %in% keep1)->dataR1
    treedata(phy,dataR,warnings = FALSE)->td2
    treedata(phy,dataR1,warnings = FALSE)->td3
    NewTax<-length(td3$phy$tip.label)-length(td2$phy$tip.label)
    sum(td3$phy$edge.length)-sum(td2$phy$edge.length)->NewBr
    NewBr/NewTax->BRav
    TLD<-cbind(0,0,0,0,0,0)
    colnames(TLD)<-c("NewBr","BRav","NewTax","meanEdge","MaxEdge","MinEdge")
    
    setdiff(td3$phy$tip.label,td2$phy$tip.label)->taxlist
    length(taxlist)->taxlength
    if(taxlength!=0){
      nodes<-sapply(taxlist,function(x,y) which(y==x),y=td3$phy$tip.label)
      edge.lengths<-setNames(td3$phy$edge.length[sapply(nodes,function(x,y) which(y==x),y=td3$phy$edge[,2])],names(nodes))
      max(edge.lengths)->MaxEdge
      min(edge.lengths)->MinEdge
      mean(edge.lengths)->meanEdge
      TLD<-data.frame(cbind(NewBr,BRav,NewTax,meanEdge,MaxEdge,MinEdge))
    }
    TLDelta<-data.frame(rbind(TLDelta,TLD))
    length(td2$phy$tip.label)->NTAX
    missingTAX<-fullTax-NTAX}
  Dates<-(seq(from = firstdate, to = lastdate, by = Steps)+Steps)
  cbind(Dates,TLDelta)->OUTPUT
  lastdate-30->prev30yrs
  keep3<-c(prev30yrs:lastdate)
  subset(OUTPUT[,1:2], OUTPUT$Date %in% keep3)->OutRange30
  lm(OutRange30[,2]~OutRange30[,1])->lmOutRange30
  lmOutRange30$coefficients->yr30Slope
  keep4<-c(firstdate:1984)
  subset(OUTPUT[,1:2], OUTPUT$Date %in% keep4)->OutRangefirst
  lm(OutRangefirst)->lmOutRangeFirst
  lmOutRangeFirst$coefficients[2]->FullSlope
  colnames(OUTPUT)<-c("Date","total_new_branch_length","Average_new_branchlength","new_taxa","MeanEdge","MaxEdge","MinEdge")
  OUTPUT[OUTPUT == 0] <- NA 
  subset(OUTPUT, OUTPUT$Average_new_branchlength != "NaN")->OUTPUT
  max(OUTPUT$MaxEdge)->maxValue
  
  plot(OUTPUT$Date,OUTPUT$MeanEdge,xlab="Date", ylab="average new branchlength",col=Color,pch=16,ylim=c(0,maxValue))
  arrows(OUTPUT$Date,OUTPUT$MaxEdge,OUTPUT$Date,OUTPUT$MinEdge,length=0, col=Color, angle=90, code=3)
  points(OUTPUT$Date,OUTPUT$MeanEdge,pch=1)    
  
  OutList<-list('output'=OUTPUT,'slope for last 30 years'=yr30Slope,'1757 to 1984 slope'=FullSlope)
  
  return (OutList) 
  }
  
  
  ##################################################
  
  #GammaPredict 
  #predict the range of gama values of "final tree" run forward from various starting points throughout the taxonomic history 
  
  ##################################################
  
  
  GammaPredict<-function(phy,data,nsims=1000,TITLE=NA,Steps=25,StartDate=1900){
  require(phytools)
  require(geiger)
  GTT<-data.frame()
  simGAM<-data.frame()
  treedata(phy,data,warnings = FALSE)->td
  td$phy->phy
  #apparently non-bifurcating trees fail to produce the simulated gamma, so constrain the tree to be bifurcating 
  phy<-multi2di(phy)
  sort(td$data)[3]->firstdate
  length(phy$tip.label)->fullTax
  nrow(td$data)->fulldate
  sort(td$data)[fulldate]->lastdate
  as.numeric(firstdate)->firstdate
  as.numeric(lastdate)->lastdate
  ltt(td$phy,plot = F)->GamNow
  GamNow$gamma->GamNow
   pb <- txtProgressBar(min = firstdate, max = lastdate, style = 3)
  for (i in seq(from = StartDate, to = lastdate, by = Steps)){
    Sys.sleep(0.1)
    message(paste("simulating present gamma range from year", i, sep=" "))
    keep<-c(1757:i)
    subset(data, data$year %in% keep)->dataR
    treedata(phy,dataR,warnings = FALSE)->td2
    td2$phy->redTree
    length(redTree$tip.label)->redTax
    fullTax-redTax->MissingTax
    ltt(td2$phy,plot = F)->LTTdata
    LTTdata$gamma->gamma
    GTT<-data.frame(rbind(GTT,gamma))
    length(td2$phy$tip.label)->NTAX
    missingTAX<-fullTax-NTAX
    simgam<-data.frame()
    repeat 
    {
      add.random(redTree,n = MissingTax) ->tr2
      ltt(tr2, plot = F)->simLTT
      as.numeric(simLTT$gamma)->simGamma
      simgam<-data.frame(rbind(simgam,simGamma))
      if(nrow(simgam)>nsims) break;} 
    colMeans(simgam)->simGamma
    quantile(simgam[,1],0.025,na.rm = TRUE)->minGamma
    as.numeric(minGamma)->minGamma
    quantile(simgam[,1],0.975,na.rm = TRUE)->maxGamma
    as.numeric(maxGamma)->maxGamma
    gamma-simGamma->DiffGamma
    as.numeric(DiffGamma)->DiffGamma
    cbind(simGamma,minGamma,maxGamma,DiffGamma)->MinMaxMeanDiff
    simGAM<-rbind(simGAM,MinMaxMeanDiff)}
  Dates<-seq(from = StartDate, to = lastdate, by = Steps)
  cbind(Dates,GTT,simGAM)->OUTPUT
  colnames(OUTPUT)<-c("Date","Gamma","SimGamma","min","max","diff")
  Min<-min(c(OUTPUT$Gamma,OUTPUT$min,OUTPUT$max,GamNow))-1
  Max<-max(c(OUTPUT$Gamma,OUTPUT$min,OUTPUT$max,GamNow))+1
  plot(OUTPUT$Date,OUTPUT$min, xlab = "Year", ylab = "Gamma value",pch = 1, ylim=c(Min,Max))
  points(OUTPUT$Date,OUTPUT$max,pch=1)
  abline(h=GamNow, lty=3)
  return(OUTPUT)
  }
  
  
  ##################################################
  
  ##################################################
  
  
  
  #####GammaRandomPlot
    #this adds X random tips to a tree n, where X is the differece in taxa between n and n+1
    GammaRandomPlot<-function(phy,data,nsims=10,TITLE=NA,Steps=1,slope_window=30){
      require(phytools)
      require(geiger)
      GTT<-data.frame()
      simGAM<-data.frame()
      treedata(phy,data,warnings = FALSE)->td
      td$phy->phy
      #apparently non-bifurcating trees fail to produce the simulated gamma, so constrain the tree to be bifurcating 
      phy<-multi2di(phy)
      sort(td$data)[3]->firstdate
      length(phy$tip.label)->fullTax
      nrow(td$data)->fulldate
      sort(td$data)[fulldate]->lastdate
      as.numeric(firstdate)->firstdate
      as.numeric(lastdate)->lastdate
      
      
      for (i in seq(from = firstdate, to = lastdate, by = Steps)){
        keep<-c(0:i)
        i+Steps->Steps2
        keep1<-c(0:Steps2)
        subset(data, data[,1] %in% keep)->dataR
        treedata(phy,dataR,warnings = FALSE)->td2
        ltt(td2$phy,plot = F)->LTTdata
        LTTdata$gamma->gamma 
        
        subset(data, data[,1] %in% keep1)->dataR1
        treedata(phy,dataR1,warnings = FALSE)->td3
        NTips1<-length(td3$phy$tip.label)-length(td2$phy$tip.label)
        Randgam<-data.frame()
  if(NTips1>0){      
        repeat
        {
          add.random(td2$phy,n=NTips1)->RandomTree1
          ltt(RandomTree1,plot=F)->RandomGamma
          RandomGamma$gamma->RandomGamma
          Randgam<-data.frame(rbind(RandomGamma, Randgam))
          if(nrow(Randgam)>nsims) break;} 
        }else{rbind(gamma,gamma)->Randgam}
        
        colMeans(Randgam)->RandGamma
        quantile(Randgam[,1],0.025,na.rm = TRUE)->minRandGamma
        as.numeric(minRandGamma)->minRandGamma
        quantile(Randgam[,1],0.975,na.rm = TRUE)->maxRandGamma
        as.numeric(maxRandGamma)->maxRandGamma
        gamma-RandGamma->DiffRandGamma
        
        cbind(gamma,RandGamma,minRandGamma,maxRandGamma,DiffRandGamma)->MinMaxMeanDiff
        simGAM<-rbind(simGAM,MinMaxMeanDiff)}
      
      Dates<-(seq(from = firstdate, to = lastdate, by = Steps)+Steps)
      cbind(Dates,simGAM)->OUTPUT
      
      colnames(OUTPUT)<-c("Date","Gamma","RandomGamma","MinRandomGamma","MaxRandomGamma","DifferenceRandomGamma")
      max(OUTPUT$MaxRandomGamma,OUTPUT$Gamma)->maxValue
      min(OUTPUT$MinRandomGamma,OUTPUT$Gamma)->minValue
      
      
      plot(OUTPUT$Date,OUTPUT$Gamma,xlab="Date", ylab="average new branchlength",pch=16,ylim=c(minValue,maxValue),type = "l")
      lines(OUTPUT$Date,OUTPUT$MinRandomGamma,xlab="Date",lty=2)
      lines(OUTPUT$Date,OUTPUT$MaxRandomGamma,xlab="Date",lty=2)   
      return (OUTPUT) 
    }
  
  
  ##################################################
  
  ##################################################
  
  
    #GammaTaxRandomPlot
    #This function randomizes the taxonomic positions, and runs multiple Gammaplot itterations, comparing that to the itterated randomly removed specimens.
    GammaTaxRandomPlot<-function(phy,data,nsims=10,TITLE=NA,Steps=1,slope_window=30){
      require(phytools)
      require(geiger)
      GTT<-data.frame()
      simGAM<-data.frame()
      treedata(phy,data,warnings = FALSE)->td
      td$phy->phy
      #apparently non-bifurcating trees fail to produce the simulated gamma, so constrain the tree to be bifurcating 
      phy<-multi2di(phy)
      sort(td$data)[3]->firstdate
      length(phy$tip.label)->fullTax
      nrow(td$data)->fulldate
      sort(td$data)[fulldate]->lastdate
      as.numeric(firstdate)->firstdate
      as.numeric(lastdate)->lastdate
      Dates<-(seq(from = firstdate, to = lastdate, by = Steps)+Steps)
      RandTree<-data.frame(Dates)
      
      
      for (i in seq(from = firstdate, to = lastdate, by = Steps)){
        keep<-c(0:i)
        subset(data, data$year %in% keep)->dataR
        treedata(phy,dataR,warnings = FALSE)->td2
        ltt(td2$phy,plot = F)->LTTdata
        LTTdata$gamma->gamma
        GTT<-data.frame(rbind(GTT,gamma))
        length(td2$phy$tip.label)->NTAX
        missingTAX<-fullTax-NTAX
        simgam<-data.frame()
        repeat 
        {
          tr2<-drop.tip(phy,sample(phy$tip.label)[1:missingTAX])
          ltt(tr2, plot = F)->simLTT
          as.numeric(simLTT$gamma)->simGamma
          simgam<-data.frame(rbind(simgam,simGamma))
          if(nrow(simgam)>nsims) break;} 
        colMeans(simgam)->simGamma
        quantile(simgam[,1],0.025,na.rm = TRUE)->minGamma
        as.numeric(minGamma)->minGamma
        quantile(simgam[,1],0.975,na.rm = TRUE)->maxGamma
        as.numeric(maxGamma)->maxGamma
        gamma-simGamma->DiffGamma
        as.numeric(DiffGamma)->DiffGamma
        cbind(simGamma,minGamma,maxGamma,DiffGamma)->MinMaxMeanDiff
        simGAM<-rbind(simGAM,MinMaxMeanDiff)}
      
      repeat{ 
        sample(data[,1])->RandomizedYears
        cbind(RandomizedYears,data[,1])->RANDOM_YEARS
        row.names(data)->row.names(RANDOM_YEARS)
        RandGTT<-data.frame()
        for (i in seq(from = firstdate, to = lastdate, by = Steps)){
          keep<-c(0:i)
          subset(RANDOM_YEARS, RANDOM_YEARS[,1] %in% keep)->dataR
          treedata(phy,dataR,warnings = FALSE)->td2
          ltt(td2$phy,plot = F)->LTTdata
          LTTdata$gamma->gamma
          RandGTT<-data.frame(rbind(RandGTT,gamma))}
        cbind(data.frame(RandTree,RandGTT))->RandTree
        if(ncol(RandTree)>nsims) break;} 
      
      ncol(RandTree)->c
      nrow(RandTree)->r
      Low5<-data.frame()
      High95<-data.frame()
      
      for (i in 1:r)
      {
        quantile(RandTree[i,3:c],0.025,na.rm = TRUE)->Low
        rbind(Low5,Low)->Low5
      }
      
      for (i in 1:r)
      {
        quantile(RandTree[i,3:c],0.975,na.rm = TRUE)->High
        rbind(High95,High)->High95
      }
      cbind(Dates,GTT,Low5,High95,simGAM)->OUTPUT
      colnames(OUTPUT)<-c("Dates","Gamma","RandomTaxa5","RandomTaxa95","simGamma","minGamma","maxGamma","DiffGamma")
      
      max(OUTPUT$Gamma,OUTPUT$RandomTaxa5,OUTPUT$RandomTaxa95, OUTPUT$maxGamma)->MAX
      min(OUTPUT$Gamma,OUTPUT$RandomTaxa5,OUTPUT$RandomTaxa95,OUTPUT$minGamma)->MIN
      plot(OUTPUT$Dates,OUTPUT$Gamma,ylim = c(MIN,MAX),type="l")
      points(OUTPUT$Dates,OUTPUT$RandomTaxa5)
      points(OUTPUT$Dates,OUTPUT$RandomTaxa95)
      lines(OUTPUT$Dates,OUTPUT$minGamma,lty=2)
      lines(OUTPUT$Dates,OUTPUT$minGamma,lty=2)
      lines(OUTPUT$Dates,OUTPUT$maxGamma,lty=2)
      return(OUTPUT)
    }
    
#################################
#GammaModify
#Modify branchlength distribution to match a particular value of Gamma
#################################    
    GammaModify<-function(tree,r=c(-10,10),g=0)
{ ebTree<-function(tree,r){
  if(r!=0){
    H<-nodeHeights(tree)
    e<-(exp(r*H[,2])-exp(r*H[,1]))/r
    tree$edge.length<-e
  }
  tree
}
EB<-function(tree,r){
    d<-max(nodeHeights(tree))
    tree<-ebTree(tree,r)
    tree$edge.length<-tree$edge.length/max(nodeHeights(tree))*d
    tree
}
gamma<-function(r,tree,g) 
    {(g-ltt(EB(tree,r),plot=FALSE)$gamma)^2}
    
    fit<-optimize(gamma,c(-10,10),tree=tree,g=g,tol=1e-12)
EB(tree,fit$minimum)->ModifiedGTree
return(ModifiedGTree)
    }
    
    #Contmap and barplot showing dates of description
    TreeDates<-function(phy,data){
      treedata(phy,data)->TD
      cbind(as.numeric(TD$data[,1]),as.numeric(TD$data[,1]))->TD_Data
      row.names(TD_Data)<-row.names(TD$data)
      par(mfrow=c(2,1))
      contMap(TD$phy,TD_Data[,1],fsize = Text)
      plotTree.wBars(TD$phy,TD_Data[,1]-min(TD_Data))
      par(mfrow=c(1,1))
      
    }
    
    
###    LambdaPredict = Idenitifies the phylogenetic signal (lambda) of X characters, simulated under Brownian motion, across a tree's taxon
    
    lambdaPredict<-function(phy,data,nsims=1000,TITLE=NA,COL="black", Sigma = 1, Steps=25,StartDate=1758){
      
      BMsims<-fastBM(phy)
      repeat
      {
        fastBM(phy,sig2 = Sigma)->BMout
        BMsims<-cbind(BMsims, BMout)
        if(ncol(BMsims)>nsims) break;} 
      
      Lambda1<-data.frame()
      Lambda2<-data.frame(c(1:nsims))
      LambdaDetails<-data.frame(c("Date","mean","Max","Min"))
      
      
      
      treedata(phy,data,warnings = FALSE)->td
      data.frame(td$data)->data
      sort(td$data)[1]->firstdate
      length(phy$tip.label)->fullTax
      nrow(td$data)->fulldate
      sort(td$data)[fulldate]->lastdate
      as.numeric(firstdate)->firstdate
      as.numeric(lastdate)->lastdate
      
      for (i in seq(from = StartDate, to = lastdate, by = Steps)){
        keep<-c(1757:i)
        subset(data, data$year %in% keep)->dataR
        treedata(phy,dataR,warnings = FALSE)->td2
        treedata(td2$phy,BMsims,warnings = FALSE)->BMsims1
        BMsims1$data->BMsims1
        Sys.sleep(0.1)
        message(paste("Lambda for year", i, sep=" "))
        
        
        Lambda1<-data.frame()
        
        
        
        
        for (j in 1:nsims){
          phylosig(td2$phy,BMsims1[,j],method = "lambda")->PhyloSigOut
          PhyloSigOut$lambda->Lambda
          rbind(Lambda1,Lambda)->Lambda1
          colnames(Lambda1)<-i}
        
        cbind(Lambda2,Lambda1)->Lambda2
        
        quantile(Lambda1,0.025,na.rm = TRUE)->minLambda
        as.numeric(minLambda)->minLambda
        quantile(Lambda1,0.975,na.rm = TRUE)->maxLambda
        as.numeric(maxLambda)->maxLambda
        mean(Lambda1[,1])->AverageLambda
        rbind(i,AverageLambda,maxLambda,minLambda)->LambdaDetails1
        cbind(LambdaDetails,LambdaDetails1)->LambdaDetails
      }
      t(LambdaDetails[,2:ncol(LambdaDetails)])->Details
      plot(Details[,1],Details[,2], ylim = c(min(Details[,4]),max(Details[,3])),col=COL,pch=16, xlab = "Dates",ylab = "Lamda")
      arrows(Details[,1],Details[,3],Details[,1],Details[,4],length=0, angle=90, code=3,col=COL)
      points(Details[,1],Details[,2])
      abline(h=Sigma,lty=2)
      OutList<-list('Lambda'=Lambda2,'Lambda_details'=Details,'BM_sims'=BMsims1)
      return(OutList)  
    }
    
    
##### RatePredict = Idenitifies the evolutioniary rate (mean PIC^2) of X characters, simulated under Brownian motion, across a tree's taxonomic history.
 
    
    
    
    RatePredict<-function(phy,data,nsims=1000,TITLE=NA,Sigma = 1,COL="black", Steps=25,StartDate=1757){
      
      BMsims<-fastBM(phy)
      repeat
      {
        fastBM(phy,sig2 = Sigma)->BMout
        BMsims<-cbind(BMsims, BMout)
        if(ncol(BMsims)>nsims) break;} 
      
      Lambda1<-data.frame()
      Lambda2<-data.frame(c(1:nsims))
      LambdaDetails<-data.frame(c("Date","mean","Max","Min"))
      
      pic(BMout,as.phylo(phy))->ix
      squared_contrasts <- ix^2
      nn <- length(squared_contrasts)
      sigma_sq_pic_out <- mean(squared_contrasts)
      
      treedata(phy,data,warnings = FALSE)->td
      data.frame(td$data)->data
      sort(td$data)[1]->firstdate
      length(phy$tip.label)->fullTax
      nrow(td$data)->fulldate
      sort(td$data)[fulldate]->lastdate
      as.numeric(firstdate)->firstdate
      as.numeric(lastdate)->lastdate
      Lambda1<-data.frame()
      
      
      for (i in seq(from = StartDate, to = lastdate, by = Steps)){
        keep<-c(1750:i)
        subset(data, data$year %in% keep)->dataR
        treedata(phy,dataR,warnings = FALSE)->td2
        td2$phy->phy2
        treedata(td2$phy,BMsims,warnings = FALSE)->BMsims1
        BMsims1$data->BMsims1
        Sys.sleep(0.1)
        message(paste("Rates for year", i, sep=" "))
        
        
        Lambda1<-data.frame()
        
        
        
        
        for (j in 1:nsims){
          pic(BMsims1[,j],as.phylo(phy2))->ix
          squared_contrasts <- ix^2
          nn <- length(squared_contrasts)
          sigma_sq_pic <- mean(squared_contrasts)
          rbind(Lambda1,sigma_sq_pic)->Lambda1
          colnames(Lambda1)<-i}
        cbind(Lambda2,Lambda1)->Lambda2
        quantile(Lambda1,0.025,na.rm = TRUE)->minLambda
        as.numeric(minLambda)->minLambda
        quantile(Lambda1,0.975,na.rm = TRUE)->maxLambda
        as.numeric(maxLambda)->maxLambda
        mean(Lambda1[,1])->AverageLambda
        rbind(i,AverageLambda,maxLambda,minLambda)->LambdaDetails1
        cbind(LambdaDetails,LambdaDetails1)->LambdaDetails
      }
      t(LambdaDetails[,2:ncol(LambdaDetails)])->Details
      plot(Details[,1],Details[,2], ylim = c(min(Details[,4]),max(Details[,3])),col=COL,pch=16, xlab = "Dates",ylab = "rate")
      arrows(Details[,1],Details[,3],Details[,1],Details[,4],length=0, angle=90, code=3,col=COL)
      points(Details[,1],Details[,2])
      abline(h=sigma_sq_pic_out,lty=2)
      OutList<-list('Lambda'=Lambda2,'Details'=Details,'BM_sims'=BMsims1)
      
      return(OutList)  
    }
    
    ###############################################
    ###############################################
    ###############################################
    
    ###    LambdaPredict = Idenitifies the phylogenetic signal (lambda) of a continouous character, across a tree's taxon
    
    SVLLambda<-function(phy,data,SVL,TITLE=NA,COL="black", Sigma = 1, Steps=25,StartDate=1758){
      treedata(phy,SVL)->TD1
      TD1$data->BMout
      TD1$phy->phy
      BMout[,1]->BMsims
      
      
      
      Lambda1<-data.frame()
      Lambda2<-data.frame(1)
      LambdaDetails<-data.frame(c("Date","mean","Max","Min"))
      
      
      
      treedata(phy,data,warnings = FALSE)->td
      data.frame(td$data)->data
      sort(td$data)[1]->firstdate
      length(phy$tip.label)->fullTax
      nrow(td$data)->fulldate
      sort(td$data)[fulldate]->lastdate
      as.numeric(firstdate)->firstdate
      as.numeric(lastdate)->lastdate
      
      phylosig(td$phy,as.matrix(BMsims),method = "lambda")->PhyloSigOut1
      PhyloSigOut1$lambda->LambdaInit
      
      
      for (i in seq(from = StartDate, to = lastdate, by = Steps)){
        keep<-c(1757:i)
        subset(data, data$year %in% keep)->dataR
        treedata(phy,dataR,warnings = FALSE)->td2
        treedata(td2$phy,BMsims,warnings = FALSE)->BMsims1
        BMsims1$data->BMsims1
        Sys.sleep(0.1)
        message(paste("Lambda for year", i, sep=" "))
        
        
        Lambda1<-data.frame()
        
        
        
        
        phylosig(td2$phy,BMsims1,method = "lambda")->PhyloSigOut
        PhyloSigOut$lambda->Lambda
        rbind(Lambda1,Lambda)->Lambda1
        colnames(Lambda1)<-i
        
        cbind(Lambda2,Lambda1)->Lambda2
        
        quantile(Lambda1,0.025,na.rm = TRUE)->minLambda
        as.numeric(minLambda)->minLambda
        quantile(Lambda1,0.975,na.rm = TRUE)->maxLambda
        as.numeric(maxLambda)->maxLambda
        mean(Lambda1[,1])->AverageLambda
        rbind(i,AverageLambda,maxLambda,minLambda)->LambdaDetails1
        cbind(LambdaDetails,LambdaDetails1)->LambdaDetails
      }
      t(LambdaDetails[,2:ncol(LambdaDetails)])->Details
      wtf <- lm(Details[,1] ~ Details[,2])
      summary(wtf)
      plot(Details[,1],Details[,2], ylim = c(min(Details[,4]),max(Details[,3])),col=COL,pch=16, xlab = "Dates",ylab = "Lamda")
      arrows(Details[,1],Details[,3],Details[,1],Details[,4],length=0, angle=90, code=3,col=COL)
      points(Details[,1],Details[,2])
      abline(h=LambdaInit ,lty=2)
      OutList<-list('Lambda'=Lambda2,'Lambda_details'=Details,'BM_sims'=BMsims1)
      return(OutList)  
    }
    
    ###############################################
    ###############################################
    ###############################################
    
    
    
    
    
    #SVLrates = Idenitifies the evolutioniary rate (mean PIC^2) of a continuous character, across a tree's taxonomic history.
    
    SVLrates<-function(phy,data,SVL,TITLE=NA,COL="black", Sigma = 1, Steps=25,StartDate=1758){
      treedata(phy,SVL,warnings = FALSE)->TD1
      TD1$data->BMout
      TD1$phy->phy
      BMout[,1]->BMsims
      
      
      pic(as.matrix(BMsims),as.phylo(phy))->ix
      squared_contrasts1 <- ix^2
      nn1 <- length(squared_contrasts1)
      sigma_sq_pic1 <- mean(squared_contrasts1)
      Rates1<-data.frame()
      Rates2<-data.frame(1)
      RatesDetails<-data.frame(c("Date","mean","Max","Min"))
      
      
      
      treedata(phy,data,warnings = FALSE)->td
      data.frame(td$data)->data
      sort(td$data)[1]->firstdate
      length(phy$tip.label)->fullTax
      nrow(td$data)->fulldate
      sort(td$data)[fulldate]->lastdate
      as.numeric(firstdate)->firstdate
      as.numeric(lastdate)->lastdate
      
      
      for (i in seq(from = StartDate, to = lastdate, by = Steps)){
        keep<-c(1757:i)
        subset(data, data$year %in% keep)->dataR
        treedata(phy,dataR,warnings = FALSE)->td2
        td2$phy->phy2
        treedata(td2$phy,BMsims,warnings = FALSE)->BMsims1
        BMsims1$data->BMsims1
        Sys.sleep(0.1)
        message(paste("Lambda for year", i, sep=" "))
        
        
        Rates1<-data.frame()
        
        
        
        pic(BMsims1,as.phylo(phy2))->ix
        squared_contrasts <- ix^2
        nn <- length(squared_contrasts)
        sigma_sq_pic <- mean(squared_contrasts)
        rbind(Rates1,sigma_sq_pic)->Rates1
        colnames(Rates1)<-i
        
        cbind(Rates2,Rates1)->Rates2
        
        quantile(Rates1,0.025,na.rm = TRUE)->minRates
        as.numeric(minRates)->minRates
        quantile(Rates1,0.975,na.rm = TRUE)->maxRates
        as.numeric(maxRates)->maxRates
        mean(Rates1[,1])->AverageRates
        rbind(i,AverageRates,maxRates,minRates)->RatesDetails1
        cbind(RatesDetails,RatesDetails1)->RatesDetails
      }
      t(RatesDetails[,2:ncol(RatesDetails)])->Details
      Details[,1:2]->Details
      plot(Details[,1],Details[,2], ylim = c(min(Details[,2]),max(Details[,2])),col=COL,pch=16, xlab = "Dates",ylab = "rate")
      points(Details[,1],Details[,2])
      abline(h=sigma_sq_pic1 ,lty=2)
      wtf <- lm(Details[,1] ~ Details[,2])
      summary(wtf)
      OutList<-list('Lambda'=Rates2,'RatesDetails'=Details,'BM_sims'=BMsims1,'Summary'=summary(wtf))
      return(OutList)
      return()
    }
    