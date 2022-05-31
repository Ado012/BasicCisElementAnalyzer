#crmfinder script v0.5 Albert Do 07/18/17 



#Sample Command

#arabMotifSearch(targetfile = "Rdatafiles/upregulatedwusgenes.csv", resultsOutput="Rdatafiles/coreMotifOutputUP7.txt", resultsBed="Rdatafiles/coreMotifOutputBEDUP7.bed", resultsSorted="Rdatafiles/coreMotifOutputSortedUP7.txt")


#INPUTS
#targetfile: List of genes you're going to be looking through. Derived from previous inhouse data.
#annotation file : currently hard coded to TAIR10_GFF3_genes_transposons.csv from TAIR?
#sequence files : hardcoded to chromosome files (eg Arabidopsis_thaliana.TAIR10.dna.chromosome.MT.fa) from TAIR/ENSEMBL/PlantGDB look in script for more info. 

#OUTPUTS
#resultsOutput: results
#resultsBed: results in BED format
#resultsSorted: sorted results. 

#The script currently expects a very specific set of inputs and might not work properly if for example the chromosome names are listed differently. Look in script for proper
#names. 



library(Biostrings)
library(combinat)
library(stringr)

#letter lowering function
lowerMotifLetters<-function(stringVector,numVector,hitOffset)
{
  
  
  for (i in 1:length(numVector))
  {
    
    
    stringVector[numVector[i]:(numVector[i]+hitOffset)]<-tolower(stringVector[numVector[i]:(numVector[i]+hitOffset)])
  }
  
  processedStringVector<-stringVector
}





#chromosome picking function
chromPicker <- function(targetGeneChromo)
{
  #Load Proper Chromosome File
  
  if (targetGeneChromo == 'ChrM')
  {
    retrievedSeq <-
      readDNAStringSet(filepath = "Rdatafiles/Arabidopsis_thaliana.TAIR10.dna.chromosome.MT.fa", format =
                         "fasta")
    targetGeneChromoName <- 'ChrM'
  }
  else if (targetGeneChromo == 'Chr1')
  {
    retrievedSeq <-
      readDNAStringSet(filepath = "Rdatafiles/Arabidopsis_thaliana.TAIR10.dna.chromosome.1.fa", format =
                         "fasta")
    targetGeneChromoName <- 'Chr1'
  }
  else if (targetGeneChromo == 'Chr2')
  {
    retrievedSeq <-
      readDNAStringSet(filepath = "Rdatafiles/Arabidopsis_thaliana.TAIR10.dna.chromosome.2.fa", format =
                         "fasta")
    targetGeneChromoName <- 'Chr2'
  }
  else if (targetGeneChromo == 'Chr3')
  {
    retrievedSeq <-
      readDNAStringSet(filepath = "Rdatafiles/Arabidopsis_thaliana.TAIR10.dna.chromosome.3.fa", format =
                         "fasta")
    targetGeneChromoName <- 'Chr3'
  }
  else if (targetGeneChromo == 'Chr4')
  {
    retrievedSeq <-
      readDNAStringSet(filepath = "Rdatafiles/Arabidopsis_thaliana.TAIR10.dna.chromosome.4.fa", format =
                         "fasta")
    targetGeneChromoName <- 'Chr4'
  }
  else if (targetGeneChromo == 'Chr5')
  {
    retrievedSeq <-
      readDNAStringSet(filepath = "Rdatafiles/Arabidopsis_thaliana.TAIR10.dna.chromosome.5.fa", format =
                         "fasta")
    targetGeneChromoName <- 'Chr5'
  }
  else if (targetGeneChromo == 'ChrC')
  {
    retrievedSeq <-
      readDNAStringSet(filepath = "Rdatafiles/Arabidopsis_thaliana.TAIR10.dna.chromosome.Pt.fa", format =
                         "fasta")
    targetGeneChromoName <- 'ChrC'
  }
  else
    print("NO recognized chromosome")
  

  #may fail if no recognized chromasome.
  retrievedSeqInfo <- list("sequence"= retrievedSeq, "chrname"=targetGeneChromoName)
  
  
  
}



#extracts info from gene info file
targetScanner <- function(target, geneNamesDataTS, geneChromListTS)
{
  #set range to search around gene
  searchRange = 3000
  

  
  geneNamesDataUL<- unlist(geneNamesDataTS[1])
  #print(head(unlist(geneNamesData)))
  #print(geneNamesData[1])
  
  #print(geneNamesData[1])
  
  #get the names of the gene transcripts.
  geneNamesList <- lapply(geneNamesDataUL, function(geneEntry)
  {
    geneEntryChar<-as.character(geneEntry)
    geneNameDataSplit <- strsplit(geneEntryChar, "[[:punct:][:space:]]+")[[1]]
    
    #print("Gene Data")
    #print(geneEntryChar)
    
    geneName <- geneNameDataSplit[2]
    
    
  })
  
  #find target gene in Arabidopsis gene list
  targetLocs <- grep(target, geneNamesList)

  
  #take the first hit, possibly modify later
  targetLoc <- targetLocs[1]
  
  
  #load proper chromasome file and find chrom location.
  arabSeqInfo <- chromPicker(targetGeneChromo=as.character(geneChromListTS[targetLoc]))

  
  #extract sequence and chromosome.
  arabSeq <- arabSeqInfo$sequence
  chrom <- arabSeqInfo$chrname
  
  
  targetLocInfoTS <- list("targetLoc"= targetLoc, "sequence"= arabSeq, "chrname"= chrom)
  
    
}





#finds motif hits across target flanking regions
motifPrelimfinder <- function(targetLocMPF, geneStartListMPF, geneEndListMPF,searchRangeMPF, arabSeqMPF)
{
  
 
  print("target list input check")
  print(head(targetLocMPF))
  
  print("check4")
  #load gene end and beginning and search range
  searchRangeStart <- geneStartListMPF[targetLocMPF] - searchRangeMPF
  searchRangeEnd  <- geneEndListMPF[targetLocMPF] + searchRangeMPF
  
  print("check5")
  
  #If Search range extends outside chromosome boundaries. Rein them in
  if (searchRangeStart < 1)
  {
    searchRangeStart = 1
    startadjustflag = 1 #flag might not trip if already at 1?
  }
  if (searchRangeEnd > width(arabSeqMPF[1]))
  {
    searchRangeEnd = width(arabSeqMPF[1])
    endadjustflag = 1
  }
  # Retrieve sequence
  arabSeqUnlist <- unlist(arabSeqMPF[1])

  print("check6")
  
  #retrieve sequence in search range
  #(-1,+1) don't include gene boundaries in search
  targetArabSeqLeft <-
    arabSeqUnlist[searchRangeStart:(geneStartListMPF[targetLocMPF]-1)]
  targetArabSeqRight <-
    arabSeqUnlist[(geneEndListMPF[targetLocMPF]+1):searchRangeEnd]

  print("Coordinate: 5' search start")
  print(searchRangeStart)
  print("Coordinate: 3' search start")
  print(geneStartListMPF[targetLocMPF]-1)
  print("seqwidth")
  print(width(arabSeqMPF[1]))
  
  #Search for cores to the left and right of gene
  motifRangesLeft <- str_locate_all(targetArabSeqLeft,"(?=(TAAT|ATTA))")  
  motifHitsLeft<-motifRangesLeft[[1]][,1]
  
  print("check7")
  
  
  motifRangesRight <- str_locate_all(targetArabSeqRight,"(?=(TAAT|ATTA))")  
  motifHitsRight<-motifRangesRight[[1]][,1]
  
  #adjust hits back to absolute sequence values
  motifHitsLeftAdj <-
    sapply(motifHitsLeft, function(motif)
      motif + (searchRangeStart-1))  #subtract one because adding in search range start adds an extra base
  motifHitsRightAdj <-
    sapply(motifHitsRight, function(motif)
      motif + (geneEndListMPF[targetLocMPF]))
  
  print("check8")
  
  #merge hits into one list and convert to vector
  motifHits <- c(motifHitsLeftAdj, motifHitsRightAdj)
  
  print("check9")
  
  motifHitPositions <- unlist(motifHits)
  
  
}






#determines clusters present in sequence from motif hit results 
clusterScanner <- function(motifHits, arabSeqCS)
{
  
  #initialize flags and containers to mark clusters and their boundaries
  clusterFlag = 0
  clusterStartList <- c()
  clusterEndList <- c()
  clusterSeqList <- c()
  clusterStart <- c()
  
  print("check10")
  
  #loop to determine clusters
  for (i in 1:length(motifHits))
  {
    print("entered cluster loop")
    
    if ((i + 1) <= length(motifHits)) #if still in motif hits list
    {
      print("loop is not finished, still in cluster?")
      if ((abs(motifHits[i] - motifHits[i + 1]) <= 50)) #if motif hit is 50 or less from the last moftif hit
      {
        if (clusterFlag == 0) #if potential cluster is not yet started start it
        {
          clusterStart <- motifHits[i]
          
          print("Begin cluster")
          print(motifHits[i])
        }
        
        clusterFlag = clusterFlag + 1 #extend potential cluster
        print("Extend Cluster")
        print(motifHits[i])
        
      }
      
      
      else if ((abs(motifHits[i] - motifHits[i + 1]) > 50)) #if motif is more than 50 away from last then end potential cluster
      {
        print("potential cluster ends, but loop is not finished")
        
        
        if (clusterFlag >= 4) #if potential cluster contains 4 or more motifs add to real cluster list 
        {
          clusterStartList <- append(clusterStartList,
                                     clusterStart) #needs to be adjusted like this for some reason
          clusterEndList <- append(clusterEndList,
                                   motifHits[i] + 3)
          
          print("cluster start")
          print(tail(clusterStartList, n = 1))
          print("cluster end")
          print((motifHits[i] + 3))
          
          #add cluster sequence to list
          clusterSeq <-
            DNAStringSet(arabSeqCS,
                         tail(clusterStartList, n = 1),
                         tail(clusterEndList, n = 1))
          
          print("clusterSeq")
          print(clusterSeq)
          
          
          clusterSequlchar <-
            as.character(unlist(clusterSeq))
          
          print("clusterSeq to Append")
          print(clusterSequlchar)
          
          print("cluster end detected, before loop end") 
          clusterSeqList <- append(clusterSeqList, clusterSequlchar)
        }
        
        clusterFlag = 0
        
      }
      
    }
    
    else if ((i + 1 > length(motifHits))) #if at the end of motif list. Just end cluster search and determine if you have a cluster with what you already have
    {
      print("motif loop at end")
     
       print("clusterFlag")
      print(clusterFlag)
      
      
      if (clusterFlag >= 4)
      {
        clusterStartList <- append(clusterStartList, clusterStart)
        clusterEndList <- append(clusterEndList,
                                 motifHits[i] + 3) #add 3 due to TAAT core length
        
        
        print("cluster start")
        print(tail(clusterStartList, n = 1))
        print("cluster end")
        print(tail(clusterEndList, n = 1))
        
        clusterSeq <-
          DNAStringSet(arabSeqCS,
                       tail(clusterStartList, n = 1),
                       tail(clusterEndList, n = 1))
        
        clusterSequlchar <-
          as.character(unlist(clusterSeq))
        
        print("clusterSeq to Append")
        print(clusterSequlchar)
        
        print("cluster end detected at loop end")
        clusterSeqList <- append(clusterSeqList, clusterSequlchar)
      }
      
      
      clusterFlag = 0
    }
    
    
  }
  
  print("listcheck")
  print(clusterStartList)
  print(clusterEndList)
  print(clusterSeqList)
  
  #Bundle together cluster coordinates and sequence for return
  clusterRange <-
    list(
      "clusterStartList" = clusterStartList,
      "clusterEndList" = clusterEndList,
      "clusterSeqList" = clusterSeqList
    )
  
  print("# of clusters in region")
  print(length(clusterRange$clusterStartList))
  
  
  clusterRange
}


#Determines complex cores present in sequence. Clusters in this case. 
complexCoreScanner <-
  function(clusterStartCCS,
           clusterEndCCS,
           clusterSeqCCS)
  {
    
    #loop through each cluster finding complex cores
    for (i in 1:length(clusterStartCCS))
    {
      print("CCS1")
      
      #search for cores in cluster
      clusterCoreRanges <-
        str_locate_all(clusterSeqCCS[i], "(?=(TAAT|ATTA))")
      clusterCoreHits <- clusterCoreRanges[[1]][, 1]
      
      #prepare flags and containers to mark complex cores
      clusternum = i
      compCoreFlag = 0
      compCoreStartList <- c()
      compCoreEndList <- c()
      clusterInfoList <- c()
      
      print("cluster core hits")
      print(length(clusterCoreHits))
      
      
      #complex core marking loop
      for (i in 1:length(clusterCoreHits))
      {
        print("i")
        print(i)
        
        #if loop is not over
        if (length(clusterCoreHits) >= (i + 1))
        {
          if (abs(clusterCoreHits[i] - clusterCoreHits[i + 1]) <= 4)
          {#initialize or extend complex core if cores are close enough
            if (compCoreFlag==0)
            { compCoreStart<-clusterCoreHits[i]
            print("comp core start")
            print(clusterCoreHits[i])
            }
            
            compCoreFlag = compCoreFlag+1
            print("Complex core length increased")
            print(clusterCoreHits[i])
          }
          
          #if core spacing is high and a complex core is already started
          else if ((abs(clusterCoreHits[i] - clusterCoreHits[i + 1]) > 4) &&
                   compCoreFlag >= 1)
          {#end the core and upload info 
            compCoreStartList <- append(compCoreStartList, compCoreStart)
            compCoreEndList <- append(compCoreEndList, clusterCoreHits[i]+3)
            compCoreFlag = 0
            print("Complex core end detected")
            print(clusterCoreHits[i+1])
            
          }
          
        }
        
        else if (compCoreFlag>=1)
        {#complex core end is detected and uploaded at end of loop
          compCoreStartList <- append(compCoreStartList, compCoreStart)

          compCoreEndList <-
          append(compCoreEndList, clusterCoreHits[i] + 3)


          compCoreFlag=0
          print("Complex core end detected with end of loop")
          print(clusterCoreHits[i])
          print(clusterCoreHits[i+1]+3)
          print("length of loop")
          print(length(clusterCoreHits))
        }
        
      }
      
      print("clusterSeq")
      print(clusterSeqCCS)
      print("cluster start and end")
      print(clusterStartCCS)
      print(clusterEndCCS)
     
      #process sequence to motif lowering function to improve readability
      ntStringVec<-unlist(strsplit(clusterSeqCCS, split=''))
      scanHits <- str_locate_all(clusterSeqCCS,"(?=(TAAT|ATTA))")
      stringHits<-unlist(scanHits[[1]][,1]) #start hits must be seperated from end hits
      ntStringVec<-lowerMotifLetters(numVec=stringHits,stringVec=ntStringVec,hitOffset=3)
      clusterSeqCCS<-paste(ntStringVec,collapse='')
      
      
      #Bundle complex core along with cluster information for return
      clusterInfo <-
        list(
          "clusterSeq" = clusterSeqCCS[clusternum] ,
          "clusterStart" = clusterStartCCS[clusternum] ,
          "clusterEnd" = clusterEndCCS[clusternum],
          "coreHits" = clusterCoreHits ,
          "complexCores" = length(compCoreStartList),
          "compCoreStartList" = compCoreStartList,
          "compCoreEndList" = compCoreEndList
        )
      
      
      clusterInfoList <- append(clusterInfoList, clusterInfo)
    }
    print("CCS4")
    
    clusterInfoList
  }





#Scores the clusters based upon the associated info
motifScorer <- function( clusterInfoMS  )
{
  motifPhaseScoreSum=0
  #calculate phasing score from the core positions
coreHits<-clusterInfoMS$coreHits
      phaseScore=0
      
      
      
      #repeat slowly
      #look around each core in cluster from left and right margins to see if surrounding cores are in a phasing pattern
      for (i in 1:length(coreHits))
        
      {
        j = 1
        #while distance is less than 50, test phasing among front cores
        while ((i+j <= length(coreHits)) && abs(coreHits[i] - coreHits[i + j]) < 50)
        {
          frontPhase <- abs(coreHits[i] - coreHits[i + j])
          
          
          if (frontPhase >= 9 &&
              frontPhase <= 11 ||
              frontPhase >= 20 &&
              frontPhase <= 22 ||
              frontPhase >= 30 &&
              frontPhase <= 33 || frontPhase >= 41 && frontPhase <= 44)
          {phaseScore = phaseScore + 1}
          
          
          else if (frontPhase == 8 ||
                   frontPhase == 12 ||
                   frontPhase == 19 ||
                   frontPhase == 23 ||
                   frontPhase == 29 ||
                   frontPhase == 34 || frontPhase == 40 || frontPhase == 45)
          {phaseScore = phaseScore + 0.5}
          
          j = j + 1
        }
        
        
      }
      
      #sum scores in score for entire cluster
      motifPhaseScoreSum<-phaseScore
      
      print("motifPhaseScoreSum")
      print(motifPhaseScoreSum)
      print("compCoreEndList")
      print(clusterInfoMS$compCoreEndList)
      print("compCoreStartList")
      print(clusterInfoMS$compCoreStartList)
      
      
      #determine complex core score by subtracting complex core ends and starts and summing
      motifComplexCoreScore<-mapply(function(compCoreStart,compCoreEnd){
        if (exists("compCoreStart") && exists("compCoreEnd")) #check later if this is working
        compCoreEnd-compCoreStart},clusterInfoMS$compCoreStartList,clusterInfoMS$compCoreEndList)
      
      print("motifComplexCoreScores")
      print(motifComplexCoreScore)
      
      motifComplexCoreScoreSum<-sum(unlist(motifComplexCoreScore))
      
      print("motifComplexCoreScoreSum")
      print(motifComplexCoreScoreSum)
      
      clusterSeqULchar<-as.character(unlist(clusterInfoMS$clusterSeq))
      
      phasePerBase<-motifPhaseScoreSum/nchar(clusterSeqULchar)
      phasePerBase<-round(phasePerBase,5)
      
     #Bundle scores along with other information for return
      clusterInfo<-list("clusterSeq"=clusterSeqULchar, "clusterStart"=clusterInfoMS$clusterStart ,"clusterEnd"=clusterInfoMS$clusterEnd, 
                        "coreHits"=clusterInfoMS$coreHits , "complexCores"=clusterInfoMS$complexCores, 
                        "compCoreScore"=motifComplexCoreScoreSum, "phaseScore"=motifPhaseScoreSum, "phasePerBase"=phasePerBase) 
      
  
  
}


#writes results to file
motifWriter<-function(clusterInfoMW, chromMW, targetMW, arabSeqMW, fileStreamMW, fileStreamBAMMW,targetDataMW)
{
  
  print("checkMW")
  
  #+3 adjust end to account for motif size
  #Bundle info to write to file
  motifCluster <-
    list(
      "chromosome"=chromMW,
      "gene"=targetMW,
      "chainstart"=clusterInfoMW$clusterStart,
      "chainend"=clusterInfoMW$clusterEnd,
      "phasescore"=clusterInfoMW$phaseScore,
      "phasePerBase"=clusterInfoMW$phasePerBase,
      "coreNum"=length(clusterInfoMW$coreHits),
      "complexCores"=clusterInfoMW$complexCore,
      "complexCoreScore"=clusterInfoMW$compCoreScore,
      "sequence"=clusterInfoMW$clusterSeq,
      "metadata"=targetDataMW
      
    )
  
  print("motifCluster")
  print(motifCluster)
  
#arrange info into data frame  
  motifClusterdf<-as.data.frame(do.call(cbind,motifCluster))
  
  print("motifClusterdf")
  print(motifClusterdf)
  
  #write info to file
  write.table(motifClusterdf,append=TRUE, col.names=FALSE, row.names=FALSE, quote=FALSE, sep=' ', fileStreamMW)
  print("babe")
  
  #write info to bed file 
  write.table(data.frame(chromMW,clusterInfoMW$clusterStart,clusterInfoMW$clusterEnd),append=TRUE, col.names= FALSE, row.names=FALSE, quote=FALSE, sep="\t",fileStreamBAMMW)
  
  #return to main function for further processing
  motifClusterdf
}


#main function
arabMotifSearch <- function(targetfile, resultsOutput="coreMotifOutput", resultsBed="coreMotifOutputBED", resultsSorted="coreMotifOutputSorted")
{
  searchRange = 3000

  #read in annotation file: Main source of info on the genes we'll be scanning
  geneAnnotations <-
    read.csv(
      "Rdatafiles/TAIR10_GFF3_genes_transposons.csv",
      header = TRUE,
      sep = ','
    )
  
  
  #Break up gene annotation data into relevant variables.
  geneChromList <- unlist(geneAnnotations[1]) #chromosome of gene feature
  geneStartList <- unlist(geneAnnotations[4]) #start position of gene features
  geneEndList <- unlist(geneAnnotations[5]) #end position of gene feature
  geneNamesData <- geneAnnotations[9] #a column listingt gene names the feature is associated with
  
  #start file streams that will be used for output
  fileStream <- file(resultsOutput, 'w')
  fileStreamBAM <- file(resultsBed, 'w')
  
  #start dataframe to store results
  resultHeader<-c("chromosome gene chainstart chainend phasescore phasePerBase corNum complexCores complexCoreScore sequence metadata")
  write.table(resultHeader,append=TRUE, col.names=FALSE, row.names=FALSE, quote=FALSE, fileStream)
  
  #read in file of gene targets that will be examined
  targetdata <-
    read.csv(targetfile,
             header = TRUE,
             sep = ',')
  
  targets <- targetdata[2]
  targetsul <- unlist(targets)
  targetEntry <- as.character(targetsul)
  
  targetMeta<- unlist(targetdata[3])
  targetMetachar<-as.character(targetMeta)
  
  
  
  #extract targets one by one
  for(j in 1:length(targetEntry))
  {
    
    print("main loop check")
    print(j)
    print(length(targetEntry))
    
    #TODO See whetehr activated and repressedtargets are properly passed.
    #grabs target and associates it with the sequence and chromosome from the annotation data
    targetLocInfo <-
      targetScanner(
        target = targetEntry[j],
        geneNamesDataTS = geneNamesData,
        geneChromListTS = geneChromList
      )
    
    targetLoc <- targetLocInfo$targetLoc
    arabSeq <- targetLocInfo$sequence
    chrom <- targetLocInfo$chrname
    
    
    
    print("check3")
    
    #scan target regions for TAAT core hits
    motifHitPositions <-
      motifPrelimfinder(
        targetLocMPF = targetLoc,
        geneStartListMPF = geneStartList,
        geneEndListMPF = geneEndList,
        searchRangeMPF = searchRange,
        arabSeqMPF = arabSeq
      )
    
    #detect clusters in from TAAT core hits lists and the sequence
    clusterList<-clusterScanner(motifHits=motifHitPositions, arabSeqCS=arabSeq)
    
    print("cluster List")
    print(clusterList)
    print(length(clusterList))
    print("cluster seq list")
    print(length(clusterList$clusterSeqList))
    
    #if clusters are present
    if (length(clusterList$clusterSeqList>=1))#for loop may be executed even with zero length lists
    {
    #for each cluster
    for(i in 1:length(clusterList$clusterSeqList))#using elements because the list itself has an extra empty for some reason
    {
      

      print("check13")

      #scan for complex cores
      clusterInfo<-complexCoreScanner(clusterStartCCS=clusterList$clusterStartList[i], clusterEndCCS=clusterList$clusterEndList[i], 
                                      clusterSeqCCS=clusterList$clusterSeq[i])
        
      print("check14")
      
      #score 
      clusterAnnotated <-motifScorer(clusterInfoMS=clusterInfo)
      
      #write
      result<-motifWriter(clusterInfoMW=clusterAnnotated, chromMW=chrom, targetMW=targetEntry[j], arabSeqMW=arabSeq, 
                          fileStreamMW=fileStream, fileStreamBAMMW=fileStreamBAM, targetDataMW=targetMetachar[j])
      
      
      
      }
      
    }
    }
  
  
#read in results
  resultList = read.table(resultsOutput, header=TRUE)
  
  #sort results 
  attach(resultList)
  resultsList<-resultList[order(phasescore),]
  detach(resultList)
  fileStreamSorted <- file(resultsSorted, 'w')
  write.table(resultList, sep=' ', quote=FALSE, row.names=FALSE, fileStreamSorted)
  
  #close filestreams
  print("check15")
  
  close(fileStream)
  close(fileStreamBAM)
  close(fileStreamSorted)
  

  
}





