
library(readxl)

targetdata <-
  read.csv("Rdatafiles/downregulatedwusgenes_cyclo.csv",
           header = TRUE,
           sep = ',')

targets <- targetdata[2]
targetsul <- unlist(targets)
targetEntry <- as.character(targetsul)

resultsdata <-
  read_excel("Rdatafiles/coreMotifOutputDOWN_cyclo.xlsx")

results <- resultsdata[2]
resultsul <- unlist(results)
resultsEntry <- as.character(resultsul)




foundlist <-list()
notfoundlist <-list()


for(j in 1:length(targetEntry))
{

  
  searchresult<-grep(targetEntry[j], resultsEntry)
  
  #find target gene in Arabidopsis gene list
  if (length(searchresult > 0))
    foundlist<-append(foundlist,targetEntry[j])
  
  else 
    notfoundlist<-append(notfoundlist,targetEntry[j])
  
  
}

lapply(foundlist, write, "foundlistdown.txt", append=TRUE, ncolumns=1000)

lapply(notfoundlist, write, "notfoundlistdown.txt", append=TRUE, ncolumns=1000)

