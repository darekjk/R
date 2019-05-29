# Zachowaj kolejność ładowania skryptów:
# library(RoughSets)
# source("DJ.RoughSets.bugFixes.R") 

require(RoughSets)

### BUG FIX funkcji computeCore
# Funkcja z pakietu RoughSets, używana przez funkcję comuteCore, dostepna tylko z poziomu źródła
convertCNFtoDNF <- function(CNFclauses){
  
  if(length(CNFclauses) > 1) {
    ## sort the clauses by their length
    CNFlengths = sapply(CNFclauses, length)
    if(any(CNFlengths == 0)) {
      CNFclauses = CNFclauses[CNFlengths > 0]
      CNFlengths = CNFlengths[CNFlengths > 0]
    }
    CNFclauses = CNFclauses[order(CNFlengths)]
    
    ## eliminate unnecessary clauses for efficiency
    tmpCNFlength = length(CNFclauses)
    j = 2
    while(j <= tmpCNFlength){
      tmpIdx = sapply(CNFclauses[j:tmpCNFlength], function(x,y) all(y %in% x), CNFclauses[[j-1]])
      CNFclauses = CNFclauses[c(rep(T, j-1), !tmpIdx)]
      j = j + 1
      tmpCNFlength = length(CNFclauses)
    }
    
    ## convert to DNF form - start from the first CNF clause
    tmpDNF = CNFclauses[[1]]
    for(i in 2:length(CNFclauses))  {
      ## expand two clauses into possible DNFs
      tmpDNF = expand.grid(tmpDNF, CNFclauses[[i]],
                           KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
      tmpDNF = split(tmpDNF, 1:nrow(tmpDNF))
      ## take only those which are unique
      tmpDNF = lapply(tmpDNF, function(x) unique(unlist(x)))
      tmpLengths = sapply(tmpDNF, length)
      tmpDNF = tmpDNF[order(tmpLengths)]
      ## eliminate unnecessary clauses for efficiency
      tmpDNFlength = length(tmpDNF)
      j = 2
      while(j <= tmpDNFlength)  {
        tmpIdx = sapply(tmpDNF[j:tmpDNFlength], function(x,y) all(y %in% x), tmpDNF[[j-1]])
        tmpDNF = tmpDNF[c(rep(T, j-1), !tmpIdx)]
        j = j + 1
        tmpDNFlength = length(tmpDNF)
      }
    }
  } else {
    tmpDNF = CNFclauses
  }
  
  DNFclauses = lapply(tmpDNF, function(x) x[order(x)])
  names(DNFclauses) = paste("reduct", 1:length(DNFclauses), sep = "")
  return(DNFclauses)
}

DJ.computeCore = function(reductList) {
  
  if(length(reductList) > 0) {
    stopFlag = FALSE
    i = 2
    core = reductList[[1]]
    
    # BUG FIXED: 27.05.2018 (Dariusz Jankowski)
    while((!stopFlag) && (i<=length(reductList))) {
    #  print(i)
     # paste0("przed:",core)
     # paste0("%in%", core[core %in% reductList[[i]]] )
      core = core[core %in% reductList[[i]]]
      #paste0("po:",core)
      i = i + 1
      if(length(core) < 1) stopFlag = TRUE
    }
  } else stop("empty reduct list")
  
  return(core)
}

FS.all.reducts.computation = function(discernibilityMatrix, ...)
{
  reducts = RoughSets::FS.all.reducts.computation(discernibilityMatrix,...)
  reductList <- convertCNFtoDNF(discernibilityMatrix$disc.list)
  # BUG FIXED: wyliczanie core
  reducts$core = DJ.computeCore(reductList)
  reducts
}
