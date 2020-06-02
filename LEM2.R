#an auxiliary function for computing laplace estimate of rule's confidence
laplaceEstimate <- function(rule, dataS, clsVec, uniqueCls, suppIdx = NULL) {
  if (is.null(suppIdx)) {
    if (length(rule$idx) > 1) {
      dataS <- do.call(paste, dataS[rule$idx])
    }
    else {
      dataS <- as.character(dataS[[rule$idx]])
    }
    
    tmpValues = paste(rule$values, collapse = " ")
    suppIdx = which(dataS == tmpValues)
  }

  clsFreqs = table(clsVec[suppIdx])
  maxIdx = which.max(clsFreqs)
  nOfCorrect = clsFreqs[maxIdx]
  
  rule$consequent = names(clsFreqs)[maxIdx]
  
  real_cls = clsFreqs[clsFreqs>0]  # usuń klasy, które nie są powiązane z regułą
  rule$all_consequents = names(real_cls) 
  rule$all_clsFreqs = real_cls

  rule$support = as.integer(suppIdx)
  rule$laplace = (nOfCorrect + 1)/(length(suppIdx) + length(uniqueCls))

  return(rule)
}


# Computation of a covering of a lower approximation of a concept by decision rules using the LEM2 algorithm.
# It is an auxiliary function.
computeLEM2covering <- function(concept, attributeValuePairs, decisionValues, uniqueCls) {
  
  if(length(concept) == 0) stop("Empty lower approximation of a decision class.")
  uncoveredConcept = concept
  rules = list()
  
  while(length(uncoveredConcept) > 0) {
    support = 1:length(decisionValues)
    rules[[length(rules) + 1]] = list(idx = integer(0), values = character(), support = support)
    tmpAttributeValuePairs = attributeValuePairs
    selectedAttributeValuePairs = list()
    
    while(!all(support %in% concept)) {
      tmpSupp = intersect(support, uncoveredConcept)
      correctLaplaceVec = sapply(tmpAttributeValuePairs,
                                 function(rule) {
                                   nOfCorrect = length(intersect(rule$support, tmpSupp));
                                   if(nOfCorrect > 0) laplace = nOfCorrect + rule$laplace  #(nOfCorrect + 1)/(length(tmpSupp) + length(uniqueCls))   #(length(rule$support) + length(uniqueCls))
                                   else laplace = 0;
                                   return(laplace)
                                 })
      tmpIdx = which.max(correctLaplaceVec)
      tmpRule = tmpAttributeValuePairs[[tmpIdx]]
      selectedAttributeValuePairs[[length(selectedAttributeValuePairs) + 1]] = tmpRule
      tmpAttributeValuePairs = tmpAttributeValuePairs[-tmpIdx]
      support = intersect(support, tmpRule$support)
    }
    rm(tmpSupp, tmpIdx, tmpRule, tmpAttributeValuePairs)
    
    toRmIdx = integer(0)
    if(length(selectedAttributeValuePairs) > 1) {
      for(i in 1:length(selectedAttributeValuePairs)) {
        suppList = lapply(selectedAttributeValuePairs[-c(toRmIdx, i)], function(x) x$support)
        if(length(suppList) > 0) {
          tmpSupport = Reduce(intersect, suppList)
        } else {
          tmpSupport = -1
        }
        if(all(tmpSupport %in% concept)) toRmIdx = c(toRmIdx, i)
      }
    }
    
    if(length(toRmIdx) > 0) {
      selectedAttributeValuePairs = selectedAttributeValuePairs[-toRmIdx]
      suppList = lapply(selectedAttributeValuePairs, function(x) x$support)
      support = Reduce(intersect, suppList)
    }
    
    idxVec = sapply(selectedAttributeValuePairs, function(x) x$idx)
    valuesVec = sapply(selectedAttributeValuePairs, function(x) x$values)
    rules[[length(rules)]]$support = support
    rules[[length(rules)]]$values = valuesVec
    rules[[length(rules)]]$idx = idxVec
    
    uncoveredConcept = setdiff(uncoveredConcept, support)
  }
  
  toRmIdx = integer(0)
  for(i in 1:length(rules)) {
    suppList = lapply(rules[-c(toRmIdx, i)], function(x) x$support)
    tmpSupport = Reduce(intersect, suppList)
    if(all(concept %in% tmpSupport)) toRmIdx = c(toRmIdx, i)
  }
  
  if(length(toRmIdx) > 0) {
    rules = rules[-toRmIdx]
  }
  
  return(rules)
}


# Modified version of function RI.LEM2Rules.RST from original RoughSets library
DJ.RI.LEM2Rules.byConcepts.RST <- function(decision.table, concepts)  {
  
  if (!inherits(decision.table, "DecisionTable")) {
    stop("Provided data should inherit from the \'DecisionTable\' class.")
  }
  
  if(is.null(attr(decision.table, "decision.attr"))) {
    stop("A decision attribute is not indicated.")
  } else decIdx = attr(decision.table, "decision.attr")
  
  if(!all(attr(decision.table, "nominal.attrs"))) {
    stop("Some of the attributes are numerical.
         Discretize attributes before calling RST-based rule induction methods.")
  }
  
  clsVec <- decision.table[,decIdx]
  uniqueCls <- unique(clsVec)
  decisionName = colnames(decision.table)[decIdx]
  clsFreqs <- table(clsVec)
  
  descriptorsList = attr(decision.table, "desc.attrs")[-decIdx]
  descriptorCandidates = list()
  for (i in 1:length(descriptorsList)) {
    descriptorCandidates = c(descriptorCandidates,
                             lapply(descriptorsList[[i]],
                                    function(v, x) return(list(idx = x, values = v)), i))
  }
  
  attributeValuePairs = lapply(descriptorCandidates, laplaceEstimate,
                               decision.table[,-decIdx], clsVec, uniqueCls)
  rm(descriptorsList, descriptorCandidates)
  
  rules = list()
  
  for(i in 1:length(concepts)) {
    concept = as.integer(concepts[[i]])
    if(is.null(concept) | length(concept)==0) {next}
    rules[[i]] = computeLEM2covering(concept, attributeValuePairs, clsVec, uniqueCls)
  }
  
  if (length(rules)==0){return(NULL)}
  
  
  rules = unlist(rules, recursive = FALSE)
  
  rules = lapply(rules, function(x) laplaceEstimate(list(idx = x$idx, values = x$values),
                                                    decision.table, clsVec, uniqueCls, suppIdx = x$support))
  
  attr(rules, "uniqueCls") <- as.character(sort(uniqueCls))
  attr(rules, "supportDenominator") <- nrow(decision.table)
  attr(rules, "clsProbs") <- clsFreqs/sum(clsFreqs)
  attr(rules, "majorityCls") <- as.character(sort(uniqueCls)[which.max(clsFreqs)])
  attr(rules, "method") <- "LEM2Rules"
  attr(rules, "dec.attr") <- decisionName
  attr(rules, "colnames") <- colnames(decision.table)[-decIdx]
  
  rules = ObjectFactory(rules, classname = "RuleSetRST")
  
  return(rules);
}

