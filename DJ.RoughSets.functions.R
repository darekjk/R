
# Inne funkcje pomocnicze
# Author: Dariusz Jankowski
# Licenacja użytkowania: na zajęciach SI - Politechnika Białostocka - semestr 2017/2018

require(caret)
require(data.table)
require(sets)
require(RoughSets)
require(knitr)
require(kableExtra)

# Konwersja zagnieżdzownej nazwanej listy na data.frame
DJ.namedlistOfLists.as.data.frame <- function(list, colnames = c())
{
  rel = list
  rel.names = names(rel)
  
  d1 <- data.frame(L1=c() ,stringsAsFactors=FALSE) 
  for(i in 1:length(rel.names))
  {
    d1 <- rbind(d1, data.frame(L1=rel.names[i])) 
  }
  
  #d1
  
  d2 <- data.frame(L2=c() ,stringsAsFactors=FALSE) 
  #d2
  for(i in 1:length(rel) ) 
  {
    s <- ""
    rel2 <- rel[[i]]
    for(j in 1:length(rel2) )
    { 
      s <- paste0(s, " ", rel2[j]  ) 
    }
    
    d2<- rbind(d2, data.frame(L2=s)) 
    
  }
  #d2
  d3 <- cbind(d1, d2)
  colnames(d3) <- colnames
  d3
}

DJ.IND.relation.as.data.frame <- function(rel, attr.cond.set.name ="B", colnames)
{
  rel.colnames <- paste( unlist(colnames), collapse=', ')
  rel.colnames <- c(paste0(attr.cond.set.name, "={", rel.colnames,"}"), "objects")
  
  DJ.namedlistOfLists.as.data.frame(rel, rel.colnames)
}

#DM=DM.A; indx.rows = c(100:101); indx.columns = c(1:2); format = "markdown"; format.args = NULL  
DJ.discernibility.mat.RST.as.ktable <- function(DM, indx.rows = NULL, indx.columns = NULL, format= 'markdown', format.args = NULL )
{
  if(is.null(indx.columns))indx.columns = 1:ncol(DM$disc.mat)
  if(is.null(indx.rows))indx.rows = 1:nrow(DM$disc.mat)
  library(knitr)
  dm <- NULL
  dm <- DM$disc.mat[indx.rows, indx.columns]
  
  for(col in 1:ncol(dm))
  {
    for(row in 1:nrow(dm))  
    {
      item <- dm[row,col]
      item <- item[[1]]
      #print(paste(row,col,length(item),item))
      if(length(item)>1)
      {
        dm[row,col] <- paste0("{",item[1])
        for(i in 2:length(item)) #print(cat(dm[row,col], ",", item[i]))
        {
          s = dm[row,col]
          dm[row,col] <-  paste0(s, ",", item[i])
        }
        s = dm[row,col]
        dm[row,col] <- paste0(s, "}")
      }
      else
      {
        # item = item[[1]]
        
        if(!is.null(item) && (is.na(item) || length(item)==0))
        { item <- "" }
        dm[row,col] <- item #paste0("{",item,"}")
      }
    }
    
    #dm[,col] = as.factor(dm[,col])
    dm[col,col] <- ""
  }
  dm <- as.data.frame(dm, row.names = 1:nrow(dm))
  colnames(dm) <- c(1:ncol(dm))
  kable(dm, row.names = T, format = format, format.args = format.args)
}

DJ.computeCore.from.discernibility.mat.RST <- function(DM)
{
  cat("Core:\n")
  cat(paste(unique(unlist(DM$disc.list[which(lengths(DM$disc.list)==1)]))))
}

DJ.LU.extract.listOfObjects <- function(lowerOrUpperApproximation, decision.name)
{
  sort(unname(lowerOrUpperApproximation[[decision.name]]))
}

DJ.aggregate.data.table.by.column <- function(data)
{
  stop("not implemented")
}

DJ.set.toString <-function(data, sep = ",")
{
  if(!any(class(data) %in% c("set", "gset", "cset"))) { stop("Invalid class of set") }
  data=ad; sep=","
  
  list <- as.list(data)
  str = ""
  for(i in 1:length(list))
  {
    str <- paste(str, as.character(list[i]), sep = sep)
  }
  substr(x = str, start = length(collapse)+1, stop = nchar(str))
}

DJ.attr_and_dec <- function(decisionTable, listOf.cond.attr)
{
  c(listOf.cond.attr, ncol(decisionTable))
}

DJ.rule.RST.toString <- function(rule)
{
  if(as.set(class(rule)) != set("RuleSetRST","list")) { stop("Błędny parametr") }
  if(length(rule)>1) { stop("Błędna liczba reguł")}
  #rule = rules[3] # only to debug
  
  cols <- attr(rule, "colnames")
  dec <- attr(rule, "dec.attr")
  l <- rule[[1]]
  str <- "IF "
  if(is.vector(l$idx))
  {
    str <- paste0(str, cols[ l$idx[1]  ],"=", l$values[1])
    
    if(length(l$idx)>1)
    {
      for(i in 2:length(l$idx))
      {
        str <- paste0(str, " AND ", cols[ l$idx[i] ],"=", l$values[i])
      }
    }
  }else{
    #  str <- paste0(str, cols[l$idx],"=", l$values)
    stop("rule $idx is not a vector")
  }
  
  str <- paste0(str, " THEN ", dec,"=", l$consequent)
  
  #c("rule"= str, "support"= l$support, "laplace"= l$laplace )
  df <- as.data.frame(str)
  
  df <- cbind(df,  length(l$support) )
  df <- cbind(df,  l$laplace[1])
  colnames(df) <- c("rule", "supportSize", "laplace")
  rownames(df) <- c()
  
  df
}

DJ.rules.RST.toString <- function(rules)
{
  if(as.set(class(rules)) != set("RuleSetRST","list")) { stop("Błędny parametr")}
  
  
  str <- DJ.rule.RST.toString(rules[1])
  if(length(rules)>1)
    for(i in 2:length(rules)) 
    {
      rule <- DJ.rule.RST.toString(rules[i])
      str <- rbind(str, rule)
    }
  
  df <- cbind(str, "RI.support"=RI.support(rules))
  df <- cbind(df, "RI.confidence"=RI.confidence(rules))
  rownames(df) <- c()
  
  df
}
