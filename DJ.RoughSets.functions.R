
# Extend RoughSets library

# Author: Dariusz Jankowski
# License: For academic use only

require(RoughSets)

# Konwersja zagnieżdzownej nazwanej listy na data.frame
DJ.namedlistOfLists.as.data.frame <- function(list, colnames = c())
{
  rel = list
  rel.names = names(rel)
  
  d1 <- data.frame(L1=c(), stringsAsFactors=FALSE) 
  for(i in 1:length(rel.names))
  {
    d1 <- rbind(d1, data.frame(L1=rel.names[i])) 
  }
  
  d2 <- data.frame(L2=c(), stringsAsFactors=FALSE) 
  
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

# Example of input values: DM=DM.A; indx.rows = c(100:101); indx.columns = c(1:2); 
DJ.discernibility.mat.as.ktable.RST <- function(DM, indx.rows = NULL, indx.columns = NULL, format= 'markdown', format.args = NULL )
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
      # print(paste(row,col,length(item),item))
      if(length(item)>1)
      {
        dm[row,col] <- paste0("{",item[1])
        for(i in 2:length(item))  # print(cat(dm[row,col], ",", item[i]))
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
        dm[row,col] <- item  # paste0("{",item,"}")
      }
    }
    
    # dm[,col] = as.factor(dm[,col])
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

DJ.LU.extract.listOfObjects <- function(lowerOrUpperApproximation, decision.name, sort_result=TRUE, sort_decreasing=FALSE)
{
  result = unname(lowerOrUpperApproximation[[decision.name]])
  if (sort_result)
  {
    result = sort(result, decreasing=sort_decreasing)
  }
  result
}

# concat lists like sets
DJ.combine.lists <- function(listOfLists, sort_result=TRUE, sort_decreasing=FALSE)
{
  if(!any(class(listOfLists) %in% c("list"))) { stop("Invalid class of list") }
  result = unique(do.call(c, listOfLists))
  if (sort_result)
  {
    result = sort(result, decreasing=sort_decreasing)
  }
  result
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


DJ.ruleKeyValue.format <- function(key, value, output_format="LERS")
{
  
  if(output_format == "LERS")
  {
    start_c = "("
    middle_c = ", "
    end_c = ")"
  }
  else
  {
    start_c = " "
    middle_c = "="
    end_c = " "
  }
  
  str = paste0(start_c, key, middle_c, value, end_c)
  
  return(str)
}
  
  
DJ.rulePart1.format <- function(cols, l, output_format="LERS")
{
  if(output_format == "LERS")
  {
    start_rule = "" # nie trzeba, bo każda para klucz/wartość zawiera nawiasy
    start_c = "("
    middle_c = ","
    end_c = ")"
    and_c = " & "
    then_c = " -> "
  }
  else
  {
    start_rule = "IF "
    start_c = " "
    middle_c = "="
    end_c = " "
    and_c = " AND "
    then_c = " THEN "
  }
  
  str <- start_rule
  if(is.vector(l$idx))
  {
    str1 = DJ.ruleKeyValue.format(key=cols[ l$idx[1]  ],value=l$values[1], output_format = output_format)
    str <- paste0(str, str1)
    
    if(length(l$idx)>1)
    {
      for(i in 2:length(l$idx))
      {
        str1 = DJ.ruleKeyValue.format(key=cols[ l$idx[i]  ],value=l$values[i], output_format = output_format)
        str <- paste0(str, and_c, str1)
      }
    }
  }else{
    #  str <- paste0(str, cols[l$idx],"=", l$values)
    stop("rule $idx is not a vector")
  }
  
  str <- paste0(str, then_c)
  
  return(str)
  
}

DJ.rulePart2.format <- function(dec, l, output_format="LERS")
{
  if(output_format == "LERS")
  {
    or_c = " | "
  }
  else
  {
    or_c = " OR "
  }
  
  str <- ""
  if(is.vector(l$all_consequents))
  {
    str1 = DJ.ruleKeyValue.format(key=dec, value=l$all_consequents[1], output_format = output_format)
    str <- paste0(str, str1)
    
    if(length(l$all_consequents)>1)
    {
      for(i in 2:length(l$all_consequents))
      {
        str1 = DJ.ruleKeyValue.format(key=dec, value=l$all_consequents[i], output_format = output_format)
        str <- paste0(str, or_c, str1)
      }
    }
  }
  else
  {
    stop("rule $idx is not a vector")
  }
  
  return(str)
  
}


DJ.rule.toString.RST <- function(rule, show_conflict_rules=TRUE, output_format = "LERS")
{
  if(as.set(class(rule)) != set("RuleSetRST","list")) { stop("Błędny parametr") }
  if(length(rule)>1) { stop("Błędna liczba reguł")}
  #rule = rules[3] # only to debug
  
  cols <- attr(rule, "colnames")
  dec <- attr(rule, "dec.attr")
  l <- rule[[1]]
  str1 = DJ.rulePart1.format(cols, l, output_format = output_format)
  
  if(show_conflict_rules)
  {
    str2 = DJ.rulePart2.format(dec, l, output_format = output_format)
  }
  else
  {
    str2 = DJ.ruleKeyValue.format(key=dec ,value=l$consequent, output_format = output_format)
  }

  str <- paste0(str1, str2)
  
  #c("rule"= str, "support"= l$support, "laplace"= l$laplace )
  df <- as.data.frame(str)
  
  df <- cbind(df,  length(l$support) )
  df <- cbind(df,  l$laplace[1])
  colnames(df) <- c("rule", "supportSize", "laplace")
  rownames(df) <- c()
  
  return(df)
}


DJ.rules.toString.RST <- function(rules, output_format="LERS", show_conflict_rules=TRUE)
{
  if(is.null(rules)) {return(NULL)}
  if(as.set(class(rules)) != set("RuleSetRST","list")) { stop("Błędny parametr")}

  formatter=DJ.rule.toString.RST
  
  str <- formatter(rules[1], show_conflict_rules, output_format=output_format)
  if(length(rules)>1)
    for(i in 2:length(rules)) 
    {
      rule <- formatter(rules[i], show_conflict_rules, output_format=output_format)
      str <- rbind(str, rule)
    }
  
  df <- cbind(str, "RI.support"=RI.support(rules))
  df <- cbind(df, "RI.confidence"=RI.confidence(rules))
  rownames(df) <- c()
  
  return(df)
}
