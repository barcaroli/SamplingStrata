# --------------------------------------------------
# assignStrataLabel
# Function to assign the optimized strata labels
# to new sampling units in the frame
# 18 April 2019
# Giulio Barcaroli
# --------------------------------------------------

assignStrataLabel <- function(dataset, s) 
  {
  colnames(dataset) <- toupper(colnames(dataset))
  colnames(s) <- toupper(colnames(s))
  # check on dataset content
  if (!"DOMAINVALUE" %in% colnames(dataset)) stop("No DOMAINVALUE variable in dataset")
  nvarX_1 <- length(grep("LOWER",colnames(s)))
  nvarX_2 <- length(grep("X",colnames(dataset))) 
  if (nvarX_1 != nvarX_2)  stop("Number of X variables in dataset not compatible with strata structure")
  df <- NULL  
  for (dom in unique(dataset$DOMAINVALUE)) {
    for (strat in unique(s$STRATUM[s$DOMAIN==dom])) {
      cat("\nDomain ",dom,"  Stratum ",strat)
      stmt <- "strato <- dataset[dataset$DOMAINVALUE == dom & ("
      # conditions for inclusion
      for (i in (1:nvarX_1)) {
        if (i != nvarX_1) {
          stmt <- paste(stmt,"(dataset$X",i," >= s[s$DOMAIN==dom & s$STRATUM==strat,'LOWER_X",i,"'] & dataset$X",i," <= s[s$DOMAIN==dom & s$STRATUM==strat,'UPPER_X",i,"']) | ",sep="")
      }                    
        if (i == nvarX_1) {
          stmt <- paste(stmt,"(dataset$X",i," >= s[s$DOMAIN==dom & s$STRATUM==strat,'LOWER_X",i,"'] & dataset$X",i," <= s[s$DOMAIN==dom & s$STRATUM==strat,'UPPER_X",i,"']) )",sep="")
        }
      }
      # conditions for exclusion
      stmt <- paste(stmt," & (")
      for (i in (1:nvarX_1)) {
        if (i != nvarX_1) {
          stmt <- paste(stmt,"!(dataset$X",i," < s[s$DOMAIN==dom & s$STRATUM==strat,'LOWER_X",i,"']) & !(dataset$X",i," > s[s$DOMAIN==dom & s$STRATUM==strat,'UPPER_X1']) & ",sep="")
        }                    
        if (i == nvarX_1) {
          stmt <- paste(stmt,"!(dataset$X",i," < s[s$DOMAIN==dom & s$STRATUM==strat,'LOWER_X",i,"']) & !(dataset$X",i," > s[s$DOMAIN==dom & s$STRATUM==strat,'UPPER_X",i,"']) ), ]",sep="")
        }
      }
      eval(parse(text=stmt))
      strato$STRATO <- paste(dom,strat,sep="***")
      stratomin <- strato[!(strato$ID %in% df$ID),]
      df <- rbind(df,stratomin)
    }
  }
  # t <- as.data.frame(table(df$STRATO))
  # write.table(t,"strati_out.csv",sep=";",row.names = F)

  return(df)
}

