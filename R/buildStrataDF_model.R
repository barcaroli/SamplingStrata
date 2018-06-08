# ----------------------------------------------------
# Function to produce the "strata" dataframe
# starting from the available sampling frame
# taking into account anticipated variance
# Author: Giulio Barcaroli
# Date: May 2018
# ----------------------------------------------------
buildStrataDF <- function(dataset, model=NULL) {
    # stdev1 is for sampling data
    stdev1 <- function(x, w) {
        mx <- sum(x * w)/sum(w)
        sqrt(sum(w * (x - mx)^2)/(sum(w) - 1))
    }
    # stdev2 is for population data
    stdev2 <- function(x, w) {
        mx <- sum(x * w)/sum(w)
        sqrt(sum(w * (x - mx)^2)/(sum(w)))
    }
    colnames(dataset) <- toupper(colnames(dataset))
    nvarX <- length(grep("X", names(dataset)))
    nvarY <- length(grep("Y", names(dataset)))
    if (length(grep("WEIGHT", names(dataset))) == 1) {
        cat("\nComputations are being done on sampling data\n")
        stdev <- "stdev1"
    }
    if (length(grep("WEIGHT", names(dataset))) == 0) {
        dataset$WEIGHT <- rep(1, nrow(dataset))
        stdev <- "stdev2"
        cat("\nComputations are being done on population data\n")
    }
    #---------------------------------------------------------
    # Check the validity of the model
    if (!is.null(model)) {
      if (nrow(model) != nvarY) stop("A model for each Y variable must be specified")
      for (i in (1:nrow(model))) {
        if (!(model$type[i] %in% c("linear","loglinear"))) stop("Type of model for Y variable ",i,"misspecified")
        if (is.na(model$beta[i])) stop("beta for Y variable ",i,"must be specified")
        if (is.na(model$sig2[i])) stop("sig2 for Y variable ",i,"must be specified")
        if (model$type[i] == "linear" & is.na(model$gamma[i])) stop("gamma for Y variable ",i,"must be specified")
      }
    }
    #---------------------------------------------------------     
    numdom <- length(levels(as.factor(dataset$DOMAINVALUE)))
    stratatot <- NULL
    # create progress bar
    pb <- txtProgressBar(min = 0, max = numdom, style = 3)
    # begin domains cycle
    for (d in (1:numdom)) {
      Sys.sleep(0.1)
      # update progress bar
      setTxtProgressBar(pb, d)
		  dom <- levels(as.factor(dataset$DOMAINVALUE))[d]
		  domain <- dataset[dataset$DOMAINVALUE == dom, ]
        listX <- NULL
        namesX <- NULL
        for (i in 1:nvarX) {
            name <- paste("X", i, sep = "")
            namesX <- cbind(namesX, name)
            if (i < nvarX) 
                listX <- paste(listX, "domain$X", i, ",", sep = "") else listX <- paste(listX, "domain$X", i, sep = "")
        }
        listM <- NULL
        listS <- NULL
        for (i in 1:nvarY) {
            listM <- paste(listM, "M", i, ",", sep = "")
            listS <- paste(listS, "S", i, ",", sep = "")
        }
        stmt <- paste("domain$STRATO <- as.factor(paste(", listX, 
            ",sep='*'))", sep = "")
        eval(parse(text = stmt))
        for (i in 1:nvarY) {
            WEIGHT <- NULL
            STRATO <- NULL
            Y <- NULL
            stmt <- paste("Y <- domain$Y", i, "[!is.na(domain$Y", 
                i, ")]", sep = "")
            eval(parse(text = stmt))
            stmt <- paste("WEIGHT <- domain$WEIGHT[!is.na(domain$Y", 
                i, ")]", sep = "")
            eval(parse(text = stmt))
            stmt <- paste("STRATO <- domain$STRATO[!is.na(domain$Y", 
                i, ")]", sep = "")
            eval(parse(text = stmt))
            STRATO <- factor(STRATO)
            # Computation of M and S without model --------------------------
            if (is.null(model)) {
              stmt <- paste("M", i, " <- tapply(WEIGHT * Y,STRATO,sum) / tapply(WEIGHT,STRATO,sum)", sep = "")
              eval(parse(text = stmt))
              samp <- NULL
              stmt <- paste("samp <- domain[!is.na(domain$Y", i, "),]", sep = "")
              eval(parse(text = stmt))
              l.split <- split(samp, samp$STRATO, drop = TRUE)
              stmt <- paste("S", i, " <- sapply(l.split, function(df,x,w) ", 
                  stdev, "(df[,x],df[,w]), x='Y", i, "', w='WEIGHT')", 
                  sep = "")
              eval(parse(text = stmt))
            }
            # Computation of M and S with linear model --------------------------
            if (!is.null(model)) {
              if (model$type[i] == "linear") {
                stmt <- paste("M", i, " <- tapply(WEIGHT * Y * model$beta[", i, "],STRATO,sum) / tapply(WEIGHT,STRATO,sum)", sep = "")
                eval(parse(text = stmt))
                samp <- NULL
                stmt <- paste("samp <- domain[!is.na(domain$Y", i, "),]", sep = "")
                eval(parse(text = stmt))
                l.split <- split(samp, samp$STRATO, drop = TRUE)
                stmt <- paste("S", i, " <- sapply(l.split, function(df,x,w) ", 
                              stdev, "(df[,x],df[,w]), x='Y", i, "', w='WEIGHT')", 
                              sep = "")
                eval(parse(text=stmt))
                st <- paste("gammas <- tapply(Y^model$gamma[",i,"],STRATO,sum) / as.numeric(table(STRATO))",sep="")
                eval(parse(text=st))
                st <- paste("S",i," <- sqrt(S",i,"^2 * model$beta[",i,"]^2 + model$sig2[",i,"] * gammas)",sep="")
                eval(parse(text=st))   
              }
              # Computation of M and S with loglinear model --------------------------
              if (!is.null(model) & model$type[i] == "loglinear") {
                stmt <- paste("M", i, " <- tapply(WEIGHT * Y ^ model$beta[",i,"],STRATO,sum) / tapply(WEIGHT,STRATO,sum)", sep = "")
                eval(parse(text = stmt))
                ph <- 1  # remember to calculate ph as proportion of Y > 0
                st <- paste("S", i, " <- sqrt(ph * (( exp(model$sig2[", i, "])* 
                               tapply(WEIGHT * Y^(2*model$beta[", i, "]),STRATO,sum)/as.numeric(table(STRATO)) -
                               ph * (tapply(WEIGHT * Y^model$beta[", i, "],STRATO,sum)/as.numeric(table(STRATO)))^2)))",sep="")
                eval(parse(text = st))
              }
            }
            # ------------------------------------------------------------------------
            if (is.null(model)) eval(parse(text = stmt))
            stmt <- paste("stratirid <- unlist(attr(M", i, ",'dimnames'))", 
                sep = "")
            eval(parse(text = stmt))
            strati <- data.frame(X1 = levels(domain$STRATO))
            stmt <- paste("m <- data.frame(cbind(X1=stratirid,X2=M", 
                i, "))", sep = "")
            eval(parse(text = stmt))
            m <- merge(strati, m, by = c("X1"), all = TRUE)
            m$X2 <- as.character(m$X2)
            m$X2 <- as.numeric(m$X2)
            m$X2 <- ifelse(is.na(m$X2), 0, m$X2)
            stmt <- paste("M", i, " <- m$X2", sep = "")
            eval(parse(text = stmt))
            stmt <- paste("s <- data.frame(cbind(X1=stratirid,X2=S", 
                i, "))", sep = "")
            eval(parse(text = stmt))
            s <- merge(strati, s, by = c("X1"), all = TRUE)
            s$X2 <- as.character(s$X2)
            s$X2 <- as.numeric(s$X2)
            s$X2 <- ifelse(is.na(s$X2), 0, s$X2)
            stmt <- paste("S", i, " <- s$X2", sep = "")
            eval(parse(text = stmt))
        }
        N <- tapply(domain$WEIGHT, domain$STRATO, sum)
        STRATO <- domain$STRATO
        COST <- rep(1, length(levels(domain$STRATO)))
        CENS <- rep(0, length(levels(domain$STRATO)))
        DOM1 <- rep(as.character(dom), length(levels(domain$STRATO)))
        stmt <- paste("strata <- as.data.frame(cbind(STRATO=levels(STRATO),N,", 
            listM, listS, "COST,CENS,DOM1))")
        eval(parse(text = stmt))
        for (i in 1:nvarX) {
            stmt <- paste("strata$X", i, " <- rep(0, length(levels(domain$STRATO)))", 
                sep = "")
            eval(parse(text = stmt))
        }
        strata$STRATO <- as.character(strata$STRATO)
        for (i in 1:nrow(strata)) {
            strata[i, c(namesX)] <- unlist(strsplit(strata$STRATO[i], 
                "\\*"))
        }
        stratatot <- rbind(stratatot, strata)
    }  # end domain cycle
    close(pb)
    colnames(stratatot) <- toupper(colnames(stratatot))
    stratatot$DOM1 <- as.factor(stratatot$DOM1)
    write.table(stratatot, "strata.txt", quote = FALSE, sep = "\t", 
                dec = ".", row.names = FALSE)
    stratatot <- read.delim("strata.txt")
    options("scipen"=100)
    for (i in (1:nvarY)) {
      for (j in (1:nrow(stratatot))) {
        stmt <- paste("stratatot$M",i,"[j] <- ifelse(stratatot$M",i,"[j] == 0,0.000000000000001,stratatot$M",i,"[j])",sep="")
        eval(parse(text=stmt))
      }
    }
    write.table(stratatot, "strata.txt", quote = FALSE, sep = "\t", 
        dec = ".", row.names = FALSE)
    stratatot <- read.delim("strata.txt")
    cat("\nNumber of strata: ",nrow(stratatot))
    cat("\n... of which with only one unit: ",sum(stratatot$N==1))
    return(stratatot)
}
