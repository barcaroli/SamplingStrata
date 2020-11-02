# ----------------------------------------------------
# Function to produce the "strata" dataframe
# starting from the available sampling frame
# taking into account anticipated variance
# Author: Giulio Barcaroli
# Date: October 2019
# ----------------------------------------------------
buildStrataDF <- function(dataset, 
                          model=NULL, 
                          progress=TRUE,
                          verbose=TRUE) {
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
    # stdev3 is for spatial models (part I)
    stdev3 <- function(Y, W, beta1, beta2) {
      a <- as.matrix(t(c(beta1,beta2)))
      b <- cov(cbind(Y, W))
      c <- as.matrix(c(beta1,beta2))
      sqrt (a %*% b %*% c)
    }
    # stdev4 is for spatial models (part II)
    stdev4 <- function(df,var_eps,range,gamma,i) {
      st <- paste("Y <- df$Y",i,sep="")
      eval(parse(text=st))
      dist <- sqrt((outer(df$LON,df$LON,"-"))^2+(outer(df$LAT,df$LAT,"-"))^2)
      pred <- beta1*Y + beta2*W
      var_ntimes <- rep(var_eps,nrow(df))
      var_ntimes <- var_ntimes*Y^(2*gamma)
      sum_couples_var <- as.matrix(outer(var_ntimes,var_ntimes,"+"))
      # prod_couples_std <- as.matrix(outer(sqrt(var_ntimes),sqrt(var_ntimes),"*"))
      prod_couples_std <- sqrt(as.matrix(outer(var_ntimes,var_ntimes, "*")))
      spatial_autocovariance <- prod_couples_std * exp(-1*dist/(range+0.0000001))
      D2<-sum_couples_var-2*spatial_autocovariance
      sum(sum_couples_var) / (2*nrow(df)^2)
      sum(2*spatial_autocovariance) / (2*nrow(df)^2)
      sum(sum_couples_var-2*spatial_autocovariance) / (2*nrow(df)^2)
      var2 <- sum(D2)/(2*nrow(df)^2) 
      sqrt (var2)
    }
    # cov1 is for spatial models (part III)
    cov1 <- function(df,psill,range,Y,W,beta1,beta2) {
      preds <- beta1 * Y + beta2 * W
      dist <- sqrt((outer(df$LON, df$LON, "-"))^2 + (outer(df$LAT, df$LAT, "-"))^2)
      std_eps_ntimes <- sqrt(rep(psill, nrow(df)))
      v <- var(preds)
      v <- ifelse(is.na(v),0,v)
      std_pred_ntimes <- sqrt(rep(v,nrow(df)))
      prod_couples_std <- as.matrix(outer(std_eps_ntimes,std_pred_ntimes ),"*")
      spatial_autocovariance <- prod_couples_std * exp(-1 * dist/(range + 1e-07))
      D2 <-  2 * spatial_autocovariance
      var2 <- sum(D2)/(2 * nrow(df)^2)
      sqrt(var2)
    }
    colnames(dataset) <- toupper(colnames(dataset))
    # if (is.factor(dataset$DOMAINVALUE)) levels(dataset$DOMAINVALUE) <- levels(droplevels(dataset$DOMAINVALUE))
    nvarX <- length(grep("X", names(dataset)))
    nvarY <- length(grep("Y", names(dataset)))
    if (length(grep("WEIGHT", names(dataset))) == 1) {
        if (verbose == TRUE) {
          cat("\nComputations are being done on sampling data\n")
        }
          stdev <- "stdev1"
    }
    if (length(grep("WEIGHT", names(dataset))) == 0) {
        dataset$WEIGHT <- rep(1, nrow(dataset))
        stdev <- "stdev2"
        if (verbose == TRUE) {
          cat("\nComputations are being done on population data\n")
        }
      }
    #---------------------------------------------------------
    # Check the validity of the model
    if (!is.null(model)) {
      if (nrow(model) != nvarY) stop("A model for each Y variable must be specified")
      for (i in (1:nrow(model))) {
        if (!(model$type[i] %in% c("linear","loglinear","spatial"))) stop("Type of model for Y variable ",i,"misspecified")
        if (is.na(model$beta[i])) stop("beta for Y variable ",i,"must be specified")
        if (is.na(model$sig2[i])) stop("sig2 for Y variable ",i,"must be specified")
        if (model$type[i] == "spatial") {
          if (is.na(model$beta2[i])) stop("beta2 for Y variable ",i,"must be specified")
          if (is.na(model$range[i])) stop("range for Y variable ",i,"must be specified")
          if (is.null(dataset$LON) | is.null(dataset$LON) ) stop("Missing coordinates on sampling frame")
        }
        if (model$type[i] == "linear" & is.na(model$gamma[i])) stop("gamma for Y variable ",i,"must be specified")
      }
      
    }
    #--------------------------------------------------------- 
    dataset$DOMAINVALUE <- factor(dataset$DOMAINVALUE)
    # dataset$DOMAINVALUE <- droplevels(dataset$DOMAINVALUE)
    numdom <- length(levels(dataset$DOMAINVALUE))
#    numdom <- length(unique(dataset$DOMAINVALUE))
    stratatot <- NULL
    # create progress bar
    if (progress == TRUE) pb <- txtProgressBar(min = 0, max = numdom, style = 3)
    # begin domains cycle
    # dataset$DOMAINVALUE <- as.numeric(dataset$DOMAINVALUE)
    # for (d in unique(dataset$DOMAINVALUE)) {
    # dataset$DOMAINVALUE <- as.numeric(dataset$DOMAINVALUE)
    doms <- 0
    for (d in (levels(dataset$DOMAINVALUE))) {
      doms <- doms + 1
      if (progress == TRUE) Sys.sleep(0.1)
      # update progress bar
      if (progress == TRUE) setTxtProgressBar(pb, doms)
      # dom <- unique(dataset$DOMAINVALUE)[d]
		  # dom <- levels(as.factor(dataset$DOMAINVALUE))[d]
      dom <- d
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
        if (!is.null(dataset$COST)) {
          cost <- tapply(domain$WEIGHT * domain$COST,domain$STRATO,sum) / tapply(domain$WEIGHT,domain$STRATO,sum)
        }
        for (i in 1:nvarY) {
            WEIGHT <- NULL
            STRATO <- NULL
            Y <- NULL
            stmt <- paste("Y <- domain$Y", i, "[!is.na(domain$Y", 
                i, ")]", sep = "")
            eval(parse(text = stmt))
            W <- NULL
            stmt <- paste("W <- domain$W", i, "[!is.na(domain$W", 
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
            if (!is.null(model)) {
              # Computation of M and S with linear model --------------------------
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
                st <- paste("gammas <- tapply(Y^(model$gamma[",i,"]*2),STRATO,sum) / as.numeric(table(STRATO))",sep="")
                eval(parse(text=st))
                st <- paste("S",i," <- sqrt(S",i,"^2 * model$beta[",i,"]^2 + model$sig2[",i,"] * gammas)",sep="")
                eval(parse(text=st))   
              }
              # Computation of M and S with loglinear model --------------------------
              if (!is.null(model) & model$type[i] == "loglinear") {
                stmt <- paste("M", i, " <- tapply(WEIGHT * Y ^ model$beta[",i,"],STRATO,sum) / tapply(WEIGHT,STRATO,sum)", sep = "")
                eval(parse(text = stmt))
                #----------------------------------------------------------
                positiv <- function(x,w) {
                  sum(x > 0) / length(x)
                }
                samp <- NULL
                stmt <- paste("samp <- domain[!is.na(domain$Y", i, "),]", sep = "")
                eval(parse(text = stmt))
                l.split <- split(samp, samp$STRATO, drop = TRUE)
                stmt <- paste("ph <- sapply(l.split, function(df,x,w) positiv(df[,x],df[,w]), x='Y", i, "', w='WEIGHT')", sep = "")
                eval(parse(text=stmt))
                #----------------------------------------------------------
                # ph <- 1  
                st <- paste("S", i, " <- sqrt(ph * (( exp(model$sig2[", i, "])* 
                               tapply(WEIGHT * Y^(2*model$beta[", i, "]),STRATO,sum)/as.numeric(table(STRATO)) -
                               ph * (tapply(WEIGHT * Y^model$beta[", i, "],STRATO,sum)/as.numeric(table(STRATO)))^2)))",sep="")
                eval(parse(text = st))
              }
              # Computation of M and S with spatial model-----------------------------------------
              if (model$type[i] == "spatial") {
                stmt <- paste("M", i, " <- tapply(WEIGHT * (Y * model$beta[",i,"] + W * model$beta2[",i,"]) ,STRATO,sum) / tapply(WEIGHT,STRATO,sum)", sep = "")
                eval(parse(text = stmt))
                samp <- NULL
                stmt <- paste("samp <- domain[!is.na(domain$Y", i, "),]", sep = "")
                eval(parse(text = stmt))
                l.split <- split(samp, samp$STRATO, drop = TRUE)
                #-- PART I ---------------
                stdev = "stdev3"
                fitting <- model$fitting[i]
                beta1 <- model$beta[i]
                beta2 <- model$beta2[i]
                # st <- paste("gammas <- tapply(Y^model$gamma[",i,"],STRATO,sum) / as.numeric(table(STRATO))",sep="")
                # eval(parse(text=st))
                stdev <- "stdev3"
                stmt <- paste("sd1 <- sapply(l.split, function(df,y,w) ",
                              stdev, "(df[,y],df[,w],beta1,beta2), y = 'Y",i,"',w = 'W",i,"')",
                              sep = "")
                eval(parse(text=stmt))
                #-- PART II ---------------
                stdev <- "stdev4"
                sig2_eps <- model$sig2[i]
                range <- model$range[i]
                gamma <- model$gamma[i]
                stmt <- paste("sd2 <- sapply(l.split, function(df) ",
                              stdev, "(df,sig2_eps,range,gamma,i))",sep = "")
                eval(parse(text=stmt))
                # sd2 <- sapply(l.split, function(df) stdev4(Y,model$sig2[i],model$range[i],model$gamma[i]))
                # stdev <- "stdev4"
                # stmt2 <- paste("sd2 <- sapply(l.split, function(df) ",
                #               stdev, "(df,psill,range))",
                #               sep = "")
                # eval(parse(text=stmt2))
                #-- PART III ---------------
                # stmt <- paste("cov1 <- sapply(l.split, function(df,y,w) ",
                #                "cov1(df,psill,range,df[,y],df[,w],beta1,beta2), y = 'Y",i,"',w = 'W",i,"')",
                #                sep = "")
                # eval(parse(text=stmt))
                #-- TOTAL S ---------------
                # st <- paste("S",i," <- sqrt(sd1^2 + sd2^2 + cov1^2)",sep="")
                # st <- paste("S",i," <- sqrt(sd1^2 + sd2^2)",sep="")
                # st <- paste("S",i," <- sqrt(sd1^2 + sd2^2)",sep="")
                # psill2 <- model$sig2_2[i]
                # range2 <- model$range_2[i]
                # stmt <- paste("cov1 <- sapply(l.split, function(df,y,w) ",
                #                "cov1(df,psill2,range2,df[,y],df[,w],beta1,beta2), y = 'Y",i,"',w = 'W",i,"')",
                #                sep = "")
                # eval(parse(text=stmt))
                #-- TOTAL S ---------------
                # st <- paste("gammas <- tapply(Y^model$gamma[",i,"],STRATO,sum) / as.numeric(table(STRATO))",sep="")
                # eval(parse(text=st))
                # st <- paste("S",i," <- sqrt(sd1^2/fitting + (sd2^2 + cov1^2) * gammas)",sep="")
                # st <- paste("S",i," <- sqrt((sd1^2/fitting) + sd2^2 * gammas)",sep="")
                st <- paste("S",i," <- sqrt((sd1^2/fitting) + sd2^2)",sep="")
                eval(parse(text=st))
              }
            }
            # ------------------------------------------------------------------------
            if (is.null(model)) eval(parse(text = stmt))
            stmt <- paste("stratirid <- unlist(attr(M", i, ",'dimnames'))", 
                sep = "")
            eval(parse(text = stmt))
            strati <- data.frame(X1 = levels(domain$STRATO), stringsAsFactors = TRUE)
            stmt <- paste("m <- data.frame(cbind(X1=stratirid,X2=M", 
                i, "), stringsAsFactors = TRUE)", sep = "")
            eval(parse(text = stmt))
            m <- merge(strati, m, by = c("X1"), all = TRUE)
            m$X2 <- as.character(m$X2)
            m$X2 <- as.numeric(m$X2)
            m$X2 <- ifelse(is.na(m$X2), 0, m$X2)
            stmt <- paste("M", i, " <- m$X2", sep = "")
            eval(parse(text = stmt))
            stmt <- paste("s <- data.frame(cbind(X1=stratirid,X2=S", 
                i, "), stringsAsFactors = TRUE)", sep = "")
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
        if (is.null(dataset$COST)) COST <- rep(1, length(levels(domain$STRATO)))
        if (!is.null(dataset$COST)) COST <- cost
        CENS <- rep(0, length(levels(domain$STRATO)))
        DOM1 <- rep(as.character(dom), length(levels(domain$STRATO)))
        stmt <- paste("strata <- as.data.frame(cbind(STRATO=levels(STRATO),N,", 
            listM, listS, "COST,CENS,DOM1), stringsAsFactors = TRUE)")
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
    if (progress == TRUE) close(pb)
    colnames(stratatot) <- toupper(colnames(stratatot))
    stratatot$DOM1 <- as.factor(stratatot$DOM1)
    # write.table(stratatot, "strata.txt", quote = FALSE, sep = "\t",
    #             dec = ".", row.names = FALSE)
    # stratatot <- read.delim("strata.txt")
    # unlink("strata.txt")
    options("scipen"=100)
    indx <- sapply(stratatot, is.factor)
    stratatot[indx] <- lapply(stratatot[indx], function(x) as.numeric(as.character(x)))
      for (j in (1:nvarX)) {
        stmt <- paste("stratatot$X",j," <- as.numeric(stratatot$X",j,")",sep="")
        eval(parse(text=stmt))
      }
      for (j in (1:nrow(stratatot))) {
        stmt <- paste("stratatot$M",i,"[j] <- ifelse(stratatot$M",i,"[j] == 0,0.000000000000001,stratatot$M",i,"[j])",sep="")
        eval(parse(text=stmt))
      }
    # }
    # if (writeFiles == TRUE )
    # write.table(stratatot, "strata.txt", quote = FALSE, sep = "\t", 
    #     dec = ".", row.names = FALSE)
    # stratatot <- read.delim("strata.txt")
    if (verbose == TRUE) {
      cat("\nNumber of strata: ",nrow(stratatot))
      cat("\n... of which with only one unit: ",sum(stratatot$N==1))
    }
    # Eliminate rare cases
    # for (i in (1:nvarY)) {
    #   st <- paste0("stratatot$S",i," <- ifelse(stratatot$M",i," < 0.01,0,stratatot$S",i,")")
    #   eval(parse(text=st))
    # }
    return(stratatot)
}
