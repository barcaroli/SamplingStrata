# ==============================================================================
# File:      aggrStrata2.R
# Functions: aggrStrata2() [exported]
#
# Purpose:
#   Builds the "strata" dataframe of a domain directly from the sampling frame,
#   given the stratum label of each unit: population size, means and standard
#   deviations of the Y variables in each stratum. If a "model" dataframe is
#   given, the anticipated variance (linear or loglinear model, with
#   heteroscedasticity) is used instead of the observed one. Used by the
#   "continuous" method (optimizeStrata2()).
# ==============================================================================

# ----------------------------------------------------
# Function to produce the "strata" dataframe
# starting from the available sampling frame
# taking into account anticipated variance
# Author: Giulio Barcaroli
# Date: August 2019
# ----------------------------------------------------
#' @title Builds the "strata" dataframe containing information on target variables Y's distributions in the different strata, starting from a frame
#'
#' @description
#' This function builds the dataframe "strata" considering as input a given domain in the sampling frame. In case a dataframe "model" is given, the anticipated variance in the strata for each target variable is calculated
#'
#' @param dataset This is the name of the dataframe containing the sample data, or the frame data.
#'   It is strictly required that auxiliary information is organised in variables named
#'   as X1, X2, ... , Xm (there should be at least one of them) and the target variables
#'   are denoted by Y1, Y2, ... , Yn.
#'   In addition, in case of sample data, a variable named 'WEIGHT' must be present in the dataframe,
#'   containing   the weigths associated to each sampling unit
#' @param model Dataframe with the parameters of the model(s) to be used to calculate
#'   anticipated variance.
#'   Default is null.
#' @param vett vector of values indicating how the units in the dataset must
#'   be aggregated in strata.
#' @param dominio Value indicating the domain in the dataset to be processed.
#'
#' @return
#' A dataframe containing strata
#'
#' @author
#' Giulio Barcaroli
#'
#' @examples
#' data(swissmunicipalities)
#' swissmunicipalities$id <- seq_len(nrow(swissmunicipalities))
#' frame <- buildFrameDF(swissmunicipalities, id = "id", X = "POPTOT",
#'                       Y = c("Pop020", "Pop2040"), domainvalue = "REG")
#' # strata of domain 1: three classes of POPTOT
#' dom1 <- frame[frame$domainvalue == 1, ]
#' vett <- var.bin(dom1$X1, bins = 3)
#' strata <- aggrStrata2(dataset = dom1,
#'                       model = NULL,
#'                       vett = vett,
#'                       dominio = 1)
#' strata
#'
#' @keywords survey
#' @export
aggrStrata2 <- function(dataset, 
                          model=NULL, 
                          vett,
                          dominio) {
    # numbers are converted to character and back (stratum labels, X values):
    # avoid scientific notation only inside this function, restoring the
    # user's option on exit
    op <- options(scipen = 100)
    on.exit(options(op), add = TRUE)
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
    stdev4 <- function(df,psill,range) {
      dist <- sqrt((outer(df$LON,df$LON,"-"))^2+(outer(df$LAT,df$LAT,"-"))^2)
      var_ntimes <- rep(psill,nrow(df))
      prod_couples_std <- as.matrix(outer(sqrt(var_ntimes),sqrt(var_ntimes),"*"))
      sum_couples_var <- as.matrix(outer(var_ntimes,var_ntimes,"+"))
      spatial_autocovariance <- prod_couples_std * exp(-1*dist/(range+0.0000001))
      D2<-sum_couples_var-2*spatial_autocovariance
      var2 <- sum(D2)/(2*nrow(df)^2)
      sqrt (var2)
    }
    # covar is for spatial models (part III)
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
           stdev <- "stdev1"
    }
    if (length(grep("WEIGHT", names(dataset))) == 0) {
        dataset$WEIGHT <- rep(1, nrow(dataset))
        stdev <- "stdev2"
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
          if (is.null(dataset$LON) | is.null(dataset$LAT) ) stop("Missing coordinates on sampling frame")
        }
        if (model$type[i] == "linear" & is.na(model$gamma[i])) stop("gamma for Y variable ",i,"must be specified")
      }
    }
    #---------------------------------------------------------     
#    numdom <- length(levels(as.factor(dataset$DOMAINVALUE)))
    strata <- NULL
    dataset$DOMAINVALUE <- as.numeric(dataset$DOMAINVALUE)
    domain <- dataset[dataset$DOMAINVALUE == dominio, ]
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
    domain$STRATO <- as.factor(vett)
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
            # heteroscedastic residual variance sig2 * Y^(2*gamma), consistent with
            # computeGamma() (sd of residuals proportional to x^gamma) and with
            # the C++ implementation in buildStrataDF()
            st <- paste("gammas <- tapply(Y^(2*model$gamma[",i,"]),STRATO,sum) / as.numeric(table(STRATO))",sep="")
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
            psill <- model$sig2[i]
            range <- model$range[i]
            sd2 <- sapply(l.split, function(df) stdev4(df,psill,range))
            # stdev <- "stdev4"
            # stmt2 <- paste("sd2 <- sapply(l.split, function(df) ",
            #               stdev, "(df,psill,range))",
            #               sep = "")
            # eval(parse(text=stmt2))
            #-- PART III ---------------
            stmt <- paste("cov1 <- sapply(l.split, function(df,y,w) ",
                           "cov1(df,psill,range,df[,y],df[,w],beta1,beta2), y = 'Y",i,"',w = 'W",i,"')",
                           sep = "")
            eval(parse(text=stmt))
            #-- TOTAL S ---------------
            st <- paste("S",i," <- sqrt(sd1^2 + sd2^2 + cov1^2)",sep="")
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
    DOM1 <- rep(as.character(dominio), length(levels(domain$STRATO)))
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
    colnames(strata) <- toupper(colnames(strata))
    strata$DOM1 <- as.factor(strata$DOM1)
    indx <- sapply(strata, is.factor)
    strata[indx] <- lapply(strata[indx], function(x) as.numeric(as.character(x)))
      for (j in (1:nvarX)) {
        stmt <- paste("strata$X",j," <- as.numeric(strata$X",j,")",sep="")
        eval(parse(text=stmt))
      }
      # replace null means with a tiny positive value, for every Y variable
      for (i in 1:nvarY) {
        mcol <- paste0("M", i)
        strata[[mcol]][strata[[mcol]] == 0] <- 0.000000000000001
      }
    # cat("\nNumber of strata: ",nrow(strata))
    # cat("\n... of which with only one unit: ",sum(strata$N==1))
    return(strata)
}
