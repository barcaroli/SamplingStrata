evalSolution <- function (frame, 
                          outstrata, 
                          nsampl = 100, 
                          cens = NULL, 
                          writeFiles = TRUE,
                          outputFolder = file.path(getwd(),"simulation"),
                          progress = TRUE) 
{
  if ( !requireNamespace("formattable", quietly = TRUE) ){
    install.packages("formattable")
  }
  if (writeFiles == TRUE) {
  	if(dir.exists(outputFolder)){
  	  warning("Folder ", outputFolder," exists and will be deleted.")
  	  unlink(outputFolder)
	  } 
    if(!dir.exists(outputFolder)) dir.create(outputFolder)
  }
  colnames(frame) <- toupper(colnames(frame))
  numY <- length(grep("Y", toupper(colnames(frame))))
  numdom <- length(levels(as.factor(frame$DOMAINVALUE)))
  param <- array(0, c(numY, numdom))
  for (i in 1:numY) {
    stmt <- paste("param[i,] <- aggregate(frame$Y", i, ",by=list(frame$DOMAINVALUE),FUN=mean)[,2]", 
                  sep = "")
    eval(parse(text = stmt))
  }
  estim <- array(0, c(numdom, nsampl, numY))
  differ <- array(0, c(numdom, nsampl, numY))
  # create progress bar
  if (progress == TRUE) pb <- txtProgressBar(min = 0, max = nsampl, style = 3)
  for (j in (1:nsampl)) {
    if (progress == TRUE) Sys.sleep(0.1)
    # update progress bar
    if (progress == TRUE) setTxtProgressBar(pb, j)
    samp <- selectSample(frame, outstrata, verbatim = FALSE)
    if (!is.null(cens)) 
      samp <- rbind(cens, samp)
    for (k in 1:numY) {
      stmt <- paste("estim[,j,k] <- aggregate(samp$Y", 
                    k, "*samp$WEIGHT,by=list(samp$DOMAINVALUE),FUN=sum)[,2] / aggregate(samp$WEIGHT,by=list(samp$DOMAINVALUE),FUN=sum)[,2]", 
                    sep = "")
      eval(parse(text = stmt))
    }
  }
  if (progress == TRUE) close(pb)
  for (j in (1:nsampl)) {
    for (i in (1:numdom)) {
      for (k in 1:numY) {
        differ[i, j, k] <- estim[i, j, k] - param[k, 
                                                  i]
      }
    }
  }
  cv <- array(0, c(numdom, numY))
#  bias <- array(0, c(numdom, numY))
  for (k in 1:numY) {
    for (i in (1:numdom)) {
      cv[i, k] <- sd(estim[i, , k])/mean(estim[i, , k])
#      bias[i, k] <- mean(differ[i, , k]/mean(estim[i, , k]))
    }
  }
  cv <- as.data.frame(cv,stringsAsFactors = TRUE)
  for (i in 1:numY) {
    stmt <- paste("colnames(cv)[", i, "] <- c('CV", i, "')", 
                  sep = "")
    eval(parse(text = stmt))
  }
  cv$dom <- paste("DOM", c(1:numdom), sep = "")
  if (writeFiles == TRUE) 
    write.table(cv, file.path(outputFolder, "expected_cv.csv"), sep = ",", row.names = FALSE, 
                col.names = TRUE, quote = FALSE)
  cv1 <- NULL
  cv1$domainvalue <- rep((1:numdom), numY)
  cv1$cv <- rep(0, (numY * numdom))
  cv1$val <- rep(0, (numY * numdom))
  cv1 <- as.data.frame(cv1,stringsAsFactors = TRUE)
  for (k in (1:numY)) {
    for (j in (1:numdom)) {
      cv1$cv[(k - 1) * numdom + j] <- k
      stmt <- paste("cv1$val[(k-1)*numdom + j] <- cv$CV", 
                    k, "[j]", sep = "")
      eval(parse(text = stmt))
    }
  }
  # if (writeFiles == TRUE) 
  #   write.table(cv1, "expected_cv1.csv", sep = ",", row.names = FALSE, 
  #               col.names = TRUE, quote = FALSE)
  # bias <- as.data.frame(bias)
  # for (i in 1:numY) {
  #   stmt <- paste("colnames(bias)[", i, "] <- c('RelBias", 
  #                 i, "')", sep = "")
  #   eval(parse(text = stmt))
  # }
  # bias$dom <- paste("DOM", c(1:numdom), sep = "")
  # if (writeFiles == TRUE) 
  #   write.table(bias, "expected_rel_bias.csv", sep = ",", 
  #               row.names = FALSE, col.names = TRUE, quote = FALSE)
  # bias1 <- NULL
  # bias1$domainvalue <- rep((1:numdom), numY)
  # bias1$bias <- rep(0, (numY * numdom))
  # bias1$val <- rep(0, (numY * numdom))
  # bias1 <- as.data.frame(bias1)
  # for (k in (1:numY)) {
  #   for (j in (1:numdom)) {
  #     bias1$bias[(k - 1) * numdom + j] <- k
  #     stmt <- paste("bias1$val[(k-1)*numdom + j] <- bias$RelBias", 
  #                   k, "[j]", sep = "")
  #     eval(parse(text = stmt))
  #   }
  # }
  # if (writeFiles == TRUE) 
  #   write.table(bias1, "expected_rel_bias1.csv", sep = ",", 
  #               row.names = FALSE, col.names = TRUE, quote = FALSE)
  diff <- NULL
  diff$dom <- rep(0, numdom * nsampl)
  diff$samp <- rep(0, numdom * nsampl)
  for (i in (1:numY)) {
    stmt <- paste("diff$diff", i, " <- rep(0,numdom*nsampl)", 
                  sep = "")
    eval(parse(text = stmt))
  }
  diff <- as.data.frame(diff,stringsAsFactors = TRUE)
  for (i in (1:numdom)) {
    for (j in (1:nsampl)) {
      diff$dom[(i - 1) * nsampl + j] <- i
      diff$samp[(i - 1) * nsampl + j] <- j
      for (k in (1:numY)) {
        stmt <- paste("diff$diff", k, "[(i-1)*nsampl+j] <- differ[i,j,k]", 
                      sep = "")
        eval(parse(text = stmt))
      }
    }
  }
  if (writeFiles == TRUE) 
    write.table(diff, file.path(outputFolder, "differences.csv"), sep = ",", row.names = FALSE, 
                col.names = TRUE, quote = FALSE)
  ############################################   
# New code for bias  
  bias <- NULL
  for (i in (3:ncol(diff))) {
    stmt <- paste("bias$y",i-2," <- tapply(diff$diff",i-2,",diff$dom,mean)",sep="")
    eval(parse(text=stmt))
  } 
  bias$dom <- paste("DOM", c(1:numdom), sep = "")
  bias <- as.data.frame(bias,stringsAsFactors = TRUE)
  Y <- aggregate(frame[,grep("Y",colnames(frame))],by=list(frame$DOMAINVALUE),mean)
  numY <- sum(grepl("Y",colnames(frame)))
  bias[,c(1:numY)] <- round(bias[,c(1:numY)]/Y[,c(2:(numY+1))],4)
  if (writeFiles == TRUE)
    write.table(bias, file.path(outputFolder, "expected_rel_bias.csv"), sep = ",",
                row.names = FALSE, col.names = TRUE, quote = FALSE)
  if (numdom > 1) {
    if (writeFiles == TRUE) 
      # pdf("cv.pdf", width = 7, height = 5)
      png(file.path(outputFolder,"cv.png"))
    boxplot(val ~ cv, data = cv1, col = "orange", main = "Distribution of CV's in the domains", 
            xlab = "Variables Y", ylab = "Value of CV")
    if (writeFiles == TRUE) 
      dev.off()
    if (writeFiles == TRUE) 
      # pdf("rel_bias.pdf", width = 7, height = 5)
      png(file.path(outputFolder,"rel_bias.png"))
    boxplot(bias[,-ncol(bias)], col = "orange", main = "Distribution of relative bias in the domains",
            xlab = "Variables Y", ylab = "Relative bias")
    # boxplot(bias, col = "orange", main = "Distribution of relative bias in the domains", 
    #         xlab = "Variables Y", ylab = "Relative bias")
    if (writeFiles == TRUE) 
      dev.off()
    
  }
  # if (writeFiles == TRUE) 
  # # pdf("differences.pdf", width = 14, height = 10)
  # png("differences.png")
  # k <- ceiling(numY/4)
  # for (j in 1:k) {
  #   split.screen(c(2, 2))
  #   for (i in 1:4) {
  #     if (i + 4 * (j - 1) <= numY) {
  #       stmt <- paste("screen(", i, ")", sep = "")
  #       eval(parse(text = stmt))
  #       stmt <- paste("boxplot(diff", i, "~dom,data=diff,ylab='Differences',xlab='Domain',col = 'orange')", 
  #                     sep = "")
  #       eval(parse(text = stmt))
  #       stmt <- paste("mtext(expression(Y", i, "), side=3, adj=0, cex=1.0, line=1)", 
  #                     sep = "")
  #       eval(parse(text = stmt))
  #     }
  #   }
  #   close.screen(all.screens = TRUE)
  # }
  # if (writeFiles == TRUE) 
  #   dev.off()
  # results <- list(coeff_var = cv1, bias = bias1)
  est <- matrix(NA,nrow=numdom*nsampl,ncol=numY)
  est <- as.data.frame(est,stringsAsFactors = TRUE) 
  colnames(est) <- c(paste("Y",c(1:numY),sep=""))
  est$dom <- rep(c(1:numdom),each=nsampl)
  for (i in (1:numdom)) {
    est[est$dom == i,c(1:(numY))] <- estim[i,,]
  }
  if (writeFiles == TRUE) {
    write.table(est,file.path(outputFolder,"estimates.csv"),sep=",",row.names=F,col.names=F)
  }

  cv[,c(1:numY)] <- round(cv[,c(1:numY)],4)
  results <- list(coeff_var = cv, rel_bias = bias, est = est)
  cv <- cbind(c(1:nrow(cv)),cv)
  colnames(cv) <- c("domain",paste("cv(Y",c(1:numY),")",sep=""))
  bias <- cbind(c(1:nrow(bias)),bias)
  colnames(bias) <- c("domain",paste("bias(Y",c(1:numY),")",sep=""))
  cv <- formattable::formattable(cv,list(area(col = 2:(numY+1)) ~ color_tile("#DeF7E9", "#71CA97")))
  bias <- formattable::formattable(bias,list(area(col = 2:(numY+1)) ~ color_tile("#DeF7E9", "#71CA97")))
  return(results)
}
