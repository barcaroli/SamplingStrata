plotStrata2d <- function (x, 
                          outstrata,
                          domain, 
                          vars, 
                          labels = NULL) 
{ 
  x1.max <- x1.min <- x1_max <- x1_min <- x2.max <- x2.min <- x2_max <- x2_min <- NULL
  colnames(x) <- toupper(colnames(x))
  if (!domain %in% levels(as.factor(x$DOMAINVALUE)))
    stop("Domain out of bounds")
  if (length(vars) != 2) stop("Indicate just two variables...")
  if (is.null(labels)) labels=vars
  outstrata <- outstrata[outstrata$DOM1 == domain,]
  outstrata <- outstrata[order(as.numeric(outstrata$STRATO)),]
  out <- NULL
  out$Stratum <- outstrata$STRATO
  out$Population <- outstrata$N
  out$Allocation <- round(outstrata$SOLUZ)
  out$'Sampling rate' <- outstrata$SOLUZ / outstrata$N
  
  x <- x[x$DOMAINVALUE == domain,]
  nstrata <- length(levels(as.factor(x$LABEL)))
  stringa <- paste("x1_min <- tapply(x$",vars[1],",x$LABEL,min)",sep="")
  eval(parse(text=stringa))
  stringa <- paste("x2_min <- tapply(x$",vars[2],",x$LABEL,min)",sep="")
  eval(parse(text=stringa))
  stringa <- paste("x1_max <- tapply(x$",vars[1],",x$LABEL,max)",sep="")
  eval(parse(text=stringa))
  stringa <- paste("x2_max <- tapply(x$",vars[2],",x$LABEL,max)",sep="")
  eval(parse(text=stringa))
  xcuts <- c(c(x1_min[2:(length(x1_min)-1)]),x1_max[length(x1_min)],x1_max[length(x1_max)])
  # out$bounds_X1 <- paste(c(x1_min[1:length(x1_min)-1],x1_min[length(x1_max)]),
  #                        c(c(x1_min[2:(length(x1_min)-1)]),x1_max[length(x1_min)],x1_max[length(x1_max)])
  #                        ,sep="-")
  # out$bounds_X2 <- paste(c(x2_min[1:length(x2_min)-1],x2_min[length(x2_max)]),
  #                        c(c(x2_min[2:(length(x2_min)-1)]),x2_max[length(x2_min)],x2_max[length(x2_max)])
  #                        ,sep="-")
  out$bounds_X1 <- paste(x1_min,x1_max,sep="-")
  out$bounds_X2 <- paste(x2_min,x2_max,sep="-")
  out <- as.data.frame(out,stringsAsFactors = TRUE) 
  lab1 <- paste("Bounds",labels[1])
  lab2 <- paste("Bounds",labels[2])
  colnames(out) <- c("Stratum","Population",
                     "Allocation","SamplingRate",
                     lab1,
                     lab2)
  
  cuts <- list(x1_max,x2_max)
  m <- length(cuts[[1]])

  stringa <- paste("x1.min <- round(min(x$",vars[1],", na.rm = TRUE),2)",sep="")
  eval(parse(text=stringa)) 
  stringa <- paste("x1.max <- round(max(x$",vars[1],", na.rm = TRUE),2)",sep="")
  eval(parse(text=stringa)) 
  stringa <- paste("x2.min <- round(min(x$",vars[2],", na.rm = TRUE),2)",sep="")
  eval(parse(text=stringa)) 
  stringa <- paste("x2.max <- round(max(x$",vars[2],", na.rm = TRUE),2)",sep="")
  eval(parse(text=stringa)) 
  cols <- rainbow(m + 1, alpha = 0.3)
  xcuts <- round(cuts[[1]],2)
  ycuts <- round(cuts[[2]],2)
  xcuts <- c(ifelse(x1.min > 0, 0, x1.min), xcuts, x1.max)
  ycuts <- c(ifelse(x2.min > 0, 0, x1.min), ycuts, x2.max)
  id <- c()
  value <- c()
  xs <- c()
  ys <- c()
  for (i in 1:(m + 1)) {
    if (i == 1) {
      xs <- c(xs, xcuts[i], xcuts[i], xcuts[i + 1], xcuts[i + 
                                                            1])
      ys <- c(ys, ycuts[i], ycuts[i + 1], ycuts[i + 1], 
              ycuts[i])
      id <- c(id, rep(i, 4))
      value <- c(value, rep(i, 4))
    }
    else {
      xs <- c(xs, xcuts[1], xcuts[1], xcuts[i + 1], xcuts[i + 
                                                            1], xcuts[i], xcuts[i])
      ys <- c(ys, ycuts[i], ycuts[i + 1], ycuts[i + 1], 
              ycuts[1], ycuts[1], ycuts[i])
      id <- c(id, rep(i, 6))
      value <- c(value, rep(i, 6))
    }
  }
  poly <- data.frame(id = as.factor(id), value = as.factor(value), 
                     x = xs, y = ys, stringsAsFactors = TRUE)

  stringa <- paste("plot(x$",vars[1],",x$",vars[2],",type='n',cex=0.01,xlab=labels[1],ylab=labels[2])",sep="")
  eval(parse(text=stringa)) 
  # plot(x$1,x$2,type="n",cex=0.01,xlab=labels[1],ylab=labels[2])
  cl <- c("orange","coral","sandybrown","cyan2","darkturquoise",
          "green","yellow","orange","coral","sandybrown","cyan2","darkturquoise",
          "green","yellow","orange","coral","sandybrown","cyan2","darkturquoise",
          "green","yellow")
  # cl <- gray(c(1:(nstrata+1)/(nstrata+1),alpha=NULL))
  for (i in (1:nstrata)) {
    j = i - 1
    m <- j - length(cl) * floor(j/length(cl)) + 1
    eval(parse(text=paste("polycorr <- poly[poly$value==",i,",]",sep="")))
    eval(parse(text=paste("polygon(polycorr$x,polycorr$y,col=cl[",m,"])",sep="")))
  }
  legend("topright", 
         # inset=c(-0.2,0),
         title="Strata",
         legend = c(as.character(c(1:(nstrata)))), 
         col = cl,
         ncol = 1, cex = 0.7, lwd = 3, text.font = 1, 
         text.col ="black",
         box.lty=1)
  title(paste("Strata boundaries in domain ",domain,sep=""),
        font.main=1,
        # col.main="red",
        cex.main = 1)
  stringa <- paste("points(x$",vars[1],",x$",vars[2],",cex=0.4)",sep="")
  eval(parse(text=stringa))
  # points(x$X1,x$X2,cex=0.4)
  
  t <- formattable(out,
                   list(
                     area(col = 2) ~ color_tile("#DeF7E9", "#71CA97"), 
                     area(col = 3) ~ color_tile("#DeF7E9", "#71CA97"),
                     'SamplingRate' = color_bar("#FA614B")))
  
  
  return(t)
  
}
