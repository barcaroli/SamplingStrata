plotStrata2d <- function (x, domain, vars, labels = NULL) 
{ 
  if (domain < 1 | domain > length(levels(as.factor(x$domainvalue))))
    stop("Domain out of bounds")
  if (is.null(labels)) labels=vars
  x <- x[x$domainvalue == domain,]
  nstrata <- length(levels(as.factor(x$LABEL)))
  if (length(vars) != 2) stop("Indicate just two variables...")
  stringa <- paste("x1 <- tapply(x$",vars[1],",x$LABEL,summary)",sep="")
  eval(parse(text=stringa))
  stringa <- paste("x2 <- tapply(x$",vars[2],",x$LABEL,summary)",sep="")
  eval(parse(text=stringa))
  stringa1 <- paste("cuts_x1 <- as.numeric(c(x1[[1]][6]")
  stringa2 <- paste("cuts_x2 <- as.numeric(c(x2[[1]][6]")

  for (i in (2:nstrata)) {
    if (i == nstrata) {
      stringa1 <- paste(stringa1,"))",sep="")
      stringa2 <- paste(stringa2,"))",sep="")
    }
    if (i < nstrata) {
      stringa1 <- paste(stringa1,",x1[[",i,"]][6]",sep="")
      stringa2 <- paste(stringa2,",x2[[",i,"]][6]",sep="")
    }
  }
  eval(parse(text=stringa1))
  eval(parse(text=stringa2))
  cuts <- list(cuts_x1,cuts_x2)
  m <- length(cuts[[1]])
  x1.min <- min(x$X1, na.rm = TRUE)
  x1.max <- max(x$X1, na.rm = TRUE)
  x2.min <- min(x$X2, na.rm = TRUE)
  x2.max <- max(x$X2, na.rm = TRUE)
  cols <- rainbow(m + 1, alpha = 0.3)
  xcuts <- cuts[[1]]
  ycuts <- cuts[[2]]
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
                     x = xs, y = ys)
  # p <- ggplot(poly, aes(x = x, y = y)) + geom_polygon(aes(fill = value, 
  #                                                         group = id))
  # stringa <- paste("p <- p + geom_point(data = x, aes_string(x = x$",vars[1],", y = x$",vars[2],"), 
  #         colour = rgb(0, 0, 0, 0.3)) +
  #         guides(fill = guide_legend(title = 'Strata'))", sep="")
  # eval(parse(text=stringa))
  # p + labs(x = vars[1]) + labs(y = vars[2])
  
  plot(x$X1,x$X2,type="n",cex=0.01,
       xlab=labels[1],ylab=labels[2])
  cl <- c("yellow","red","salmon","green","orange")
  # cl <- gray(c(1:(nstrata+1)/(nstrata+1),alpha=NULL))
  for (i in (1:nstrata)) {
    j = i - 1
    m <- j - length(cl) * floor(j/length(cl)) + 1
    eval(parse(text=paste("polycorr <- poly[poly$value==",i,",]",sep="")))
    eval(parse(text=paste("polygon(polycorr$x,polycorr$y,col=cl[",m,"])",sep="")))
  }
  legend("topright", 
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
  points(x$X1,x$X2,cex=0.4)
  
  cat("\n--------------------------------")
  cat("\nBoundaries for the 1st variable")
  cat("\n",as.numeric(xcuts))
  cat("\nBoundaries for the 2nd variable")
  cat("\n",as.numeric(ycuts))
  cat("\n--------------------------------")
  boundaries <- list(x1_boundaries = as.numeric(xcuts),
              x2_boundaries = as.numeric(ycuts))
  return(boundaries)
}
