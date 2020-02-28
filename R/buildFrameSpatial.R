buildFrameSpatial <- function(df,id,X,Y,variance,lon,lat,domainvalue) {
  vars <- colnames(df)
  if (!id %in% vars) stop("Id name not in the frame variables")
  if (!lon %in% vars) stop("lon name not in the frame variables")
  if (!lat %in% vars) stop("lat name not in the frame variables")
  
  for (i in (1:length(X))) {
    var = X[i]
    if (!var %in% vars) {
      msg = paste("X variable ",var," not in the frame variables",sep="")
      stop(msg)
    }
  }
  for (i in (1:length(Y))) {
    var = Y[i]
    if (!var %in% vars) {
       msg = paste("Y variable ",var," not in the frame variables",sep="")
       stop(msg)
    }
  }
  if (length(variance) != length(Y)) stop("Number of variances not equal to number of Y's")
  for (i in (1:length(variance))) {
    var = variance[i]
    if (!var %in% vars) {
      msg = paste("Variance variable ",var," not in the frame variables",sep="")
      stop(msg)
    }
  }
  if (!domainvalue %in% vars) stop("Domain name not in the frame variables")  
  dframe <- NULL
  stmt <- paste("dframe$id <- df$",as.character(id),sep="")
  eval(parse(text=stmt))
  for (i in (1:length(X))) {
    aux <- paste("dframe$X",i," <- df$",as.character(X[i]),sep="")
    eval(parse(text=aux))
  }
  for (i in (1:length(Y))) {
    target <- paste("dframe$Y",i," <- df$",as.character(Y[i]),sep="")
    eval(parse(text=target))
  }
  for (i in (1:length(variance))) {
    target <- paste("dframe$var",i," <- df$",as.character(variance[i]),sep="")
    eval(parse(text=target))
  }
  target <- paste("dframe$lon <- df$",as.character(lon),sep="")
  eval(parse(text=target))
  target <- paste("dframe$lat <- df$",as.character(lat),sep="")
  eval(parse(text=target))
  stmt <- paste("dframe$domainvalue <- df$",as.character(domainvalue),sep="")
  eval(parse(text=stmt))
  dframe <- as.data.frame(dframe,stringsAsFactors = TRUE)
  dframe$domainvalue <- as.integer(dframe$domainvalue)
  return(dframe)
}