buildErrorsDF <- function(frame=NULL,strata=NULL) {
  # if (is.null(frame & is.null(strata))) stop("Indicate a frame or a strata dataframe")
  # if (!is.null(frame & !is.null(strata))) stop("Indicate only one between frame and strata dataframe")
  if (!is.null(frame)) {
    ndom <- length(unique(frame$DOMAINVALUE))
    nvars <- length(grep("Y",names(frame)))
    st <- "cv <- as.data.frame(list(DOM=rep('DOM1',ndom)"
    for (k in c(1:nvars)) {
      if (k < nvars) {
        st <- paste0(st,",CV",k,"=rep(NA,ndom)")
      }
      if (k == nvars) {
        st <- paste0(st,",CV",k,"=rep(NA,ndom),domainvalue=c(1:ndom)))")
      }
    }
    eval(parse(text=st))
  }
  if (!is.null(strata)) {
    ndom <- length(unique(strata$DOM1))
    m_cols <- grep("^M\\d+$", names(strata), value = TRUE)
    nvars <- length(m_cols)      
    st <- "cv <- as.data.frame(list(DOM=rep('DOM1',ndom)"
    for (k in c(1:nvars)) {
      if (k < nvars) {
        st <- paste0(st,",CV",k,"=rep(NA,ndom)")
      }
      if (k == nvars) {
        st <- paste0(st,",CV",k,"=rep(NA,ndom),domainvalue=c(1:ndom)))")
      }
    }
    eval(parse(text=st))
  }
  return(cv)
}


