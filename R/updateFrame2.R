updateFrame2 <- function(frame,solution,writeFiles=FALSE) 
{
  indices <- as.data.frame(solution$indices)
  colnames(indices) <- c("id","LABEL")
  indices$STRATO <- indices$LABEL
  framenew <- merge(frame,indices,by=c("id"))
  if (writeFiles == TRUE) {
    write.table(framenew, "framenew.txt", row.names = FALSE, 
                sep = "\t", quote = FALSE)
  }
  return(framenew)
}