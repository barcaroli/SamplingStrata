plotSamprate <- function (solution, dom) {
  df <- as.data.frame(solution$aggr_strata[solution$aggr_strata$DOM1 == dom,])
#   plot(df$SOLUZ/df$N,
#        ylim = c(0,1),
#        ylab="sampling units / total units",
#        xlab="strata",
# 	   type="h")
  barplot(SOLUZ/N ~ STRATO,data=df, ylab="sampling units / total units",xlab="strata",)
  title(paste("Sampling rate per stratum in domain ",dom,sep=""))
}    