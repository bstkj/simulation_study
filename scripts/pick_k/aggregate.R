# Process results in:
  # pick_k/sikora_uv/, pick_k/sikora_mv/
  # pick_k/lipson_uv/, pick_k/lipson_mv/
  # pick_k/muller_uv/, pick_k/muller_mv/
# Create:
  # sikora_uv.csv, sikora_uv_mean.csv, sikora_mv.csv, sikora_mv_mean.csv
  # lipson_uv.csv, lipson_uv_mean.csv, lipson_mv.csv, lipson_mv_mean.csv
  # muller_uv.csv, muller_uv_mean.csv, muller_mv.csv, muller_mv_mean.csv

################################################################################
sikora_ks <- 3:5
lipson_ks <- 3:7
muller_ks <- c(3:20, seq(25, 50, 5), 51:54)

for (netname in c("sikora", "lipson", "muller")) {
  ks <- eval(parse(text=paste0(netname, "_ks")))
  numks <- length(ks)
  numrows <- 100*numks
  absresidrel <- rep(0, numrows) # abs(LL-FE)/abs(LL)
  absresidrelmean <- rep(0, numks)
  runtimeb <- rep(0, numrows) # 1st run
  runtimebmean <- rep(0, numks)
  runtimea <- rep(0, numrows) # 2nd run
  runtimeamean <- rep(0, numks)
  
  for (p in c("uv", "mv")) { # trait dimension
    idx <- 1
    for (j in 1:numks) {
      k <- ks[j]
      for (i in 1:100) {
        fp <- paste0("results/pick_k/", netname, "_", p, "/k", k,"/rep", i, ".csv")
        datarow <- read.table(fp, header=T)
        print(fp)
        absresidrel[idx] <- abs(datarow$residtrue / datarow$lltrue)
        runtimeb[idx] <- datarow$fetimeb
        runtimea[idx] <- datarow$fetimea
        idx <- idx + 1
      }
      absresidrelmean[j] <- mean(absresidrel[(1+100*(j-1)):(j*100)])
      runtimebmean[j] <- mean(runtimeb[(1+100*(j-1)):(j*100)])
      runtimeamean[j] <- mean(runtimea[(1+100*(j-1)):(j*100)])
    }
    write.table(data.frame(k=rep(ks, each=100), absresidrel, runtimeb, runtimea),
                paste0("results/pick_k/", netname, "_", p, ".csv"))
    write.table(data.frame(k=ks, absresidrel=absresidrelmean,
                           runtimeb=runtimebmean, runtimea=runtimeamean),
                paste0("results/pick_k/", netname, "_", p, "_mean.csv"))
  }
}