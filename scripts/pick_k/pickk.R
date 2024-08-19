# Reproduce "relative deviation and runtime vs maximum cluster size" figures in
# main text (Figs 4.1 and 4.2)

################################################################################
### Read in data

## 100 reps/datasets per k per network
# colnames: k (max cluster size), resid (LL-FE), residrel ((LL-FE)/abs(LL)),
# runtimeb (1st runtime), runtimea (2nd runtime)
sikora_uv <- read.table("results/pick_k/sikora_uv.csv", header=T, sep="")
sikora_mv <- read.table("results/pick_k/sikora_mv.csv", header=T, sep="")
lipson_uv <- read.table("results/pick_k/lipson_uv.csv", header=T, sep="")
lipson_mv <- read.table("results/pick_k/lipson_mv.csv", header=T, sep="")
muller_uv <- read.table("results/pick_k/muller_uv.csv", header=T, sep="")
muller_mv <- read.table("results/pick_k/muller_mv.csv", header=T, sep="")

## mean value over reps
# colnames: k, resid, residrel, absresidrel, runtimeb, runtimea
sikora_uv_mean <- read.table("results/pick_k/sikora_uv_mean.csv", header=T, sep="")
sikora_mv_mean <- read.table("results/pick_k/sikora_mv_mean.csv", header=T, sep="")
lipson_uv_mean <- read.table("results/pick_k/lipson_uv_mean.csv", header=T, sep="")
lipson_mv_mean <- read.table("results/pick_k/lipson_mv_mean.csv", header=T, sep="")
muller_uv_mean <- read.table("results/pick_k/muller_uv_mean.csv", header=T, sep="")
muller_mv_mean <- read.table("results/pick_k/muller_mv_mean.csv", header=T, sep="")

################################################################################
### Plot data

## Fig 4.1: Relative deviation and runtime vs maximum cluster size for the Sikora and Lipson networks
# top-left
png("figures/pickk_sikora_relresid.png", width=687, height=690, res=120)
  plot(x=sikora_uv$k-0.1, y=log10(sikora_uv$absresidrel), main="Sikora",
       xlab="k", ylab="abs((FE-LL)/LL)", ylim=c(-16, -2), xaxt="n", yaxt="n",
       xlim=c(2.8,5.2))
  axis(1, at=c(3,4,5), labels=c(3,4,5))
  axis(2, at=c(-4,-12,-15), labels=expression(10^-4,10^-12,10^-15),
       las=2, cex.axis=0.8)
  lines(x=sikora_uv_mean$k-0.1, y=log10(sikora_uv_mean$absresidrel), type="o")
  points(x=sikora_mv$k+0.1, y=log10(sikora_mv$absresidrel), pch=0)
  lines(x=sikora_mv_mean$k+0.1, y=log10(sikora_mv_mean$absresidrel), pch=0, type="o")
  legend(x=2.8, y=-12, legend=c("univariate", "multivariate"), pch=c(1,0), bty="n")
dev.off()
# top-right
png("figures/pickk_lipson_relresid.png", width=687, height=690, res=120)
  plot(x=lipson_uv$k-0.1, y=log10(lipson_uv$absresidrel), main="Lipson",
       xlab="k", ylab="abs((FE-LL)/LL)", ylim=c(-16, -2), yaxt="n", xlim=c(2.8, 7.2))
  axis(2, at=c(-3,-12,-15), labels=expression(10^-3,10^-12,10^-15), las=2,
       cex.axis=0.8)
  lines(x=lipson_uv_mean$k-0.1, y=log10(lipson_uv_mean$absresidrel), type="o")
  points(x=lipson_mv$k+0.1, y=log10(lipson_mv$absresidrel), pch=0)
  lines(x=lipson_mv_mean$k+0.1, y=log10(lipson_mv_mean$absresidrel), type="o", pch=0)
dev.off()
# bottom-left
png("figures/pickk_sikora_latency.png", width=687, height=690, res=120)
  plot(x=sikora_uv$k-0.1, y=log10(sikora_uv$runtimeb), main="Sikora",
       xlab="k", ylab="Runtime (s)", ylim=c(-2.7,1.1), xaxt="n", yaxt="n",
       xlim=c(2.8,5.2), col=rgb(1,0,0,0.5))
  axis(1, at=c(3,4,5), labels=c(3,4,5))
  axis(2, at=c(-2,-1,0,1), labels=expression(0.01,0.1,1,10), las=2)
  lines(x=sikora_uv_mean$k-0.1, y=log10(sikora_uv_mean$runtimeb), type="o",
        col="red")
  points(x=sikora_uv$k-0.1, y=log10(sikora_uv$runtimea), col=rgb(0,0,1,0.5))
  lines(x=sikora_uv_mean$k-0.1, y=log10(sikora_uv_mean$runtimea), type="o",
        col="blue")
  points(x=sikora_mv$k+0.1, y=log10(sikora_mv$runtimeb), pch=0, col=rgb(1,0,0,0.5))
  lines(x=sikora_mv_mean$k+0.1, y=log10(sikora_mv_mean$runtimeb), type="o", pch=0,
        col="red")
  points(x=sikora_mv$k+0.1, y=log10(sikora_mv$runtimea), pch=0, col=rgb(0,0,1,0.5))
  lines(x=sikora_mv_mean$k+0.1, y=log10(sikora_mv_mean$runtimea), type="o", pch=0,
        col="blue")
  legend(x=2.8, y=0.1, legend=c("1st runtime", "2nd runtime"),
         col=c("red", "blue"), seg.len=0.8, lty=c("solid", "solid"), bty="n")
dev.off()
# bottom-right
png("figures/pickk_lipson_latency.png", width=687, height=690, res=120)
  plot(x=lipson_uv$k-0.1, y=log10(lipson_uv$runtimeb), main="Lipson",
       xlab="k", ylab="Runtime (s)", ylim=c(-2.7,1.1), yaxt="n",
       col=rgb(1,0,0,0.5), xlim=c(2.8,7.2))
  axis(2, at=c(-2,-1,0,1), labels=expression(0.01,0.1,1,10), las=2)
  lines(x=lipson_uv_mean$k-0.1, y=log10(lipson_uv_mean$runtimeb), type="o",
        col="red")
  points(x=lipson_uv$k-0.1, y=log10(lipson_uv$runtimea), col=rgb(0,0,1,0.5))
  lines(x=lipson_uv_mean$k-0.1, y=log10(lipson_uv_mean$runtimea), type="o",
        col="blue")
  points(x=lipson_mv$k+0.1, y=log10(lipson_mv$runtimeb), pch=0, col=rgb(1,0,0,0.5))
  lines(x=lipson_mv_mean$k+0.1, y=log10(lipson_mv_mean$runtimeb), type="o", pch=0,
        col="red")
  points(x=lipson_mv$k+0.1, y=log10(lipson_mv$runtimea), col=rgb(0,0,1,0.5), pch=0)
  lines(x=lipson_mv_mean$k+0.1, y=log10(lipson_mv_mean$runtimea),
        col="blue", type="o", pch=0)
dev.off()

## Fig 4.2: Relative deviation and runtime vs maximum cluster size for the Müller network
# top-left
png("figures/pickk_muller_uv_relresid.png", width=687, height=690, res=120)
  plot(x=muller_uv$k, y=log10(muller_uv$absresidrel), main="Muller (univariate)",
       xlab="k", ylab="abs((FE-LL)/LL)", ylim=c(-14,7), cex=0.5, yaxt="n")
  lines(x=muller_uv_mean$k, y=log10(muller_uv_mean$absresidrel), col="red",
        type="o", pch=16, cex=0.8)
  axis(2, at=c(-12,-2,-1,0,5), labels=expression(10^-12,0.01,0.1,1,10^5), las=2,
       cex.axis=1)
  abline(h=c(-1,-2), lty="dashed")
  abline(v=c(11,25), lty="dashed")
dev.off()
# top-right
png("figures/pickk_muller_mv_relresid.png", width=687, height=690, res=120)
  plot(x=muller_mv$k, y=log10(muller_mv$absresidrel), main="Muller (multivariate)",
       xlab="k", ylab="abs((FE-LL)/LL)", ylim=c(-14,7), cex=0.5, yaxt="n")
  lines(x=muller_mv_mean$k, y=log10(muller_mv_mean$absresidrel), col="red",
        type="o", pch=16, cex=0.8)
  axis(2, at=c(-12,-2,-1,0,5), labels=expression(10^-12,0.01,0.1,1,10^5), las=2,
       cex.axis=1)
  abline(h=c(-1,-2), lty="dashed")
  abline(v=c(11,20), lty="dashed")
dev.off()
# bottom-left
png("figures/pickk_muller_uv_latency.png", width=687, height=690, res=120)
  plot(x=muller_uv$k, y=log10(muller_uv$runtimeb), main="Muller (univariate)",
       xlab="k", ylab="Runtime (s)", ylim=c(-0.3,1.1), yaxt="n",
       cex=0.5, col=rgb(1,0,0,0.3))
  axis(2, at=log10(c(0.5,1,2,5,10)), labels=c(0.5,1,2,5,10), las=2)
  abline(v=35, lty="dashed")
  lines(x=muller_uv_mean$k, y=log10(muller_uv_mean$runtimeb), type="o", cex=0.1)
  points(x=muller_uv$k, y=log10(muller_uv$runtimea), cex=0.5, col=rgb(0,0,1,0.3))
  lines(x=muller_uv_mean$k, y=log10(muller_uv_mean$runtimea), type="o", cex=0.2)
dev.off()
# bottom-right
png("figures/pickk_muller_mv_latency.png", width=687, height=690, res=120)
  plot(x=muller_mv$k, y=log10(muller_mv$runtimeb), main="Muller (multivariate)",
       xlab="k", ylab="Runtime (s)", ylim=c(0.5,2.7), yaxt="n",
       cex=0.5, col=rgb(1,0,0,0.3))
  axis(2, at=c(log10(5),1,log10(20),log10(50),2,log10(200),log10(500)),
       labels=expression(5,10,20,50,100,200,500), las=2)
  abline(v=c(15,35), lty="dashed")
  lines(x=muller_mv_mean$k, y=log10(muller_mv_mean$runtimeb), type="o", cex=0.1)
  points(x=muller_mv$k, y=log10(muller_mv$runtimea), cex=0.5, col=rgb(0,0,1,0.3))
  lines(x=muller_mv_mean$k, y=log10(muller_mv_mean$runtimea), type="o", cex=0.2)
dev.off()