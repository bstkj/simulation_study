# Reproduce "Residuals for MFE and ML estimates for the Müller network" figures
# in main text (Fig 4.4)

################################################################################
### Read in data

optfe_muller_mv <- readRDS("results/estimation/muller_mv_fe.rds")
# names(optfe_muller_mv): "FEopt" "Θopt"
# names(optfe_muller_mv$Θopt[[1]]): "μ" "R"
optll_muller_mv <- readRDS("results/estimation/muller_mv_ml.rds")
# names(optll_muller_mv): "LLopt" "Θopt"  "Θhat"
# names(optll_muller_mv$Θopt[[1]]): "μ" "R"
relresid <- ((unlist(optll_muller_mv$LLopt) - unlist(optfe_muller_mv$FEopt))/
               unlist(optll_muller_mv$LLopt))

grp1 <- which(abs(relresid) <= 0.2) # good partition
grp2 <- which(abs(relresid) > 0.2) # bad partition
diagR_true <- c(0.8,0.8,1.1,0.5) # diagonal values of true rate matrix

################################################################################
### Mean

## Fig 4.4, bottom right
png("figures/opt_muller_mv_mu.png", width=687, height=690, res=120)
  boxplot(
    vapply(optfe_muller_mv$Θopt, function(x) x$μ[1], FUN.VALUE=0.0),
    vapply(optll_muller_mv$Θopt, function(x) x$μ[1], FUN.VALUE=0.0),
    vapply(optfe_muller_mv$Θopt, function(x) x$μ[2], FUN.VALUE=0.0),
    vapply(optll_muller_mv$Θopt, function(x) x$μ[2], FUN.VALUE=0.0),
    vapply(optfe_muller_mv$Θopt, function(x) x$μ[3], FUN.VALUE=0.0),
    vapply(optll_muller_mv$Θopt, function(x) x$μ[3], FUN.VALUE=0.0),
    vapply(optfe_muller_mv$Θopt, function(x) x$μ[4], FUN.VALUE=0.0),
    vapply(optll_muller_mv$Θopt, function(x) x$μ[4], FUN.VALUE=0.0),
    xaxt="n", yaxt="n", cex=0.5, outline=F, col=NA,
    main="Muller (multivariate)", ylab=""
  )
  axis(2, at=c(-100,-50,0,50,100), labels=c(-100,-50,0,50,100), las=2)
  axis(1, at=c(1.5, 3.5, 5.5, 7.5), tick=F, line=0.5,
       labels=expression(widehat(mu)[1], widehat(mu)[2], widehat(mu)[3], widehat(mu)[4]))
  axis(1, at=1:8, tick=F, labels=rep(c("FE", "LL"), 4), line=-0.9)
  mtext(text=expression(widehat(mu)[i]-mu[i]), side=2, line=2.5)
  for (i in 1:4) {
    # bad partition
    points(x=rep((2*i-1)+0.2,33), pch=1, col=rgb(0,0,1,0.8),
           y=vapply(optfe_muller_mv$Θopt[grp2], function(x) x$μ[i], FUN.VALUE=0.0), cex=0.5)
    points(x=rep(2*i+0.2,33), pch=1, col=rgb(0,0,1,0.8),
           y=vapply(optll_muller_mv$Θopt[grp2], function(x) x$μ[i], FUN.VALUE=0.0), cex=0.5)
    # good partition
    points(x=rep((2*i-1)-0.2,67), pch=1, col=rgb(1,0,0,0.8),
           y=vapply(optfe_muller_mv$Θopt[grp1], function(x) x$μ[i], FUN.VALUE=0.0), cex=0.5)
    points(x=rep(2*i-0.2,67), pch=1, col=rgb(1,0,0,0.8),
           y=vapply(optll_muller_mv$Θopt[grp1], function(x) x$μ[i], FUN.VALUE=0.0), cex=0.5)
    # means for each partition
    lines(x=c((2*i-1)-0.2,(2*i-1)+0.2),
          y=c(mean(vapply(optfe_muller_mv$Θopt[grp1], function(x) x$μ[i], FUN.VALUE=0.0)),
              mean(vapply(optfe_muller_mv$Θopt[grp2], function(x) x$μ[i], FUN.VALUE=0.0))),
          lwd=2, type="o", pch=0)
    lines(x=c(2*i-0.2,2*i+0.2),
          y=c(mean(vapply(optll_muller_mv$Θopt[grp1], function(x) x$μ[i], FUN.VALUE=0.0)),
              mean(vapply(optll_muller_mv$Θopt[grp2], function(x) x$μ[i], FUN.VALUE=0.0))),
          lwd=2, type="o", pch=0)
  }
dev.off()

################################################################################
### Variance

## Fig 4.4, top left
png("figures/opt_muller_mv_rate_zoomedout.png", width=687, height=690, res=120)
  boxplot(
    vapply(optfe_muller_mv$Θopt, function(x) log10(diag(x$R)[1]/diagR_true[1]), FUN.VALUE=0.0),
    vapply(optll_muller_mv$Θopt, function(x) log10(diag(x$R)[1]/diagR_true[1]), FUN.VALUE=0.0),
    vapply(optfe_muller_mv$Θopt, function(x) log10(diag(x$R)[2]/diagR_true[2]), FUN.VALUE=0.0),
    vapply(optll_muller_mv$Θopt, function(x) log10(diag(x$R)[2]/diagR_true[2]), FUN.VALUE=0.0),
    vapply(optfe_muller_mv$Θopt, function(x) log10(diag(x$R)[3]/diagR_true[3]), FUN.VALUE=0.0),
    vapply(optll_muller_mv$Θopt, function(x) log10(diag(x$R)[3]/diagR_true[3]), FUN.VALUE=0.0),
    vapply(optfe_muller_mv$Θopt, function(x) log10(diag(x$R)[4]/diagR_true[4]), FUN.VALUE=0.0),
    vapply(optll_muller_mv$Θopt, function(x) log10(diag(x$R)[4]/diagR_true[4]), FUN.VALUE=0.0),
    yaxt="n", xaxt="n", cex=0.5, outline=F, col=NA, main="Muller (multivariate)",
    ylab="", ylim=c(0,12)
  )
  axis(2, at=seq(0,12,2), labels=seq(0,12,2), las=2)
  mtext(text=expression(paste("log", widehat(Sigma)["i,i"]) - paste("log", Sigma["i,i"])), side=2, line=2)
  axis(1, at=c(1.65, 3.65, 5.65, 7.65), tick=F, line=0.5,
       labels=expression(widehat(Sigma)["1,1"], widehat(Sigma)["2,2"],
                         widehat(Sigma)["3,3"], widehat(Sigma)["4,4"]))
  axis(1, at=1:8, tick=F, labels=rep(c("FE", "LL"), 4), line=-0.9)
  for (i in 1:4) {
    # bad partition
    points(x=rep((2*i-1)+0.2,33), pch=1, col=rgb(0,0,1,0.8),
           y=vapply(optfe_muller_mv$Θopt[grp2], function(x) log10(diag(x$R)[i]/diagR_true[i]), FUN.VALUE=0.0), cex=0.5)
    points(x=rep(2*i+0.2,33), pch=1, col=rgb(0,0,1,0.8),
           y=vapply(optll_muller_mv$Θopt[grp2], function(x) log10(diag(x$R)[i]/diagR_true[i]), FUN.VALUE=0.0), cex=0.5)
    # good partition
    points(x=rep((2*i-1)-0.2,67), pch=1, col=rgb(1,0,0,0.8),
           y=vapply(optfe_muller_mv$Θopt[grp1], function(x) log10(diag(x$R)[i]/diagR_true[i]), FUN.VALUE=0.0), cex=0.5)
    points(x=rep(2*i-0.2,67), pch=1, col=rgb(1,0,0,0.8),
           y=vapply(optll_muller_mv$Θopt[grp1], function(x) log10(diag(x$R)[i]/diagR_true[i]), FUN.VALUE=0.0), cex=0.5)
    # means for each partition
    lines(x=c((2*i-1)-0.2,(2*i-1)+0.2), y=c(mean(vapply(optfe_muller_mv$Θopt[grp1], function(x) log10(diag(x$R)[i]/diagR_true[i]), FUN.VALUE=0.0)),
                                mean(vapply(optfe_muller_mv$Θopt[grp2], function(x) log10(diag(x$R)[i]/diagR_true[i]), FUN.VALUE=0.0))),
          lwd=2, type="o", pch=0)
  }
dev.off()

## Fig 4.4, top right
png("figures/opt_muller_mv_rate_zoomedin.png", width=687, height=690, res=120)
  boxplot(
    vapply(optfe_muller_mv$Θopt, function(x) log10(diag(x$R)[1]/diagR_true[1]), FUN.VALUE=0.0),
    vapply(optll_muller_mv$Θopt, function(x) log10(diag(x$R)[1]/diagR_true[1]), FUN.VALUE=0.0),
    vapply(optfe_muller_mv$Θopt, function(x) log10(diag(x$R)[2]/diagR_true[2]), FUN.VALUE=0.0),
    vapply(optll_muller_mv$Θopt, function(x) log10(diag(x$R)[2]/diagR_true[2]), FUN.VALUE=0.0),
    vapply(optfe_muller_mv$Θopt, function(x) log10(diag(x$R)[3]/diagR_true[3]), FUN.VALUE=0.0),
    vapply(optll_muller_mv$Θopt, function(x) log10(diag(x$R)[3]/diagR_true[3]), FUN.VALUE=0.0),
    vapply(optfe_muller_mv$Θopt, function(x) log10(diag(x$R)[4]/diagR_true[4]), FUN.VALUE=0.0),
    vapply(optll_muller_mv$Θopt, function(x) log10(diag(x$R)[4]/diagR_true[4]), FUN.VALUE=0.0),
    yaxt="n", xaxt="n", cex=0.5, outline=F, col=NA,
    main="Muller (multivariate)",
    ylab="", ylim=c(-0.3,0.3)
  )
  axis(2, at=c(-0.3,-0.2,-0.1,0,0.1,0.2), labels=c(-0.3,-0.2,-0.1,0,0.1,0.2), las=2)
  mtext(text=expression(paste("log", widehat(Sigma)["i,i"]) - paste("log", Sigma["i,i"])), side=2, line=2.5)
  axis(1, at=c(1.65, 3.65, 5.65, 7.65), tick=F, line=0.5,
       labels=expression(widehat(Sigma)["1,1"], widehat(Sigma)["2,2"],
                         widehat(Sigma)["3,3"], widehat(Sigma)["4,4"]))
  axis(1, at=1:8, tick=F, labels=rep(c("FE", "LL"), 4), line=-0.9)
  for (i in 1:4) {
    # bad partition
    points(x=rep((2*i-1)+0.2,33), pch=1, col=rgb(0,0,1,0.8),
           y=vapply(optfe_muller_mv$Θopt[grp2], function(x) log10(diag(x$R)[i]/diagR_true[i]), FUN.VALUE=0.0), cex=0.5)
    points(x=rep(2*i+0.2,33), pch=1, col=rgb(0,0,1,0.8),
           y=vapply(optll_muller_mv$Θopt[grp2], function(x) log10(diag(x$R)[i]/diagR_true[i]), FUN.VALUE=0.0), cex=0.5)
    # good partition
    points(x=rep((2*i-1)-0.2,67), pch=1, col=rgb(1,0,0,0.8),
           y=vapply(optfe_muller_mv$Θopt[grp1], function(x) log10(diag(x$R)[i]/diagR_true[i]), FUN.VALUE=0.0), cex=0.5)
    points(x=rep(2*i-0.2,67), pch=1, col=rgb(1,0,0,0.8),
           y=vapply(optll_muller_mv$Θopt[grp1], function(x) log10(diag(x$R)[i]/diagR_true[i]), FUN.VALUE=0.0), cex=0.5)
    # means for each partition
    lines(x=c(2*i-0.2,2*i+0.2), y=c(mean(vapply(optll_muller_mv$Θopt[grp1], function(x) log10(diag(x$R)[i]/diagR_true[i]), FUN.VALUE=0.0)),
                                mean(vapply(optll_muller_mv$Θopt[grp2], function(x) log10(diag(x$R)[i]/diagR_true[i]), FUN.VALUE=0.0))),
          lwd=2, type="o", pch=0)
  }
dev.off()

################################################################################
### Correlation

rho12 <- rep(0,100); rho12_mle <- rep(0,100)
rho13 <- rep(0,100); rho13_mle <- rep(0,100)
rho14 <- rep(0,100); rho14_mle <- rep(0,100)
rho23 <- rep(0,100); rho23_mle <- rep(0,100)
rho24 <- rep(0,100); rho24_mle <- rep(0,100)
rho34 <- rep(0,100); rho34_mle <- rep(0,100)
for (i in 1:100) {
  uppertri <- cov2cor(optfe_muller_mv$Θopt[[i]]$R)[upper.tri(diag(4))]
  uppertri_mle <- cov2cor(optll_muller_mv$Θopt[[i]]$R)[upper.tri(diag(4))]
  rho12[i] <- uppertri[1]; rho12_mle[i] <- uppertri_mle[1]
  rho13[i] <- uppertri[2]; rho13_mle[i] <- uppertri_mle[2]
  rho23[i] <- uppertri[3]; rho23_mle[i] <- uppertri_mle[3]
  rho14[i] <- uppertri[4]; rho14_mle[i] <- uppertri_mle[4]
  rho24[i] <- uppertri[5]; rho24_mle[i] <- uppertri_mle[5]
  rho34[i] <- uppertri[6]; rho34_mle[i] <- uppertri_mle[6]
}

rho_fe <- list(rho12, rho13, rho23, rho14, rho24, rho34)
rho_ll <- list(rho12_mle, rho13_mle, rho23_mle, rho14_mle, rho24_mle, rho34_mle)

Rtrue <- matrix(c(0.8, -0.71, -0.8, 0.49,
                 -0.71, 0.8, 0.81, -0.41,
                 -0.8, 0.81, 1.1, -0.4,
                 0.49, -0.41, -0.4, 0.5), nrow=4, ncol=4)
uppertri_true <- cov2cor(Rtrue)[upper.tri(diag(4))]
rho12_true <- uppertri_true[1]
rho13_true <- uppertri_true[2]
rho23_true <- uppertri_true[3]
rho14_true <- uppertri_true[4]
rho24_true <- uppertri_true[5]
rho34_true <- uppertri_true[6]

## Fig 4.4, bottom left
png("figures/opt_muller_mv_cor.png", width=687, height=690, res=120)
  boxplot(
    rho12-rho12_true, rho12_mle-rho12_true,
    rho13-rho13_true, rho13_mle-rho13_true,
    rho23-rho23_true, rho23_mle-rho23_true,
    rho14-rho14_true, rho14_mle-rho14_true,
    rho24-rho24_true, rho24_mle-rho24_true,
    rho34-rho34_true, rho34_mle-rho34_true,
    yaxt="n", xaxt="n", cex=0.5, outline=F, col=NA,
    main="Muller (multivariate)",
    ylab="", ylim=c(-2,2)
  )
  axis(2, at=c(-2,-1,0,1,2), labels=c(-2,-1,0,1,2), las=2)
  mtext(text=expression(widehat(rho)[ij] - rho[ij]), side=2, line=2)
  axis(1, at=c(1.75, 3.75, 5.75, 7.75, 9.75, 11.75), tick=F, line=0.5,
       labels=expression(widehat(rho)["1,2"], widehat(rho)["1,3"],
                         widehat(rho)["2,3"], widehat(rho)["1,4"],
                         widehat(rho)["2,4"], widehat(rho)["3,4"]))
  axis(1, at=1:12, tick=F, labels=rep(c("FE", "LL"), 6), line=-0.9)
  
  for (i in 1:6) {
    # bad partition
    points(x=rep((2*i-1)+0.2,33), pch=1, col=rgb(0,0,1,0.8),
           y=rho_fe[[i]][grp2]-uppertri_true[i], cex=0.5)
    points(x=rep(2*i+0.2,33), pch=1, col=rgb(0,0,1,0.8),
           y=rho_ll[[i]][grp2]-uppertri_true[i], cex=0.5)
    # good partition
    points(x=rep((2*i-1)-0.2,67), pch=1, col=rgb(1,0,0,0.8),
           y=rho_fe[[i]][grp1]-uppertri_true[i], cex=0.5)
    points(x=rep(2*i-0.2,67), pch=1, col=rgb(1,0,0,0.8),
           y=rho_ll[[i]][grp1]-uppertri_true[i], cex=0.5)
    # means for each partition
    lines(x=c((2*i-1)-0.2,(2*i-1)+0.2),
          y=c(mean(rho_fe[[i]][grp1]), mean(rho_fe[[i]][grp2]))-uppertri_true[i],
          lwd=2, type="o", pch=0)
    lines(x=c(2*i-0.2,2*i+0.2),
          y=c(mean(rho_ll[[i]][grp1]), mean(rho_ll[[i]][grp2]))-uppertri_true[i],
          lwd=2, type="o", pch=0)
  }
dev.off()