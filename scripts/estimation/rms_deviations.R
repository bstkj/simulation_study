# Reproduce "RMS deviations for MFE and ML estimates for the Müller network" figures
# in main text (Fig 4.3)

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

# rmse for maximum factored energy (MFE) estimates of root mean, mu
rmse_fe_mu_p1 <- rep(0,4) # good partition (|LL-FE|/|LL| <= 0.2)
rmse_fe_mu_p2 <- rep(0,4) # bad partition (|LL-FE|/|LL| > 0.2)
rmse_fe_mu_all <- rep(0,4) # overall
# rmse for ML estimates of mu
rmse_ll_mu_all <- rep(0,4) # overall
# rms deviation between MFE and ML estimates of mu
rmse_mu_p1 <- rep(0,4) # good partition
rmse_mu_p2 <- rep(0,4) # bad partition
rmse_mu_all <- rep(0,4) # overall
for (i in 1:4) {
  v1 <- vapply(optfe_muller_mv$Θopt[grp1], function(x) x$μ[i], FUN.VALUE=0.0) # good
  v2 <- vapply(optfe_muller_mv$Θopt[grp2], function(x) x$μ[i], FUN.VALUE=0.0) # bad
  v12 <- vapply(optfe_muller_mv$Θopt, function(x) x$μ[i], FUN.VALUE=0.0) # all
  rmse_fe_mu_p1[i] <- sqrt(mean(v1^2))
  rmse_fe_mu_p2[i] <- sqrt(mean(v2^2))
  rmse_fe_mu_all[i] <- sqrt(mean(v12^2))
  v3 <- vapply(optll_muller_mv$Θopt[grp1], function(x) x$μ[i], FUN.VALUE=0.0)
  v4 <- vapply(optll_muller_mv$Θopt[grp2], function(x) x$μ[i], FUN.VALUE=0.0)
  v34 <- vapply(optll_muller_mv$Θopt, function(x) x$μ[i], FUN.VALUE=0.0)
  rmse_ll_mu_all[i] <- sqrt(mean(v34^2))
  rmse_mu_p1[i] <- sqrt(mean(v1-v3)^2)
  rmse_mu_p2[i] <- sqrt(mean(v2-v4)^2)
  rmse_mu_all[i] <- sqrt(mean(v12-v34)^2)
}

## Fig 4.3, bottom left
png("figures/opt_muller_mv_mu_rmse.png", width=687, height=690, res=120)
  plot(1, type="n", xlim=c(0,4), ylim=c(0,40), yaxt="n", ylab="", xlab="", xaxt="n")
  axis(1, at=0.5+0:3, tick=F, line=0.5,
       labels=expression(widehat(mu)[1], widehat(mu)[2], widehat(mu)[3], widehat(mu)[4]))
  axis(2, at=seq(0,40,10), labels=seq(0,40,10), las=2)
  mtext(text="RMS deviation", line=2.5, side=2)
  for (i in 1:4) {
    # good partition
    lines(x=c(0.3,0.5)+(i-1), pch=c(16,1), type="p", col=rgb(0,0,1,0.5),
          y=c(rmse_mu_p1[i], rmse_fe_mu_p1[i]))
    # bad partition
    points(x=c(0.3,0.5)+(i-1), pch=c(16,1),
           y=c(rmse_mu_p2[i], rmse_fe_mu_p2[i]), col=rgb(1,0,0,0.5))
    # all
    lines(x=c(0.3,0.5,0.7)+(i-1), type="p", col=c(rgb(0,0,0,0.5), rgb(0,0,0,0.5), "black"),
          y=c(rmse_mu_all[i], rmse_fe_mu_all[i], rmse_ll_mu_all[i]), pch=c(16,1,10))
    # vertical lines to create visual separation
    lines(x=c(0.5,0.5)+(i-1), type="l", col=rgb(0,0,0,0.5),
          y=c(rmse_fe_mu_p1[i], rmse_fe_mu_p2[i]))
  }
dev.off()

################################################################################
### Variance

# rmse for MFE estimates of trait variances, sigma
rmse_fe_R_p1 <- rep(0,4)
rmse_fe_R_p2 <- rep(0,4)
rmse_fe_R_all <- rep(0,4)
# rmse for ML estimates of sigma
rmse_ll_R_all <- rep(0,4)
# rms deviation between MFE and ML estimates of sigma
rmse_R_p1 <- rep(0,4)
rmse_R_p2 <- rep(0,4)
rmse_R_all <- rep(0,4)

for (i in 1:4) {
  v1 <- vapply(optfe_muller_mv$Θopt[grp1], 
               function(x) log10(diag(x$R)[i]/diagR_true[i]), FUN.VALUE=0.0)
  v2 <- vapply(optfe_muller_mv$Θopt[grp2], 
               function(x) log10(diag(x$R)[i]/diagR_true[i]), FUN.VALUE=0.0)
  v12 <- vapply(optfe_muller_mv$Θopt, 
                function(x) log10(diag(x$R)[i]/diagR_true[i]), FUN.VALUE=0.0)
  rmse_fe_R_p1[i] <- sqrt(mean(v1^2))
  rmse_fe_R_p2[i] <- sqrt(mean(v2^2))
  rmse_fe_R_all[i] <- sqrt(mean(v12^2))
  v3 <- vapply(optll_muller_mv$Θopt[grp1],
               function(x) log10(diag(x$R)[i]/diagR_true[i]), FUN.VALUE=0.0)
  v4 <- vapply(optll_muller_mv$Θopt[grp2],
               function(x) log10(diag(x$R)[i]/diagR_true[i]), FUN.VALUE=0.0)
  v34 <- vapply(optll_muller_mv$Θopt,
                function(x) log10(diag(x$R)[i]/diagR_true[i]), FUN.VALUE=0.0)
  rmse_ll_R_all[i] <- sqrt(mean(v34^2))
  rmse_R_p1[i] <- sqrt(mean((v1-v3)^2))
  rmse_R_p2[i] <- sqrt(mean((v2-v4)^2))
  rmse_R_all[i] <- sqrt(mean((v12-v34)^2))
}

## Fig 4.3, bottom left
png("figures/opt_muller_mv_rate_rmse.png", width=687, height=690, res=120)
  plot(1, type="n", xlim=c(0,4), ylim=c(0,5), yaxt="n", ylab="", xlab="", xaxt="n")
  axis(1, at=0.5+0:3, tick=F, line=0.5,
       labels=expression("log"*widehat(Sigma)["1,1"], "log"*widehat(Sigma)["2,2"],
                         "log"*widehat(Sigma)["3,3"], "log"*widehat(Sigma)["4,4"]))
  axis(2, at=0:5, labels=0:5, las=2)
  mtext(text="RMS deviation", line=2, side=2)
  for (i in 1:4) {
    # good partition
    lines(x=c(0.3,0.5)+(i-1), type="p", pch=c(16,1), col="blue",
          y=c(rmse_R_p1[i], rmse_fe_R_p1[i]))
    # bad partition
    lines(x=c(0.3,0.5)+(i-1), type="p", pch=c(16,1), col="red",
          y=c(rmse_R_p2[i], rmse_fe_R_p2[i]))
    # all
    lines(x=c(0.3,0.5,0.7)+(i-1), type="p", pch=c(16,1,10),
          y=c(rmse_R_all[i], rmse_fe_R_all[i], rmse_ll_R_all[i]))
    # vertical lines to create visual separation
    lines(x=c(0.3,0.3)+(i-1), type="l", y=c(rmse_R_p1[i], rmse_R_p2[i]), col=rgb(0,0,0,0.5))
    lines(x=c(0.5,0.5)+(i-1), type="l", y=c(rmse_fe_R_p1[i], rmse_fe_R_p2[i]), col=rgb(0,0,0,0.5))
  }
dev.off()

################################################################################
### Legend

## Fig 4.3, bottom right
png("figures/opt_muller_mv_legend.png", width=687, height=690, res=120)
  plot(x=0, type="n", xlim=c(0,4), ylim=c(0,5), yaxt="n", ylab="", xlab="", xaxt="n", frame.plot=F)
  legend(bty="n", fill=c("blue", "red", "black"), border=NA, x.intersp=0.6, title.adj=0.5,
         x=0.76, y=5, title=expression(underline("Partition of datasets")),
         legend=expression(textstyle(abs(frac(widehat(MFE)-MLL,MLL)) <= 0.05),
                           textstyle(abs(frac(widehat(MFE)-MLL,MLL)) > 0.27), "all"), cex=1.2, box.lwd=2)
  
  legend(x=0.86, y=2, bty="n", title=expression(underline("Deviation")), title.adj=0.5,
         legend=expression(widehat(theta)["MFE"]-widehat(theta)["ML"],
                           widehat(theta)["MFE"]-theta["true"],
                           widehat(theta)["ML"]-theta["true"]),
         pch=c(16,1,10), pt.cex=2, cex=1.2)
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

Rtrue <- matrix(c(0.8, -0.71, -0.8, 0.49,
                  -0.71, 0.8, 0.81, -0.41,
                  -0.8, 0.81, 1.1, -0.4,
                  0.49, -0.41, -0.4, 0.5), nrow=4, ncol=4)
uppertri_true <- cov2cor(Rtrue)[upper.tri(diag(4))]
rho_fe <- list(rho12, rho13, rho23, rho14, rho24, rho34)
rho_ll <- list(rho12_mle, rho13_mle, rho23_mle, rho14_mle, rho24_mle, rho34_mle)

# rmse for MFE estimates of trait correlations, rho
rmse_fe_rho_p1 <- rep(0,6)
rmse_fe_rho_p2 <- rep(0,6)
rmse_fe_rho_all <- rep(0,6)
# rmse for ML estimates of rho
rmse_ll_rho_all <- rep(0,6)
# rms deviation for MFE and ML estimates of rho
rmse_rho_p1 <- rep(0,6)
rmse_rho_p2 <- rep(0,6)
rmse_rho_all <- rep(0,6)
for (i in 1:6) {
  v1 <- (rho_fe[[i]]-uppertri_true[i])[grp1]
  v2 <- (rho_fe[[i]]-uppertri_true[i])[grp2]
  v12 <- rho_fe[[i]]-uppertri_true[i]
  rmse_fe_rho_p1[i] <- sqrt(mean(v1^2))
  rmse_fe_rho_p2[i] <- sqrt(mean(v2^2))
  rmse_fe_rho_all[i] <- sqrt(mean(v12^2))
  v3 <- (rho_ll[[i]]-uppertri_true[i])[grp1]
  v4 <- (rho_ll[[i]]-uppertri_true[i])[grp2]
  v34 <- rho_ll[[i]]-uppertri_true[i]
  rmse_ll_rho_all[i] <- sqrt(mean(v34^2))
  
  rmse_rho_p1[i] <- sqrt(mean((v1-v3)^2))
  rmse_rho_p2[i] <- sqrt(mean((v2-v4)^2))
  rmse_rho_all[i] <- sqrt(mean(v12-v34)^2)
}

## Fig 4.3, top right
png("figures/opt_muller_mv_cor_rmse.png", width=687, height=690, res=120)
  plot(1, type="n", xlim=c(0,6), ylim=c(0,1.5), yaxt="n", ylab="", xlab="", xaxt="n")
  axis(1, at=0.6+0:5, tick=F, line=0.5,
       labels=expression(widehat(rho)["1,2"], widehat(rho)["1,3"],
                         widehat(rho)["2,3"], widehat(rho)["1,4"],
                         widehat(rho)["2,4"], widehat(rho)["3,4"]))
  axis(2, at=seq(0,1.5,0.5), labels=seq(0,1.5,0.5), las=2)
  mtext(text="RMS deviation", side=2, line=2.5)
  for (i in 1:6) {
    # good partition
    lines(x=c(0.3,0.5)+(i-1), pch=c(16,1), col="blue", type="p",
          y=c(rmse_rho_p1[i], rmse_fe_rho_p1[i]))
    # bad partition
    lines(x=c(0.3,0.5)+(i-1), pch=c(16,1), col="red", type="p",
          y=c(rmse_rho_p2[i], rmse_fe_rho_p2[i]))
    # all
    lines(x=c(0.3,0.5,0.7)+(i-1), pch=c(16, 1, 10), type="p",
          y=c(rmse_rho_all[i], rmse_fe_rho_all[i], rmse_ll_rho_all[1]))
    # vertical lines to create visual separation
    lines(x=c(0.3,0.3)+(i-1), y=c(rmse_rho_p2[i],rmse_rho_p1[i]),
          col=rgb(0,0,0,0.5))
    lines(x=c(0.5,0.5)+(i-1), y=c(rmse_fe_rho_p2[i],rmse_fe_rho_p1[i]),
          col=rgb(0,0,0,0.5))
  }
dev.off()