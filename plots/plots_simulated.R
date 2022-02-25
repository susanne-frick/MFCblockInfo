####----------- example plots for paper, with simulated questionnaire ------####

library(MFCblockInfo)
# devtools::load_all()

# design.load.all <- readRDS("simulation/design_load_all_234.rds")
design.load.all <- readRDS("../simulation/design_load_all_234.rds")

nb <- 3

design.load <- design.load.all[[as.character(nb)]][["12"]]

K <- nrow(design.load)/nb

# trait correlations (Big 5 from meta-analysis van der Linden et al.)
trait.cov <- matrix(c(1,-.36,-.17,-.36,-.43,
                      -.36,1,.43,.26,.29,
                      -.17,.43,1,.21,.20,
                      -.36,.26,.21,1,.43,
                      -.43,.29,.20,.43,1),
                    nrow=5,ncol=5)

####------------------ simulate item parameters -------------------####
set.seed(1204+1)
items <- sim.items(design.load=design.load, K=K, nb=nb,
                   load.range=c(.65, .95),
                   int.range=c(-1,1))
load.mat <- items$loads * design.load
gamma.true <- create.design.mat(K=K, nb=nb) %*% items$u.mean

####----------------------- info for plots -------------------------####
# blockinfo
# info.block.plot <- calc.info.plot(tr.levels=seq(-2,2,.5), fix.levels=c(-1,0,1), fix.level.others=0,
#                                   K=K, which.blocks = 1:K,
#                                   FUN=lhb.mplus,
#                                   int=gamma.true, loads=load.mat, uni=diag(items$uni), nb=nb)
# saveRDS(info.block.plot, "info_block_plot.rds")

# testinfo
# info.test.plot <- calc.testinfo.plot(tr.levels=seq(-2,2,.5), fix.level.others=0,
#                                      K=K, which.blocks = 1:K,
#                                      FUN=lhb.mplus,
#                                      int=gamma.true, loads=load.mat, uni=diag(items$uni))
# saveRDS(info.test.plot, "info_test_plot.rds")

#only 1 trait varied, others randomly drawn from MVN, n for each trait level
# info.test.plot.1d <- calc.testinfo.1d(tr.levels=seq(-2,2,.5), sigma=trait.cov, n=100, which.blocks = 1:K,
#                                       seed=42,
#                                       FUN=lhb.mplus,
#                                       int=gamma.true, loads=load.mat, uni=diag(items$uni),
#                                       K=K, nb=nb)
# saveRDS(info.test.plot.1d, "info_test_plot_1d.rds")

####------------------------- plots --------------------------####
library(Cairo)
#item parameters for text (caption)
blocks <- create.block.ind(K, nb)
loads.blocks <- t(apply(blocks, 1, function(b, dl) colSums(dl[b,]), dl=load.mat))
loads.blocks

loads.blocks[4,]
items$u.mean[blocks[4,]]

#blockinfo
info.block.plot <- readRDS("plots/info_block_plot.rds")
CairoPDF("plots/plot_blockinfo.pdf", width=8, height=7, pointsize=12)
plot.block(which.blocks = 4, info=info.block.plot, K=20, loads=load.mat)
dev.off()

#testinfo
info.test.plot <- readRDS("plots/info_test_plot.rds")
CairoPDF("plots/plot_testinfo.pdf", width=12, height=22, pointsize=20)
plot.testinfo(info.test.plot, loads=load.mat, par.mfrow=c(5,2), cex=2)
dev.off()

#testinfo 1D
info.test.plot.1d <- readRDS("plots/info_test_plot_1d.rds")
se.lower <- info.test.plot.1d$ses[,1] - info.test.plot.1d$SDses[,1]
se.upper <- info.test.plot.1d$ses[,1] + info.test.plot.1d$SDses[,1]
#plot with base r or anything else
y.lim <- range(info.test.plot.1d$ses, na.rm=T)
CairoPDF("plots/plot_testinfo_1d.pdf", width=8, height=8, pointsize=20)
plot(info.test.plot.1d$variedlevels, info.test.plot.1d$ses[,1], type="l", ylim=y.lim, ylab="SE",xlab=expression(theta), xlim=c(-1.8,1.8))
grid(col="lightgrey", lty="solid")
polygon(x=c(info.test.plot.1d$variedlevels,rev(info.test.plot.1d$variedlevels)), c(info.test.plot.1d$ses[,1], rev(se.upper)), col="skyblue1")
polygon(x=c(info.test.plot.1d$variedlevels,rev(info.test.plot.1d$variedlevels)), c(info.test.plot.1d$ses[,1], rev(se.lower)), col="skyblue1")
lines(info.test.plot.1d$variedlevels, info.test.plot.1d$ses[,1], lwd=2)
dev.off()
