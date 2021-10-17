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

####------------------------- plots --------------------------####
library(Cairo)
#item parameters for text (caption)
blocks <- create.block.ind(K, nb)
loads.blocks <- t(apply(blocks, 1, function(b, dl) colSums(dl[b,]), dl=load.mat))
loads.blocks

loads.blocks[18,]
items$u.mean[blocks[18,]]

#blockinfo
info.block.plot <- readRDS("plots/info_block_plot.rds")
CairoPDF("plots/plot_blockinfo_Diss.pdf", width=8, height=7, pointsize=12)
plot.block(which.blocks = 18, info=info.block.plot, K=20, loads=load.mat)
dev.off()

file.copy(from="plots/plot_blockinfo_Diss.pdf", to="../../Dissertation/figures/plot_blockinfo_Diss.pdf", overwrite = TRUE)

