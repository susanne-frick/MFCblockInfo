####--------------- analysis simulation testinfo ---------------####

library(psych)
devtools::load_all()
devtools::load_all("../DataAnalysisSimulation/")

####------------------ read in and check results ---------------####

# combine both versions
res_ML <- readRDS("simulation/results_simulation_opt_grid1.rds")
res_posterior <- readRDS("simulation/results_simulation_opt_grid1_posterior_r2-avg.rds")

# remove rows with r2-avg and r2-avg-sd
res_posterior <- res_posterior[!(res_posterior$algorithm %in% c("r2-avg", "r2-avg-sd", "greedy-a", "greedy-d", "mean-a")),]
table(res_posterior$algorithm)
table(res_ML$algorithm)
# recode opt-t as opt
res_posterior$algorithm[res_posterior$algorithm == "opt-t"] <- "opt"

res <- rbind(res_ML, res_posterior)
res$likelihood <- c(rep("ML", nrow(res_ML)), rep("posterior", nrow(res_posterior)))
rm(res_ML, res_posterior)

head(res)

anyNA(res)

#all conditions included?
table(res$blocksize, useNA="always")
table(res$keying)
table(res$length)
table(res$intercepts)
table(res$loads)
table(res$constraints)
table(res$target)
table(res$ntraits)
table(res$likelihood)

####--------------- data preparation ----------------------####

#remove keying, intercepts, loads and length (they were not varied)
res <- res[, ! colnames(res) %in% c("keying","intercepts","loads","length","constraints","ntraits")]
head(res)

#remove rows with convergence issues
res <- res[rowSums(is.na(res))==0,]

#add true reliability
res$rel <- res$rec^2

#add Fisher Z of r(true,est)
res$fisherz.r <- fisherz(res$rec)

head(res)
tail(res)

#make algorithm a factor
str(res)
head(res$algorithm)
res$algorithm <- factor(res$algorithm, levels=c("opt","r2","loads","random"))
head(res$algorithm)

####------------------ descriptives in text -------------####

mean.frame(dvs = c("rec", "RMSE", "MAB", "rel", "fisherz.r"), ivs = c("algorithm","blocksize"), results = res, rn = 3)
mean.frame(dvs = c("sens", "spec"), ivs = c("algorithm","blocksize"), results = res, rn = 3)
mean.frame(dvs = c("A","D","T.opt","Frob"), ivs = c("algorithm","blocksize"), results = res, rn = 3)

mean.frame(dvs = c("rec", "RMSE", "MAB", "rel", "fisherz.r"), ivs = c("algorithm"), results = res, rn = 3)
mean.frame(dvs = c("sens", "spec"), ivs = c("algorithm"), results = res, rn = 3)
mean.frame(dvs = c("A","D","T.opt","Frob"), ivs = c("algorithm"), results = res, rn = 3)


# algorithm vs. random
# round(mean(res$MAB[res$algorithm %in% c("opt","r2","loads")]), 2)
# round(mean(res$MAB[res$algorithm %in% c("random")]), 2)
# 
# # Info vs. Loadings
# round(mean(res$MAB[res$algorithm %in% c("opt","r2")]), 2)
# round(mean(res$MAB[res$algorithm %in% c("loads")]), 2)
# 
# # Equal vs. Weighted
# round(mean(res$MAB[res$target %in% c("equal")]), 2)
# round(mean(res$MAB[res$target %in% c("weighted")]), 2)
# 
# round(mean(res$rec[res$target %in% c("equal")]), 2)
# round(mean(res$rec[res$target %in% c("weighted")]), 2)
# 
# #blocksize
# round(mean(res$MAB[res$blocksize %in% 2]), 2)
# round(mean(res$MAB[res$blocksize %in% 3]), 2)
# round(mean(res$MAB[res$blocksize %in% 4]), 2)

#number of traits
# round(mean(res$MAB[res$ntraits %in% 5]), 2)
# round(mean(res$MAB[res$ntraits %in% 15]), 2)

####------------------ absolute recovery ----------------####

# #without constraints because they did not make a difference and the table gets too long otherwise
# means.paper <- mean.frame(dvs=c("rec","rel","MAB","RMSE"), ivs=c("algorithm","target","blocksize"), results=res, na.rm=T, rn=2)
# means.paper$algorithm <- recode.df(means.paper$algorithm, c("opt","r2","loads","random"), c("MIP T","Block $R^2$","Mean Loadings","Random"))
# means.paper$target <- recode.df(means.paper$target, c("equal","weighted"), c("Equal","Weighted"))
# # means.paper$constraints <- recode.df(means.paper$constraints, c("constrained","unconstrained"), c("Constrained","Unconstrained"))
# 
# #re-order columns
# factors <- 1:3
# nlevels <- apply(means.paper[,factors], 2, function(cl) length(unique(cl)))
# 
# #only every first occurence of factor level
# for(cl in factors) {
#   if(cl > 1) {
#     means.paper[,cl] <- as.character(means.paper[,cl])
#     means.paper[-c(seq(1, nrow(means.paper), prod(nlevels[1:(cl-1)]))), cl] <- ""
#   }
# }
# means.paper[,factors] <- means.paper[,rev(factors)]
# #add SDs
# means.paper <- add.brackets.SDs(means.paper)
# 
# means.paper
# 
# header <- list()
# header$pos <- list(-1, nrow(means.paper))
# header$command <- c("\\hline \n Blocksize & Algorithm & Target & \\multicolumn{2}{c}{$r(\\theta,\\hat{\\theta})$} & \\multicolumn{2}{c}{$r(\\theta,\\hat{\\theta})^2$} & \\multicolumn{2}{c}{MAB} & \\multicolumn{2}{c}{RMSE} \\\\",
#                     "\\hline \n \\multicolumn{11}{l}{\\small \\textit{Note.} MAB = Mean Absolute Bias, RMSE = Root Mean Squared Error, MIP T = Mixed Integer} \\\\
#                     \\multicolumn{11}{l}{\\small Programming with T-optimality. Standard deviations are given in parentheses.} \n")
# print(xtable::xtable(means.paper, digits=2,
#                      caption="Mean trait recovery by condition in simulation study 2 on test construction",
#                      label="tb:means_con"), include.colnames = F, include.rownames=F, hline.after=c(0,4,8,12,16,20),
#       sanitize.rownames.function=function(x){x}, sanitize.colnames.function = function(x){x},
#       sanitize.text.function = function(x){x},
#       NA.string = "", table.placement = "htp", add.to.row = header,
#       caption.placement = "top", latex.environments = NULL,
#       floating = TRUE, floating.environment = "sidewaystable",
#       file="../../Projekte/MFC_blocks/paper/Revision1_Psychometrika/files_submission_Psychometrika_Revision1/textable_means_construction.tex")
# 

####------------------ differences between algorithms -----------####

# misty::multilevel.icc(res[,c("RMSE","MAB","fisherz.r")], group = res$rep)
# misty::multilevel.icc(res[,c("RMSE","MAB","fisherz.r")], group = res$trait)
#variance due to replication and trait is negligible

contrasts(res$algorithm)
contrasts(res$algorithm) <- matrix(c(1,1,1,-3, 1,1,-2,0, 1,-1,0,0), 4, 3,
                                   dimnames=list(c("opt","r2","loads","random"),
                                                 c("algovsrandom","infovsloadings","optvsr2")))

head(res$blocksize)
res$blocksize <- as.factor(res$blocksize)
head(res$blocksize)
contrasts(res$blocksize)
contrasts(res$blocksize) <- matrix(c(2,-1,-1,0,1,-1), 3, 2,
                                   dimnames=list(c("2","3","4"),
                                                 c("2vs34", "3vs4")))

head(res$likelihood)
res$likelihood <- as.factor(res$likelihood)
head(res$likelihood)
contrasts(res$likelihood)

# lm.algo.main <- calc.lms.main(dvs=c("fisherz.r","MAB","RMSE"), ivs=c("algorithm","target"), results=res)
lm.algo.main <- calc.lms.main(dvs=c("sens","spec","A","D","T.opt","Frob"), ivs=c("algorithm", "blocksize", "likelihood"), results=res)
var.expl(lm.algo.main)

# lm.algo <- calc.lms(dvs=c("fisherz.r","MAB","RMSE"), ivs=c("algorithm","target"), results=res)
lm.algo <- calc.lms(dvs=c("sens","spec","A","D","T.opt","Frob"), ivs=c("algorithm", "blocksize", "likelihood"), results=res)
var.expl.algo <- var.expl(lm.algo)
rm.0rows(var.expl.algo)

lm.algo.rec <- calc.lms(dvs=c("fisherz.r","MAB","RMSE"), ivs=c("algorithm", "blocksize"), results=res)
var.expl.algo.rec <- var.expl(lm.algo.rec)
rm.0rows(var.expl.algo.rec)

#format for latex
# var.paper <- var.expl.algo
# var.paper <- rm.0rows(var.paper)
# colnames(var.paper) <- c("Fisher Z($r(\theta, \hat{\theta})$)","MB","RMSE")
# rownames(var.paper) <- c("Algorithm vs. Random", "Info vs. Mean Loadings","Weighted vs. Equal","2 vs. 3 and 4", "3 vs. 4","Residuals")
# var.paper
# 
# header <- list()
# header$pos <- list(-1, nrow(var.paper))
# header$command <- c("\\hline \n Factor & $r(\\theta, \\hat{\\theta})$ & MAB & RMSE\\\\",
#                     "\\hline \n \\multicolumn{4}{l}{\\small \\textit{Note.} MAB = Mean Absolute Bias, RMSE = Root} \\\\
#                     \\multicolumn{4}{l}{\\small Mean Squared Error, MIP T = Mixed Integer} \\\\
#                     \\multicolumn{4}{l}{\\small Programming with T-optimality. $r(\\theta, \\hat{\\theta})$ was Fisher \\textit{Z}}\\\\
#                     \\multicolumn{4}{l}{\\small transformed.} \n")
# 
# print(xtable::xtable(var.paper, digits=0,
#                      caption="Variance in trait recovery explained in \\% by algorithm, target and constraints in simulation study 2 on test construction",
#                      label="tb:var_sim_con"), include.colnames = F, include.rownames=T, hline.after=c(0, nrow(var.paper)-1),
#       sanitize.rownames.function=function(x){x}, sanitize.colnames.function = function(x){x},
#       sanitize.text.function = function(x){x},
#       NA.string = "", table.placement = "htp", add.to.row = header,
#       caption.placement = "top", latex.environments = NULL,
#       file="../../Projekte/MFC_blocks/paper/Revision1_Psychometrika/files_submission_Psychometrika_Revision1/textable_var_construction.tex")

####---------------- plots --------------####
library(ggplot2)
library(gridExtra)
library(colorspace)

plot.algo <- function(y, ylab, data) {
  ggplot(data=data, aes(y=get(y), x=algorithm, fill=algorithm)) +
    geom_violin(show.legend=FALSE) +
    labs(y=ylab, x="Algorithm") +
    scale_x_discrete(labels = c('opt' = "MIP T",
                                'r2'   = expression(plain(Block)~R^2),
                                "loads" = "Loadings",
                                "random" = "Random")) +
    scale_fill_manual(values=qualitative_hcl(6)[3:6]) +
    theme(axis.text=element_text(size=11),
          axis.title=element_text(size=11))
}

res.b3 <- res[res$blocksize == "3",]

plot.sens <- plot.algo("sens", "Sensitivity", res.b3)
plot.A <- plot.algo("A", "A-optimality", res.b3)
plot.Topt <- plot.algo("T.opt", "T-optimality", res.b3)

plot.spec <- plot.algo("spec", "Specificity", res.b3)
plot.D <- plot.algo("D", "D-optimality", res.b3)
plot.Frob <- plot.algo("Frob", "Frobenius Norm Testinfo", res.b3)

ggsave("simulation/plot_opt_grid1.pdf",
       grid.arrange(plot.sens, plot.A, plot.Topt,
                    plot.spec, plot.D, plot.Frob,
                    nrow=2, ncol=3),
       width=10, height=8, units="in")

plot.rec <- plot.algo("rec", expression(r(theta,hat(theta))), res.b3)
plot.MAB <- plot.algo("MAB", "MAB", res.b3)
plot.RMSE <- plot.algo("RMSE", "RMSE", res.b3)

ggsave("simulation/plot_opt_grid1_recovery.pdf",
       grid.arrange(plot.rec, plot.MAB, plot.RMSE, nrow=1, ncol=3),
       width=10, height=4, units="in")

