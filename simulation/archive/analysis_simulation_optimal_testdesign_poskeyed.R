####--------------- analysis simulation testinfo ---------------####
setwd("/home/frick/Dokumente")
setwd("/media/susanne/Volume/frick/Arbeit_TU")

library(psych)
devtools::load_all()
devtools::load_all("../DataAnalysisSimulation/")

####------------------ read in and check results ---------------####

res <- readRDS("simulation/results_simulation_opt_poskeyed.rds")
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

# describe(res[,-c(1:8)])
#looks like no difference

nrow(res)/(5*4)
#there should be 2400
#4 are missing, probably did not finish in time

#number of NAs (optimizer did not find a solution)
table(rowSums(is.na(res))>0)
res[rowSums(is.na(res))>0,]
#3 tests for block size 2 #same as with mixed item keying!

mean.frame(dvs=c("rec","MAB","RMSE"), ivs=c("algorithm","target","constraints","blocksize"), results=res)

####--------------- data preparation ----------------------####

#remove keying, intercepts, loads and length (they were not varied)
res <- res[, ! colnames(res) %in% c("keying","intercepts","loads","length")]
#also remove ntraits because they did not finish in time (actually only 5 traits)
res <- res[, ! colnames(res) %in% c("ntraits")]
head(res)

#remove rows with convergence issues
res <- res[rowSums(is.na(res))==0,]

#add true reliability
res$rel <- res$rec^2

#add Fisher Z of r(true,est)
res$fisherz.r <- fisherz(res$rec)

#move columns target and constrained to the beginning
res <- res[,c("rep","blocksize","trait","algorithm","target","constraints","rec","fisherz.r","rel","RMSE","MAB")]

head(res)
tail(res)

#make algorithm a factor
str(res)
head(res$algorithm)
res$algorithm <- factor(res$algorithm, levels=c("opt","r2","loads","random"))
head(res$algorithm)

####------------------ descriptives in text -------------####
# algorithm vs. random
round(mean(res$MAB[res$algorithm %in% c("opt","r2","loads")]), 2)
round(mean(res$MAB[res$algorithm %in% c("random")]), 2)

# Info vs. Loadings
round(mean(res$MAB[res$algorithm %in% c("opt","r2")]), 2)
round(mean(res$MAB[res$algorithm %in% c("loads")]), 2)

# Equal vs. Weighted
round(mean(res$MAB[res$target %in% c("equal")]), 2)
round(mean(res$MAB[res$target %in% c("weighted")]), 2)

round(mean(res$rec[res$target %in% c("equal")]), 2)
round(mean(res$rec[res$target %in% c("weighted")]), 2)

#blocksize
round(mean(res$MAB[res$blocksize %in% 2]), 2)
round(mean(res$MAB[res$blocksize %in% 3]), 2)
round(mean(res$MAB[res$blocksize %in% 4]), 2)

#number of traits
# round(mean(res$MAB[res$ntraits %in% 5]), 2)
# round(mean(res$MAB[res$ntraits %in% 15]), 2)

####------------------ absolute recovery ----------------####

#without constraints because they did not make a difference and the table gets too long otherwise
means.paper <- mean.frame(dvs=c("rec","rel","MAB","RMSE"), ivs=c("algorithm","target","blocksize"), results=res, na.rm=T, rn=2)
means.paper$algorithm <- recode.df(means.paper$algorithm, c("opt","r2","loads","random"), c("MIP T","Block $R^2$","Mean Loadings","Random"))
means.paper$target <- recode.df(means.paper$target, c("equal","weighted"), c("Equal","Weighted"))
# means.paper$constraints <- recode.df(means.paper$constraints, c("constrained","unconstrained"), c("Constrained","Unconstrained"))

#re-order columns
factors <- 1:3
nlevels <- apply(means.paper[,factors], 2, function(cl) length(unique(cl)))

#only every first occurence of factor level
for(cl in factors) {
  if(cl > 1) {
    means.paper[,cl] <- as.character(means.paper[,cl])
    means.paper[-c(seq(1, nrow(means.paper), prod(nlevels[1:(cl-1)]))), cl] <- ""
  }
}
means.paper[,factors] <- means.paper[,rev(factors)]
#add SDs
means.paper <- add.brackets.SDs(means.paper)

means.paper

header <- list()
header$pos <- list(-1, nrow(means.paper))
header$command <- c("\\hline \n Blocksize & Algorithm & Target & \\multicolumn{2}{c}{$r(\\theta,\\hat{\\theta})$} & \\multicolumn{2}{c}{$r(\\theta,\\hat{\\theta})^2$} & \\multicolumn{2}{c}{MAB} & \\multicolumn{2}{c}{RMSE} \\\\",
                    "\\hline \n \\multicolumn{11}{l}{\\small \\textit{Note.} MAB = Mean Absolute Bias, RMSE = Root Mean Squared Error, MIP T = Mixed Integer} \\\\
                    \\multicolumn{11}{l}{\\small Programming with T-optimality. Standard deviations are given in parentheses.} \n")
print(xtable::xtable(means.paper, digits=2,
                     caption="Mean trait recovery by condition in simulation study 2 on test construction",
                     label="tb:means_con"), include.colnames = F, include.rownames=F, hline.after=c(0,4,8,12,16,20),
      sanitize.rownames.function=function(x){x}, sanitize.colnames.function = function(x){x},
      sanitize.text.function = function(x){x},
      NA.string = "", table.placement = "htp", add.to.row = header,
      caption.placement = "top", latex.environments = NULL,
      floating = TRUE, floating.environment = "sidewaystable",
      file="Projekte/MFC_blocks/paper/Revision1_Psychometrika/files_submission_Psychometrika_Revision1/textable_means_construction.tex")


####------------------ differences between algorithms -----------####

misty::multilevel.icc(res[,c("RMSE","MAB","fisherz.r")], group = res$rep)
misty::multilevel.icc(res[,c("RMSE","MAB","fisherz.r")], group = res$trait)
#variance due to replication and trait is negligible

contrasts(res$algorithm)
contrasts(res$algorithm) <- matrix(c(1,1,1,-3, 1,1,-2,0, 1,-1,0,0), 4, 3,
                                   dimnames=list(c("opt","r2","loads","random"),
                                                 c("algovsrandom","infovsloadings","optvsr2")))
contrasts(res$target)
contrasts(res$target) <- matrix(c(1,-1), 2, 1,
                                dimnames=(list(c("weighted","equal"),
                                              c("weightedvsequal"))))
contrasts(res$constraints)
contrasts(res$constraints) <- matrix(c(1,-1), 2, 1,
                                     dimnames=list(c("unconstrained","constrained"),
                                                   c("unconstrainedvscontrained")))

head(res$blocksize)
res$blocksize <- as.factor(res$blocksize)
head(res$blocksize)
contrasts(res$blocksize)
contrasts(res$blocksize) <- matrix(c(2,-1,-1,0,1,-1), 3, 2,
                                   dimnames=list(c("2","3","4"),
                                                 c("2vs34", "3vs4")))

lm.algo.main <- calc.lms.main(dvs=c("fisherz.r","MAB","RMSE"), ivs=c("algorithm","constraints","target","blocksize"), results=res)
var.expl(lm.algo.main)

lm.algo <- calc.lms(dvs=c("fisherz.r","MAB","RMSE"), ivs=c("algorithm","constraints","target","blocksize"), results=res)
var.expl.algo <- var.expl(lm.algo)
rm.0rows(var.expl.algo)

#format for latex
var.paper <- var.expl.algo
var.paper <- rm.0rows(var.paper)
colnames(var.paper) <- c("Fisher Z($r(\theta, \hat{\theta})$)","MB","RMSE")
rownames(var.paper) <- c("Algorithm vs. Random", "Info vs. Mean Loadings","Weighted vs. Equal","2 vs. 3 and 4", "3 vs. 4","Residuals")
var.paper

header <- list()
header$pos <- list(-1, nrow(var.paper))
header$command <- c("\\hline \n Factor & $r(\\theta, \\hat{\\theta})$ & MAB & RMSE\\\\",
                    "\\hline \n \\multicolumn{4}{l}{\\small \\textit{Note.} MAB = Mean Absolute Bias, RMSE = Root} \\\\
                    \\multicolumn{4}{l}{\\small Mean Squared Error, MIP T = Mixed Integer} \\\\
                    \\multicolumn{4}{l}{\\small Programming with T-optimality. $r(\\theta, \\hat{\\theta})$ was Fisher \\textit{Z}}\\\\
                    \\multicolumn{4}{l}{\\small transformed.} \n")

print(xtable::xtable(var.paper, digits=0,
                     caption="Variance in trait recovery explained in \\% by algorithm, target and constraints in simulation study 2 on test construction",
                     label="tb:var_sim_con"), include.colnames = F, include.rownames=T, hline.after=c(0, nrow(var.paper)-1),
      sanitize.rownames.function=function(x){x}, sanitize.colnames.function = function(x){x},
      sanitize.text.function = function(x){x},
      NA.string = "", table.placement = "htp", add.to.row = header,
      caption.placement = "top", latex.environments = NULL,
      file="Projekte/MFC_blocks/paper/Revision1_Psychometrika/files_submission_Psychometrika_Revision1/textable_var_construction.tex")

####---------------- plots --------------####
library(ggplot2)
library(gridExtra)
library(colorspace)

res.equal.con <- res[res$target=="equal" & res$constraints=="constrained" & res$blocksize==3,]

plot.MAB <- ggplot(data=res.equal.con, aes(y=MAB, x=algorithm, fill=algorithm)) +
  geom_violin(show.legend=FALSE) +
  labs(y="MAB", x="Algorithm") +
  scale_x_discrete(labels = c('opt' = "MIP T",
                              'r2'   = expression(plain(Block)~R^2),
                              "loads" = "Loadings",
                              "random" = "Random")) +
  scale_fill_manual(values=qualitative_hcl(6)[3:6]) +
  theme(axis.text=element_text(size=11),
        axis.title=element_text(size=11))
plot.RMSE <- ggplot(data=res.equal.con, aes(y=RMSE, x=algorithm, fill=algorithm)) +
  geom_violin(show.legend=FALSE) +
  labs(y="RMSE", x="Algorithm") +
  scale_x_discrete(labels = c('opt' = "MIP T",
                              'r2'   = expression(plain(Block)~R^2),
                              "loads" = "Loadings",
                              "random" = "Random")) +
  scale_fill_manual(values=qualitative_hcl(6)[3:6]) +
  theme(axis.text=element_text(size=11),
        axis.title=element_text(size=11))
plot.rec <- ggplot(data=res.equal.con, aes(y=rec, x=algorithm, fill=algorithm)) +
  geom_violin(show.legend=FALSE) +
  labs(y=expression(r(theta,hat(theta))), x="Algorithm") +
  scale_x_discrete(labels = c('opt' = "MIP T",
                              'r2'   = expression(plain(Block)~R^2),
                              "loads" = "Loadings",
                              "random" = "Random")) +
  scale_fill_manual(values=qualitative_hcl(6)[3:6]) +
  theme(axis.text=element_text(size=11),
        axis.title=element_text(size=11))


ggsave("Projekte/MFC_blocks/paper/Revision1_Psychometrika/files_submission_Psychometrika_Revision1/plot_optimal_testdesign.pdf",
       grid.arrange(plot.MAB, plot.RMSE, plot.rec, nrow=1, ncol=3),
       width=10, height=4, units="in")

