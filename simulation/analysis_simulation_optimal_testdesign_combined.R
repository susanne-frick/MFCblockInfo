####--------------- analysis simulation testinfo ---------------####
library(psych)
devtools::load_all()
devtools::load_all("../../DataAnalysisSimulation/")

####------------------ read in and check results ---------------####

res.simple <- readRDS("simulation/results_simulation_opt.rds")
res.con <- readRDS("simulation/results_simulation_opt_con.rds")
res.equal <- readRDS("simulation/results_simulation_opt_equal.rds")
res.equal.con <- readRDS("simulation/results_simulation_opt_equal_con.rds")
#combine
res <- rbind(res.simple, res.con, res.equal, res.equal.con)
#add identifying columns
res$target <- c(rep("weighted", nrow(res.simple)+nrow(res.con)), rep("equal", nrow(res.equal)+nrow(res.equal.con)))
res$constraints <- c(rep("free", nrow(res.simple)), rep("constrained", nrow(res.con)),
                     rep("free", nrow(res.equal)), rep("constrained", nrow(res.equal.con)))
res$target <- factor(res$target)
res$constraints <- factor(res$constraints)

nrow(res)/(5*4*4)
head(res)
anyNA(res)
#all finished


describe(res[,-c(1:8)])
#looks like no difference

####--------------- data preparation ----------------------####

#remove keying and blocksize (they were not varied)
res <- res[, ! colnames(res) %in% c("keying","blocksize","intercepts","loads","length")]
head(res)

#remove rows with convergence issues
# res <- res[rowSums(is.na(res))==0,]

#add true reliability
res$rel <- res$rec^2

#add Fisher Z of r(true,est)
res$fisherz.r <- fisherz(res$rec)

#move columns target and constrained to the beginning
res <- res[,c("rep","trait","algorithm","target","constraints","rec","fisherz.r","rel","RMSE","MAB")]

head(res)
tail(res)

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

####------------------ absolute recovery ----------------####

means.paper <- mean.frame(dvs=c("rec","rel","MAB","RMSE"), ivs=c("algorithm","constraints","target"), results=res, na.rm=T, rn=2)
means.paper$algorithm <- recode.df(means.paper$algorithm, c("opt","r2","loads","random"), c("MIP T","Block $R^2$","Mean Loadings","Random"))
means.paper$target <- recode.df(means.paper$target, c("equal","weighted"), c("Equal","Weighted"))
means.paper$constraints <- recode.df(means.paper$constraints, c("constrained","free"), c("Constrained","Free"))
means.paper <- add.brackets.SDs(means.paper)

#move target and constraints to the beginning
means.paper <- means.paper[,c("target","constraints","algorithm",grep("SD[.]|M[.]",colnames(means.paper), value=T))]
#recode duplicated values in target to empty spaces
means.paper[,"target"][duplicated(means.paper[,"target"])] <- ""
means.paper[,"constraints"][-c(seq(1,nrow(means.paper), by=length(unique(res$algorithm))))] <- ""

means.paper

header <- list()
header$pos <- list(-1, nrow(means.paper))
header$command <- c("\\hline \n Algorithm & Constraints & Target & \\multicolumn{2}{c}{$r(\\theta,\\hat{\\theta})$} & \\multicolumn{2}{c}{$r(\\theta,\\hat{\\theta})^2$} & \\multicolumn{2}{c}{MAB} & \\multicolumn{2}{c}{RMSE} \\\\",
                    "\\hline \n \\multicolumn{11}{l}{\\small \\textit{Note.} MAB = Mean Absolute Bias, RMSE = Root Mean Squared Error, MIP T = Mixed Integer} \\\\
                    \\multicolumn{11}{l}{\\small Programming with T-optimality. Standard deviations are given in brackets.} \n")
print(xtable::xtable(means.paper, digits=2,
                     caption="Mean trait recovery by condition in simulation study 2 on test construction",
                     label="tb:means_con"), include.colnames = F, include.rownames=F, hline.after=c(0,4,8,12),
      sanitize.rownames.function=function(x){x}, sanitize.colnames.function = function(x){x},
      sanitize.text.function = function(x){x},
      NA.string = "", table.placement = "htp", add.to.row = header,
      caption.placement = "top", latex.environments = NULL,
      floating = TRUE, floating.environment = "sidewaystable",
      file="../paper/textable_means_construction.tex")


####------------------ differences between algorithms -----------####

misty::multilevel.icc(res[,c("RMSE","MAB","fisherz.r")], group = res$rep)
misty::multilevel.icc(res[,c("RMSE","MAB","fisherz.r")], group = res$trait)
#variance due to replication and trait is negligible

contrasts(res$algorithm)
contrasts(res$algorithm) <- matrix(c(1,1,1,-3, 1,1,-2,0, 1,-1,0,0), 4, 3,
                                   dimnames=list(c("opt","r2","loads","random"),
                                                 c("algovsrandom","infovsloadings","optvsr2")))
contrasts(res$target)
contrasts(res$target) <- matrix(c(-1,1), 2, 1,
                                dimnames=(list(c("equal","weighted"),
                                              c("weightedvsequal"))))
contrasts(res$constraints)
contrasts(res$constraints) <- matrix(c(-1,1), 2, 1,
                                     dimnames=list(c("constrained","free"),
                                                   c("freevscontrained")))

lm.algo.main <- calc.lms.main(dvs=c("fisherz.r","MAB","RMSE"), ivs=c("algorithm","constraints","target"), results=res)
var.expl(lm.algo.main)

lm.algo <- calc.lms(dvs=c("fisherz.r","MAB","RMSE"), ivs=c("algorithm","constraints","target"), results=res)
var.expl.algo <- var.expl(lm.algo)
var.expl.algo

#format for latex
var.paper <- var.expl.algo
var.paper <- rm.0rows(var.paper)
colnames(var.paper) <- c("Fisher Z($r(\theta, \hat{\theta})$)","MB","RMSE")
rownames(var.paper) <- c("Algorithm vs. Random", "Info vs. Mean Loadings","Weighted vs. Equal","Residuals")
var.paper

header <- list()
header$pos <- list(-1, nrow(var.paper))
header$command <- c("\\hline \n Factor & $r(\\theta, \\hat{\\theta})$ & MAB & RMSE\\\\",
                    "\\hline \n \\multicolumn{4}{l}{\\small \\textit{Note.} MAB = Mean Absolute Bias, RMSE = Root} \\\\
                    \\multicolumn{4}{l}{\\small Mean Squared Error, MIP T = Mixed Integer} \\\\
                    \\multicolumn{4}{l}{\\small Programming with T-optimality. $r(\\theta, \\hat{\\theta})$ was Fisher \\textit{Z}}\\\\
                    \\multicolumn{4}{l}{\\small transformed.} \n")

print(xtable::xtable(var.paper, digits=0,
                     caption="Variance in recovery of trait scores explained in \\% by algorithm, target and constraints in simulation study 2 on test construction",
                     label="tb:var_sim_con"), include.colnames = F, include.rownames=T, hline.after=c(0, nrow(var.paper)-1),
      sanitize.rownames.function=function(x){x}, sanitize.colnames.function = function(x){x},
      sanitize.text.function = function(x){x},
      NA.string = "", table.placement = "htp", add.to.row = header,
      caption.placement = "top", latex.environments = NULL,
      file="../paper/textable_var_construction.tex")

####---------------- plots --------------####
library(ggplot2)
library(gridExtra)
library(colorspace)

res.equal.con <- res[res$target=="equal" & res$constraints=="constrained",]

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


ggsave("../paper/plot_optimal_testdesign.pdf",
       grid.arrange(plot.MAB, plot.RMSE, plot.rec, nrow=1, ncol=3),
       width=10, height=4, units="in")

