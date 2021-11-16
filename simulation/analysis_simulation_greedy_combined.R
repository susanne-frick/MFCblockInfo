####--------------- analysis simulation testinfo ---------------####
library(psych)
devtools::load_all()
devtools::load_all("../../DataAnalysisSimulation/")

####------------------ read in and check results ---------------####

#because 1/400 replications did not converge, results were saved separately
#and combined with combine_results_optimal_testdesign.R

res.simple <- readRDS("simulation/results_simulation_greedy.rds")
res.equal <-  readRDS("simulation/results_simulation_greedy_equal.rds")

#combine
res <- rbind(res.simple, res.equal)
#add identifying columns
res$target <- c(rep("weighted", nrow(res.simple)), rep("equal", nrow(res.equal)))
res$target <- factor(res$target)

nrow(res)/(5*6)
head(res)
anyNA(res)
#all finished


describe(res[,-c(1:7)])
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

#move column target to the beginning
res <- res[,c("rep","trait","algorithm","target","rec","fisherz.r","rel","RMSE","MAB")]


head(res)
tail(res)

####------------------ descriptives in text -------------####
# algorithm vs. random
round(mean(res$MAB[! res$algorithm %in% c("random")]), 2)
round(mean(res$MAB[res$algorithm %in% c("random")]), 2)


# Greedy + R2 vs. MIP T + loads
round(mean(res$MAB[res$algorithm %in% c("greedy-a","greedy-d","r2")]), 2)
round(mean(res$MAB[res$algorithm %in% c("opt","loads")]), 2)

#var: greedy vs. others
round(sd(res$MAB[res$algorithm %in% c("greedy-a","greedy-d")]), 2)
round(sd(res$MAB[res$algorithm %in% c("opt","loads","r2")]), 2)

# Equal vs. Weighted
round(mean(res$MAB[res$target %in% c("equal")]), 2)
round(mean(res$MAB[res$target %in% c("weighted")]), 2)

round(mean(res$rec[res$target %in% c("equal")]), 2)
round(mean(res$rec[res$target %in% c("weighted")]), 2)

####------------------ absolute recovery ----------------####

means.paper <- mean.frame(dvs=c("rec","rel","MAB","RMSE"), ivs=c("algorithm","target"), results=res, na.rm=T, rn=2)

means.paper$algorithm <- recode.df(means.paper$algorithm,
                                   c("greedy-a","greedy-d","opt-t","r2","loads","random"),
                                   c("Greedy A","Greedy D","MIP T","Block $R^2$","Mean Loadings","Random"))
means.paper$target <- recode.df(means.paper$target, c("equal","weighted"), c("Equal","Weighted"))

means.paper <- add.brackets.SDs(means.paper)
#change order of target and algorithm
means.paper <- means.paper[,c("target","algorithm",grep("SD[.]|M[.]",colnames(means.paper), value=T))]
#recode duplicated values in target to empty spaces
means.paper[,"target"][duplicated(means.paper[,"target"])] <- ""

means.paper

header <- list()
header$pos <- list(-1, nrow(means.paper))
header$command <- c("\\hline \n Target & Algorithm & \\multicolumn{2}{c}{$r(\\theta,\\hat{\\theta})$} & \\multicolumn{2}{c}{$r(\\theta,\\hat{\\theta})^2$} & \\multicolumn{2}{c}{MAB} & \\multicolumn{2}{c}{RMSE} \\\\",
                    "\\hline \n \\multicolumn{10}{l}{\\small \\textit{Note.} MAB = Mean Absolute Bias, RMSE = Root Mean Squared Error, A = A-optimality,} \\\\
                    \\multicolumn{10}{l}{\\small D = D-optimality, MIP T = Mixed Integer Programming with T-optimality. Standard} \\\\
                    \\multicolumn{10}{l}{\\small deviations are given in parentheses.}\n")

print(xtable::xtable(means.paper, digits=2,
                     caption="Mean trait recovery by algorithm in simulation study 3 on test extension",
                     label="tb:means_ext"), include.colnames = F, include.rownames=F, hline.after=c(0,6),
      sanitize.rownames.function=function(x){x}, sanitize.colnames.function = function(x){x},
      sanitize.text.function = function(x){x},
      NA.string = "", table.placement = "htp", add.to.row = header,
      caption.placement = "top", latex.environments = NULL,
      file="../paper/textable_means_extension.tex")


####------------------ differences between algorithms -----------####
contrasts(res$algorithm)
contrasts(res$algorithm) <- matrix(c(rep(1,5),-5, rep(1,4),-4,0, rep(1,3),-3,rep(0,2), rep(1,2),-2,rep(0,3), 1,-1,rep(0,4)), 6, 5,
                                   dimnames=list(c("greedy-a","greedy-d","opt-t","r2","loads","random"),
                                                 c("algovsrandom","infovsloadings","optgreedyvsr2","greedyvsopt","avsd")))

contrasts(res$target)
contrasts(res$target) <- matrix(c(-1,1), 2, 1,
                                dimnames=(list(c("equal","weighted"),
                                               c("weightedvsequal"))))

lm.algo.main <- calc.lms.main(dvs=c("fisherz.r","MAB","RMSE"), ivs=c("algorithm","target"), results=res)
var.expl(lm.algo.main)

lm.algo <- calc.lms(dvs=c("fisherz.r","MAB","RMSE"), ivs=c("algorithm","target"), results=res)
var.expl.algo <- var.expl(lm.algo)
var.expl.algo

#format for latex
var.paper <- var.expl.algo
var.paper <- rm.0rows(var.paper)
rownames(var.paper) <- c("Algorithm vs. Random", "Info vs. Mean Loadings","Weighted vs. Equal","Residuals")
var.paper

header <- list()
header$pos <- list(-1, nrow(var.paper))
header$command <- c("\\hline \n Factor & $r(\\theta, \\hat{\\theta})$ & MAB & RMSE\\\\",
                    "\\hline \n \\multicolumn{4}{l}{\\small \\textit{Note.} MAB = Mean Absolute Bias, RMSE = Root} \\\\
                    \\multicolumn{4}{l}{\\small Mean Squared Error. $r(\\theta, \\hat{\\theta})$ was Fisher \\textit{Z} } \\\\
                    \\multicolumn{4}{l}{\\small transformed.}\n")

print(xtable::xtable(var.paper, digits=0,
                     caption="Variance in trait recovery explained in \\% by algorithm and target in simulation study 3 on test extension",
                     label="tb:var_sim_ext"), include.colnames = F, include.rownames=T, hline.after=c(0, nrow(var.paper)-1),
      sanitize.rownames.function=function(x){x}, sanitize.colnames.function = function(x){x},
      sanitize.text.function = function(x){x},
      NA.string = "", table.placement = "htp", add.to.row = header,
      caption.placement = "top", latex.environments = NULL,
      file="../paper/textable_var_extension.tex")

####---------------- plots --------------####
library(ggplot2)
library(gridExtra)
library(colorspace)

res.equal <- res[res$target=="equal",]

plot.MAB <- ggplot(data=res.equal, aes(y=MAB, x=algorithm, fill=algorithm)) +
  geom_violin(show.legend=FALSE) +
  labs(y="MAB", x="Algorithm") +
  scale_x_discrete(labels = c('greedy-a' = "Greedy A",
                              'greedy-d' = "Greedy D",
                              'opt-t' = "MIP T",
                              'r2'   = expression(plain(Block)~R^2),
                              "loads" = "Loadings",
                              "random" = "Random")) +
  scale_fill_manual(values=qualitative_hcl(6)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=16))
plot.RMSE <- ggplot(data=res.equal, aes(y=RMSE, x=algorithm, fill=algorithm)) +
  geom_violin(show.legend=FALSE) +
  labs(y="RMSE", x="Algorithm") +
  scale_x_discrete(labels = c('greedy-a' = "Greedy A",
                              'greedy-d' = "Greedy D",
                              'opt-t' = "MIP T",
                              'r2'   = expression(plain(Block)~R^2),
                              "loads" = "Loadings",
                              "random" = "Random")) +
  scale_fill_manual(values=qualitative_hcl(6)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=16))
plot.rec <- ggplot(data=res.equal, aes(y=rec, x=algorithm, fill=algorithm)) +
  geom_violin(show.legend=FALSE) +
  labs(y=expression(r(theta,hat(theta))), x="Algorithm") +
  scale_x_discrete(labels = c('greedy-a' = "Greedy A",
                              'greedy-d' = "Greedy D",
                              'opt-t' = "MIP T",
                              'r2'   = expression(plain(Block)~R^2),
                              "loads" = "Loadings",
                              "random" = "Random")) +
  scale_fill_manual(values=qualitative_hcl(6)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=16))
ggsave("../paper/plot_greedy.pdf",
       grid.arrange(plot.MAB, plot.RMSE, plot.rec, nrow=2, ncol=2),
       width=15, height=12, units="in")

####------------------- plot Dissertation ------------------####
res.weighted <- res[res$target=="weighted",]

plot.rec <- ggplot(data=res.weighted, aes(y=rec, x=algorithm, fill=algorithm)) +
  geom_violin(show.legend=FALSE) +
  labs(y=expression(r(eta,hat(eta))), x="Algorithm") +
  scale_x_discrete(labels = c('greedy-a' = "Greedy A",
                              'greedy-d' = "Greedy D",
                              'opt-t' = "MIP T",
                              'r2'   = expression(plain(Block)~R^2),
                              "loads" = "Loadings",
                              "random" = "Random")) +
  scale_fill_manual(values=qualitative_hcl(6)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=16))
ggsave(filename="/media/susanne/Volume/frick/ideas_tests/Dissertation/figures/plot_greedy_weighted.pdf",
       plot=plot.rec,
       width=7.5, height=5, units="in")
