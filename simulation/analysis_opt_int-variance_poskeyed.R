####--------------- analysis simulation testinfo ---------------####

library(psych)
devtools::load_all()
devtools::load_all("../DataAnalysisSimulation/")

####------------------ read in and check results ---------------####

res <- readRDS("simulation/results_opt_int-variance_poskeyed_2397.rds")
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
table(res$algorithm)

nrow(res)/(5*4)

#number of NAs (optimizer did not find a solution)
table(rowSums(is.na(res))>0)

# #identifier for which ones are missing
# id.res <- paste(res$blocksize, res$target, res$intercepts, res$rep, sep="_")
# 
# factor.blocksize <- 2:4
# factor.keying <- "12" # c("0","12","23")
# factor.int <- c("random", "ordered")
# factor.load <- "acceptable"
# factor.length <- "long"
# factor.algorithm <- c("greedy-a","greedy-d","opt-t","mean-a","r2","loads","random")
# factor.constraints <- c("unconstrained")
# factor.target <- c("weighted","equal")
# factor.ntraits <- "5" #c("5","15")
# 
# #number of replications
# R <- 200 #500
# 
# design.sim <- expand.grid("blocksize"=factor.blocksize, "keying"=factor.keying, "length"=factor.length, "intercepts"=factor.int,
#                           "loads"=factor.load, 
#                           "constraints"=factor.constraints, "target"=factor.target, "ntraits"=factor.ntraits,
#                           "rep"=1:R)
# id.res.full <- paste(design.sim$blocksize, design.sim$target, design.sim$intercepts, design.sim$rep, sep="_")
# setdiff(id.res.full, id.res)
# which(id.res.full == setdiff(id.res.full, id.res))


####--------------- data preparation ----------------------####

#remove keying, intercepts, loads and length (they were not varied)
res <- res[, ! colnames(res) %in% c("keying","constraints","loads","length","ntraits")]
head(res)

#remove rows with convergence issues
res <- res[rowSums(is.na(res))==0,] # none in this case

#add true reliability
res$rel <- res$rec^2

#add Fisher Z of r(true,est)
res$fisherz.r <- fisherz(res$rec)

#move columns to the beginning
res <- res[,c("rep","blocksize","intercepts","target","trait","algorithm","rec","fisherz.r","rel","RMSE","MAB","A","D","T.opt","Frob")]

head(res)
tail(res)

#make algorithm a factor
str(res)
head(res$algorithm)
res$algorithm <- factor(res$algorithm, levels=c("greedy-a","greedy-d","opt-t","r2","mean-a","loads","random"))
head(res$algorithm)

####------------------ descriptives in text -------------####

# algorithm vs. random
round(mean(res$MAB[res$algorithm %in% c("greedy-a","greedy-d","opt-t","r2","mean-a","loads")]), 2)
round(mean(res$MAB[res$algorithm %in% c("random")]), 2)

# Algo (Optimality) vs. means
round(mean(res$MAB[res$algorithm %in% c("greedy-a","greedy-d","opt-t")]), 2)
round(mean(res$MAB[res$algorithm %in% c("r2", "mean-a")]), 2)

# R2 vs. means-a
round(mean(res$MAB[res$algorithm %in% c("r2")]), 2)
round(mean(res$MAB[res$algorithm %in% c("mean-a")]), 2)


# Intercepts: orderd vs. random
round(mean(res$MAB[res$intercepts %in% c("ordered")]), 2)
round(mean(res$MAB[res$intercepts %in% c("random")]), 2)


# Target: Equal vs. Weighted
round(mean(res$MAB[res$target %in% c("equal")]), 2)
round(mean(res$MAB[res$target %in% c("weighted")]), 2)

round(mean(res$rec[res$target %in% c("equal")]), 2)
round(mean(res$rec[res$target %in% c("weighted")]), 2)

#blocksize
round(mean(res$MAB[res$blocksize %in% 2]), 2)
round(mean(res$MAB[res$blocksize %in% 3]), 2)
round(mean(res$MAB[res$blocksize %in% 4]), 2)

####------------------ absolute recovery ----------------####

#without blocksize because the table gets too long otherwise and block size is more the effect of total information
means.paper <- mean.frame(dvs=c("rec","rel","MAB","RMSE"), ivs=c("algorithm","target","intercepts"), results=res, na.rm=T, rn=2)
means.paper$algorithm <- recode.df(means.paper$algorithm, 
                                   c("greedy-a","greedy-d","opt-t","mean-a","r2","loads","random"), 
                                   c("Greedy Variances", "Greedy Determinant", "MIP Trace", "Mean Variances","Block $R^2$","Mean Loadings","Random"))
means.paper$target <- recode.df(means.paper$target, c("equal","weighted"), c("Equal","Weighted"))
means.paper$intercepts <- recode.df(means.paper$intercepts, c("ordered", "random"), c("Ordered", "Random"))

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
header$command <- c("\\hline \n Intercepts & Target & Algorithm & \\multicolumn{2}{c}{$r(\\theta,\\hat{\\theta})$} & \\multicolumn{2}{c}{$r(\\theta,\\hat{\\theta})^2$} & \\multicolumn{2}{c}{MAB} & \\multicolumn{2}{c}{RMSE} \\\\",
                    "\\hline \n \\multicolumn{11}{l}{\\small \\textit{Note.} MAB = Mean Absolute Bias, RMSE = Root Mean Squared Error, MIP = Mixed Integer Programming.} \\\\
                    \\multicolumn{11}{l}{\\small Standard deviations are given in parentheses.} \n")
print(xtable::xtable(means.paper, digits=2,
                     caption="Mean trait recovery by condition in the simulation on test construction with all positively keyed items for the weighted and equal target (population test)",
                     label="tb:means_rec_pos"), include.colnames = F, include.rownames=F, 
      hline.after=seq(0, nrow(means.paper)-1, by = length(unique(res$algorithm))),
      sanitize.rownames.function=function(x){x}, sanitize.colnames.function = function(x){x},
      sanitize.text.function = function(x){x},
      NA.string = "", table.placement = "htp", add.to.row = header,
      caption.placement = "top", latex.environments = NULL,
      floating = TRUE, floating.environment = "sidewaystable",
      file="../../Projekte/MFC_blocks/paper/Revision2_Psychometrika/SOM/textable_means_population_poskeyed.tex")


####------------------ differences between algorithms -----------####

# misty::multilevel.icc(res[,c("RMSE","MAB","fisherz.r")], group = res$rep)
# misty::multilevel.icc(res[,c("RMSE","MAB","fisherz.r")], group = res$trait)
#variance due to replication and trait is negligible

contrasts(res$algorithm)
contrasts(res$algorithm) <- matrix(c(1,1,1,1,1,1,-6, 
                                     1,1,1,1,1,-5,0,
                                     1,1,1,-1.5,-1.5,0,0,
                                     0,0,0,1,-1,0,0,
                                     -1,-1,2,0,0,0,0,
                                     1,-1,0,0,0,0,0), 
                                   7, 6,
                                   dimnames=list(c("greedy-a", "greedy-d","opt-t","r2","mean-a","loads","random"),
                                                 c("algovsrandom","infovsloadings","algovsmeans","r2vsmeana","optvsgreedy", "greedyavsd")))
contrasts(res$target)
contrasts(res$target) <- matrix(c(1,-1), 2, 1,
                                dimnames=(list(c("weighted","equal"),
                                              c("weightedvsequal"))))

contrasts(res$intercepts)
contrasts(res$intercepts) <- matrix(c(-1,1), 2, 1,
                                dimnames=(list(c("random","ordered"),
                                               c("orderedvsrandom"))))


head(res$blocksize)
res$blocksize <- as.factor(res$blocksize)
head(res$blocksize)
contrasts(res$blocksize)
contrasts(res$blocksize) <- matrix(c(2,-1,-1,0,1,-1), 3, 2,
                                   dimnames=list(c("2","3","4"),
                                                 c("2vs34", "3vs4")))

# save prepared data with correct contrasts on factors
saveRDS(res, file = "simulation/results_opt_int-variance_poskeyed_cleaned.rds")

# lm.algo.main <- calc.lms.main(dvs=c("fisherz.r","MAB","RMSE"), ivs=c("algorithm","intercepts","target","blocksize"), results=res)
# var.expl(lm.algo.main)

lm.algo <- calc.lms(dvs=c("fisherz.r","MAB","RMSE"), ivs=c("algorithm","intercepts","target","blocksize"), results=res)
var.expl.algo <- var.expl(lm.algo)
rm.0rows(var.expl.algo)

#format for latex
var.paper <- var.expl.algo
var.paper <- rm.0rows(var.paper)
colnames(var.paper) <- c("Fisher Z($r(\\theta, \\hat{\\theta})$)","MAB","RMSE")
rownames(var.paper) <- c("Algorithm vs. Random", 
                         "Info vs. Loadings", 
                         "Algorithm vs. Means", 
                         "$R^2$ vs. Mean Variances", 
                         "Intercepts","Target","2 vs. 3 and 4", "3 vs. 4",
                         "Algorithm vs. Random $\\times$ Intercepts", 
                         "$R^2$ vs. Mean Variances $\\times$ Intercepts", 
                         "Target $\\times$ Intercepts",
                         "Algorithm vs. Random $\\times$ 2 vs. 3 and 4",
                         "$R^2$ vs. Mean Variances $\\times$ 2 vs. 3 and 4", 
                         "2 vs. 3 and 4 $\\times$ Intercepts",
                         "$R^2$ vs. Mean Variances $\\times$ Intercepts $\\times$ 2 vs. 3 and 4", 
                         "Residuals")
var.paper

header <- list()
header$pos <- list(-1, nrow(var.paper))
header$command <- c("\\hline \n Factor & $r(\\theta, \\hat{\\theta})$ & MAB & RMSE\\\\",
                    "\\hline \n \\multicolumn{4}{l}{\\small \\textit{Note.} MAB = Mean Absolute Bias, RMSE = Root Mean Squared Error.} \\\\
                    \\multicolumn{4}{l}{\\small $r(\\theta, \\hat{\\theta})$ was Fisher \\textit{Z} transformed.} \n")

print(xtable::xtable(var.paper, digits=0,
                     caption="Variance in trait recovery explained in \\% by algorithm, target and intercepts in the simulation on test construction with all positively keyed items for the weighted and equal target (population test)",
                     label="tb:var_rec_pos"), include.colnames = F, include.rownames=T, hline.after=c(0, nrow(var.paper)-1),
      sanitize.rownames.function=function(x){x}, sanitize.colnames.function = function(x){x},
      sanitize.text.function = function(x){x},
      NA.string = "", table.placement = "htp", add.to.row = header,
      caption.placement = "top", latex.environments = NULL,
      file="../../Projekte/MFC_blocks/paper/Revision2_Psychometrika/SOM/textable_var_population_poskeyed.tex")

####---------------- plots --------------####
library(ggplot2)
library(gridExtra)
library(colorspace)

# !weighted instead of equal
res.equal.ord <- res[res$target=="equal" & res$intercepts=="ordered" & res$blocksize==3,]

plot.algo <- function(y, ylab, data) {
  ggplot(data=data, aes(y=get(y), x=algorithm, fill=algorithm)) +
    geom_violin(show.legend=FALSE) +
    labs(y=ylab, x="Algorithm") +
    scale_x_discrete(labels = c("greedy-a" = "Variances",
                                "greedy-d" = "Determinant",
                                'opt-t' = "Trace",
                                'r2' = expression(plain(Block)~R^2),
                                "mean-a" = "M(Variances)",
                                "loads" = "Loadings",
                                "random" = "Random")) +
    scale_fill_manual(values=qualitative_hcl(7)) +
    theme(axis.text=element_text(size=11),
          axis.title=element_text(size=11))
}

plot.MAB <- plot.algo("MAB", "MAB", res.equal.ord)
plot.RMSE <- plot.algo("RMSE", "RMSE", res.equal.ord)
plot.rec <- plot.algo("rec", expression(r(theta,hat(theta))), res.equal.ord)
plot.A <- plot.algo("A", "A-optimality", res.equal.ord)
plot.D <- plot.algo("D", "D-optimality", res.equal.ord)
plot.T.opt <- plot.algo("T.opt", "T-optimality", res.equal.ord)
plot.Frob <- plot.algo("Frob", "Frobenius Norm", res.equal.ord)

ggsave("../../Projekte/MFC_blocks/paper/Revision2_Psychometrika/SOM/plot_opt_recovery_poskeyed.pdf",
       grid.arrange(plot.rec, plot.MAB, plot.RMSE,
                    nrow=1, ncol=3),
       width=20, height=6, units="in")

ggsave("simulation/plot_opt_int-variance_poskeyed_equal-ordered-3.pdf",
       grid.arrange(plot.A, plot.D, plot.T.opt, plot.Frob,
                    plot.rec, plot.MAB, plot.RMSE,
                    nrow=2, ncol=4),
       width=24, height=8, units="in")


res.equal.ord.b2 <- res[res$target=="equal" & res$intercepts=="ordered" & res$blocksize==2,]
plot.b2.MAB <- plot.algo("MAB", "MAB", res.equal.ord.b2)
plot.b2.RMSE <- plot.algo("RMSE", "RMSE", res.equal.ord.b2)
plot.b2.rec <- plot.algo("rec", expression(r(theta,hat(theta))), res.equal.ord.b2)
ggsave("../../Projekte/MFC_blocks/paper/Revision2_Psychometrika/SOM/plot_opt_recovery_poskeyed_B2.pdf",
       grid.arrange(plot.b2.rec, plot.b2.MAB, plot.b2.RMSE,
                    nrow=1, ncol=3),
       width=20, height=6, units="in")

