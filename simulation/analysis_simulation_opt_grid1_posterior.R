####--------------- analysis simulation testinfo ---------------####

library(psych)
devtools::load_all()
devtools::load_all("../DataAnalysisSimulation/")

####------------------ read in and check results ---------------####

res <- readRDS("simulation/results_simulation_opt_grid1_posterior.rds")
head(res)

anyNA(res)

#all conditions included?
table(res$blocksize)
table(res$keying)
table(res$length)
table(res$intercepts)
table(res$loads)
table(res$constraints)
table(res$target)
table(res$ntraits)
table(res$algorithm)

####--------------- data preparation ----------------------####


#remove keying, loads and length (they were not varied)
res <- res[, ! colnames(res) %in% c("keying","loads","length","constraints","ntraits", "target")]
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
res$algorithm <- factor(res$algorithm, levels=c("greedy-a","greedy-d","opt-t","r2","mean-a","loads","random"))
head(res$algorithm)

####------------------ descriptives in text -------------####

# mean.frame(dvs = c("rec", "RMSE", "MAB", "rel", "fisherz.r"), ivs = c("algorithm","blocksize", "intercepts"), results = res, rn = 3)
# mean.frame(dvs = c("sens", "spec"), ivs = c("algorithm","blocksize", "intercepts"), results = res, rn = 3)
# mean.frame(dvs = c("A","D","T.opt","Frob"), ivs = c("algorithm","blocksize", "intercepts"), results = res, rn = 3)
# 
# mean.frame(dvs = c("rec", "RMSE", "MAB", "rel", "fisherz.r"), ivs = c("algorithm"), results = res, rn = 3)
# mean.frame(dvs = c("sens", "spec"), ivs = c("algorithm"), results = res, rn = 3)
# mean.frame(dvs = c("A","D","T.opt","Frob"), ivs = c("algorithm"), results = res, rn = 3)
# 
# mean.frame(dvs = c("A","D","T.opt","Frob"), ivs = c("blocksize"), results = res, rn = 3)
# mean.frame(dvs = c("A","D","T.opt","Frob"), ivs = c("intercepts"), results = res, rn = 3)


# algorithm vs. random
# round(mean(res$MAB[res$algorithm %in% c("opt","r2","loads")]), 2)
# round(mean(res$MAB[res$algorithm %in% c("random")]), 2)
# 
# # Info vs. Loadings
# round(mean(res$MAB[res$algorithm %in% c("opt","r2")]), 2)
# round(mean(res$MAB[res$algorithm %in% c("loads")]), 2)
# 

#blocksize
# round(mean(res$T.opt[res$blocksize %in% 2]), 2)
# round(mean(res$T.opt[res$blocksize %in% 3]), 2)
# round(mean(res$T.opt[res$blocksize %in% 4]), 2)
# 
# #intercepts
# round(mean(res$T.opt[res$intercepts %in% "random"]), 2)
# round(mean(res$T.opt[res$intercepts %in% "ordered"]), 2)

# mean table

#without constraints because they did not make a difference and the table gets too long otherwise
means.paper <- mean.frame(dvs=c("sens", "spec"), ivs=c("algorithm","intercepts"), results=res, na.rm=T, rn=2)
means.paper$algorithm <- recode.df(means.paper$algorithm, 
                                   c("greedy-a","greedy-d","opt-t","mean-a","r2","loads","random"), 
                                   c("Greedy Variances", "Greedy Determinant", "MIP Trace", "Mean Variances","Block $R^2$","Mean Loadings","Random"))
means.paper$intercepts <- recode.df(means.paper$intercepts, c("ordered", "random"), c("Ordered", "Random"))

#re-order columns
factors <- 1:2
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
header$command <- c("\\hline \n Intercepts & Algorithm & \\multicolumn{2}{c}{Sensitivity} & \\multicolumn{2}{c}{Specificity} \\\\",
                    "\\hline \n \\multicolumn{6}{l}{\\small \\textit{Note.} MIP = Mixed Integer Programming. Standard deviations are} \\\\ 
                    \n \\multicolumn{6}{l}{\\small given in parentheses.} \n")
print(xtable::xtable(means.paper, digits=2,
                     caption="Mean sensitivity and specificity by condition in simulation study 2 on test construction for the single target (screening test)",
                     label="tb:means_screen"), include.colnames = F, include.rownames=F, 
      hline.after=seq(0, nrow(means.paper)-1, by = length(unique(res$algorithm))),
      sanitize.rownames.function=function(x){x}, sanitize.colnames.function = function(x){x},
      sanitize.text.function = function(x){x},
      NA.string = "", table.placement = "htp", add.to.row = header,
      caption.placement = "top", latex.environments = NULL,
      #floating = TRUE, floating.environment = "sidewaystable",
      file="../../Projekte/MFC_blocks/paper/Revision3_Psychometrika/SOM/textable_means_screening.tex")

# descriptives in text

# specificity overall mean
round(mean(res$spec), 2)
round(sd(res$spec), 2)

# sensitivity overall mean
round(mean(res$sens), 2)

# within random intercepts
round(mean(res$sens[(res$intercepts %in% "random") & (res$algorithm %in% c("random"))]), 2)
round(mean(res$sens[(res$intercepts %in% "random") & (res$algorithm %in% c("loads", "mean-a"))]), 2)
round(mean(res$sens[(res$intercepts %in% "random") & (res$algorithm %in% c("greedy-a","greedy-d","r2","opt-t"))]), 2)


# within ordered intercepts
# random vs. mip.t + variances, rest
round(mean(res$sens[(res$intercepts %in% "ordered") & (res$algorithm %in% c("random"))]), 2)
round(mean(res$sens[(res$intercepts %in% "ordered") & (res$algorithm %in% c("opt-t", "mean-a"))]), 2)
round(mean(res$sens[(res$intercepts %in% "ordered") & (res$algorithm %in% c("greedy-a","greedy-d","r2", "loads"))]), 2)


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
                                                 c("algovsrandom","infovsloadings","algovsmeans","r2vsmeansa","optvsgreedy", "greedyavsd")))

head(res$blocksize)
res$blocksize <- as.factor(res$blocksize)
head(res$blocksize)
contrasts(res$blocksize)
contrasts(res$blocksize) <- matrix(c(2,-1,-1,0,1,-1), 3, 2,
                                   dimnames=list(c("2","3","4"),
                                                 c("2vs34", "3vs4")))

head(res$intercepts)
res$intercepts <- factor(res$intercepts, levels = c("random", "ordered"))
contrasts(res$intercepts) <- matrix(c(-1,1), 2, 1,
                                    dimnames=(list(c("random","ordered"),
                                                   c("orderedvsrandom"))))
contrasts(res$intercepts)

# save prepared data with correct contrasts on factors
saveRDS(res, file = "simulation/results_opt_grid1_posterior_cleaned.rds")

# only for sensitivity and specificity
lm.sens <- calc.lms(dvs=c("sens","spec"), ivs=c("algorithm", "blocksize", "intercepts"), results=res)
var.expl.sens <- var.expl(lm.sens)
rm.0rows(var.expl.sens)

#format for latex
var.sens <- var.expl.sens
var.sens <- rm.0rows(var.sens)
colnames(var.sens) <- c("Sensitivity", "Specificity")
rownames(var.sens) <- c("Algorithm vs. Random", "Block $R^2$ vs. Mean Variances",
                        "2 vs. 3 and 4", "3 vs. 4",
                        "Intercepts",
                        "Algorithm vs. Random $\\times$ Intercepts",
                        "2 vs. 3 and 4 $\\times$ Intercepts",
                        "Residuals")
var.sens

header <- list()
header$pos <- list(-1, nrow(var.sens))
header$command <- c("\\hline \n Factor & Sens. & Spec. \\\\",
                    "\\hline \n \\multicolumn{3}{l}{\\small \\textit{Note.} Sens = Sensitivity, Spec = Specificity.} \n")

print(xtable::xtable(var.sens, digits=0,
                     caption="Variance in sensitivity and specificity explained in \\% by algorithm, intercepts and block size in simulation study 2 on test construction for the single target (screening test)",
                     label="tb:var_screen"), include.colnames = F, include.rownames=T, hline.after=c(0, nrow(var.sens)-1),
      sanitize.rownames.function=function(x){x}, sanitize.colnames.function = function(x){x},
      sanitize.text.function = function(x){x},
      NA.string = "", table.placement = "htp", add.to.row = header,
      caption.placement = "top", latex.environments = NULL,
      file="../../Projekte/MFC_blocks/paper/Revision3_Psychometrika/manuscript_Psychometrika_Revision3/textable_var_screening.tex")


# lm.algo.main <- calc.lms.main(dvs=c("fisherz.r","MAB","RMSE"), ivs=c("algorithm","target"), results=res)
lm.algo.main <- calc.lms.main(dvs=c("sens","spec","A","D","T.opt","Frob"), ivs=c("algorithm", "blocksize", "intercepts"), results=res)
var.expl(lm.algo.main)

# lm.algo <- calc.lms(dvs=c("fisherz.r","MAB","RMSE"), ivs=c("algorithm","target"), results=res)
lm.algo <- calc.lms(dvs=c("sens","spec","A","D","T.opt","Frob"), ivs=c("algorithm", "blocksize", "intercepts"), results=res)
var.expl.algo <- var.expl(lm.algo)
rm.0rows(var.expl.algo)




####---------------- plots --------------####

library(ggplot2)
library(gridExtra)
library(colorspace)

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

res.b3 <- res[res$blocksize == "3" & res$intercepts == "ordered",]

plot.sens <- plot.algo("sens", "Sensitivity", res.b3)
plot.A <- plot.algo("A", "A-optimality", res.b3)
plot.Topt <- plot.algo("T.opt", "T-optimality", res.b3)

plot.spec <- plot.algo("spec", "Specificity", res.b3)
plot.D <- plot.algo("D", "D-optimality", res.b3)
plot.Frob <- plot.algo("Frob", "Frobenius Norm Testinfo", res.b3)

ggsave("../../Projekte/MFC_blocks/paper/Revision3_Psychometrika/manuscript_Psychometrika_Revision3/plot_screening.pdf",
       grid.arrange(plot.sens, plot.spec,
                    nrow=1, ncol=2),
       width=14, height=4, units="in")

ggsave("simulation/plot_opt_grid1_posterior.pdf",
       grid.arrange(plot.sens, plot.A, plot.Topt,
                    plot.spec, plot.D, plot.Frob,
                    nrow=2, ncol=3),
       width=20, height=8, units="in")

plot.rec <- plot.algo("rec", expression(r(theta,hat(theta))), res.b3)
plot.MAB <- plot.algo("MAB", "MAB", res.b3)
plot.RMSE <- plot.algo("RMSE", "RMSE", res.b3)

ggsave("simulation/plot_opt_grid1_posterior_recovery.pdf",
       grid.arrange(plot.rec, plot.MAB, plot.RMSE, nrow=1, ncol=3),
       width=18, height=4, units="in")
