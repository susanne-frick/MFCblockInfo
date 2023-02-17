####---------------- analysis simulation optimal testdesign across all three targets: weighted, equal, single -----####

#####---------------------------- combine with population target ------------------------------------####

res_pop <- readRDS("simulation/results_opt_int-variance_cleaned.rds")
res_screen <- readRDS("simulation/results_opt_grid1_posterior_cleaned.rds")


res_screen$target <- "single"
cl_names <- c("rep", "blocksize", "intercepts", "target", "trait", "algorithm", "A", "D", "T.opt", "Frob")
res_both <- rbind(res_screen[, cl_names], res_pop[, cl_names])

res_both$target <- factor(res_both$target, levels = c("weighted", "equal", "single"))
contrasts(res_both$target)
contrasts(res_both$target) <- matrix(c(1,1,-2,1,-1,0), 3, 2,
                                     dimnames=list(c("weighted", "equal", "single"),
                                                   c("populationvsscreening", "weightedvsequal")))

contrasts(res_both$algorithm)
contrasts(res_both$algorithm) <- matrix(c(1,1,1,1,1,1,-6, 
                                          1,1,1,1,1,-5,0,
                                          1,1,1,-1.5,-1.5,0,0,
                                          0,0,0,1,-1,0,0,
                                          -1,-1,2,0,0,0,0,
                                          1,-1,0,0,0,0,0), 
                                        7, 6,
                                        dimnames=list(c("greedy-a", "greedy-d","opt-t","r2","mean-a","loads","random"),
                                                      c("algovsrandom","infovsloadings","algovsmeans","r2vsmeansa","optvsgreedy", "greedyavsd")))

contrasts(res_both$blocksize)
contrasts(res_both$blocksize) <- matrix(c(2,-1,-1,0,1,-1), 3, 2,
                                        dimnames=list(c("2","3","4"),
                                                      c("2vs34", "3vs4")))

contrasts(res_both$intercepts)
contrasts(res_both$intercepts) <- matrix(c(-1,1), 2, 1,
                                         dimnames=(list(c("random","ordered"),
                                                        c("orderedvsrandom"))))

####--------------------- explained variance -----------------------------------------------####

lm.algo.both <- calc.lms(dvs=c("A","D","T.opt","Frob"), ivs=c("algorithm", "target", "blocksize", "intercepts"), results=res_both)
var.expl.algo.both <- var.expl(lm.algo.both)

rm.1rows <- function(df) {
  df[rowSums(df<=1)<ncol(df),]
}

rm.1rows(var.expl.algo.both)

#format for latex
var.paper <- var.expl.algo.both
var.paper <- rm.1rows(var.paper)
colnames(var.paper) <- c("Variances", "Determinant", "Trace", "Frobenius Norm")
rownames(var.paper) <- c("Algorithm vs. Random", "Info vs. Mean Loadings",
                         "Population vs. Screening", "Weighted vs. Equal",
                         "2 vs. 3 and 4", "3 vs. 4",
                         "Intercepts",
                         "Info vs. Mean Loadings $\\times$ Population vs. Screening",
                         "Algorithm vs. Random $\\times$ 2 vs. 3 and 4",
                         "Population vs. Screening $\\times$ 2 vs. 3 and 4",
                         "Weighted vs. Equal $\\times$ 2 vs. 3 and 4",
                         "Algorithm vs. Random $\\times$ Intercepts",
                         "Population vs. Screening $\\times$ Intercepts",
                         "Weighted vs. Equal $\\times$ Intercepts",
                         "2 vs. 3 and 4 $\\times$ Intercepts",
                         "Residuals")
var.paper

header <- list()
header$pos <- list(-1, nrow(var.paper))
header$command <- c("\\hline \n Factor & Var. & Det. & Trace & Frob.\\\\",
                    "\\hline \n \\multicolumn{4}{l}{\\small \\textit{Note.} Var = Variances, Det = Determinant, Frob = Frobenius.} \n")

print(xtable::xtable(var.paper, digits=0,
                     caption="Variance in optimization criteria explained in \\% by algorithm, target, intercepts and block size in simulation study 2 on test construction",
                     label="tb:var_opt"), include.colnames = F, include.rownames=T, hline.after=c(0, nrow(var.paper)-1),
      sanitize.rownames.function=function(x){x}, sanitize.colnames.function = function(x){x},
      sanitize.text.function = function(x){x},
      NA.string = "", table.placement = "htp", add.to.row = header,
      caption.placement = "top", latex.environments = NULL,
      file="../../Projekte/MFC_blocks/paper/Revision2_Psychometrika/manuscript_Psychometrika_Revision2/textable_var_optimization.tex")



####-------------------------------- descriptives in text ------------####

mean.frame(dvs = c("A","D","T.opt","Frob"), ivs = c("algorithm"), results = res_both)

# algorithm vs. random
round(mean(res_both$A[res_both$algorithm %in% c("greedy-a","greedy-d","opt-t","r2","mean-a","loads")]), 2)
round(mean(res_both$A[res_both$algorithm %in% c("random")]), 2)

# algorithm vs. random in ordered vs. random
round(mean(res_both$A[(res_both$intercepts %in% "random") & (res_both$algorithm %in% c("greedy-a","greedy-d","opt-t","r2","mean-a","loads"))]), 2) -
  round(mean(res_both$A[(res_both$intercepts %in% "random") & (res_both$algorithm %in% c("random"))]), 2)

round(mean(res_both$A[(res_both$intercepts %in% "ordered") & (res_both$algorithm %in% c("greedy-a","greedy-d","opt-t","r2","mean-a","loads"))]), 2) -
  round(mean(res_both$A[(res_both$intercepts %in% "ordered") & (res_both$algorithm %in% c("random"))]), 2)

res_both$algovsrandom <- ifelse(res_both$algorithm == "random", "random", "algorithm")
mean.frame(dvs = c("A","D","T.opt","Frob"), ivs = c("algovsrandom", "intercepts"), results = res_both)

# algorithm vs. random in 2 vs. 3 and 4
mean.frame(dvs = c("A","D","T.opt","Frob"), ivs = c("algovsrandom", "blocksize"), results = res_both)

round(mean(res_both$Frob[(res_both$blocksize %in% 2) & (res_both$algovsrandom == "random")]), 2) -
  round(mean(res_both$Frob[(res_both$blocksize %in% 2) & (res_both$algovsrandom == "algorithm")]), 2)

round(mean(res_both$Frob[(res_both$blocksize %in% 3:4) & (res_both$algovsrandom == "random")]), 2) -
  round(mean(res_both$Frob[(res_both$blocksize %in% 3:4) & (res_both$algovsrandom == "algorithm")]), 2)


# Info vs. Loadings
round(mean(res_both$A[res_both$algorithm %in% c("greedy-a","greedy-d","opt-t","r2","mean-a")]), 2)
round(mean(res_both$A[res_both$algorithm %in% c("loads")]), 2)

####----------------- table of means -------------------------------####

# summarize target as population (weighted, equal) vs. screening (single)
res_both$target2 <- ifelse(res_both$target == "single", "Screening", "Population")
res_both$target2 <- factor(res_both$target2, levels = c("Population", "Screening"))

#without constraints because they did not make a difference and the table gets too long otherwise
means.paper <- mean.frame(dvs=c("A","D","T.opt","Frob"), ivs=c("algorithm","target2","intercepts"), results=res_both, na.rm=T, rn=2)
means.paper$algorithm <- recode.df(means.paper$algorithm, 
                                   c("greedy-a","greedy-d","opt-t","mean-a","r2","loads","random"), 
                                   c("Greedy Variances", "Greedy Determinant", "MIP Trace", "Mean Variances","Block $R^2$","Mean Loadings","Random"))
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
header$command <- c("\\hline \n Intercepts & Target & Algorithm & \\multicolumn{2}{c}{Variances} & \\multicolumn{2}{c}{Determinant} & \\multicolumn{2}{c}{Trace} & \\multicolumn{2}{c}{Frobenius Norm} \\\\",
                    "\\hline \n \\multicolumn{11}{l}{\\small \\textit{Note.} MIP = Mixed Integer Programming. Standard deviations are given in parentheses.} \n")
print(xtable::xtable(means.paper, digits=2,
                     caption="Mean optimization criteria by condition in simulation study 2 on test construction",
                     label="tb:means_opt"), include.colnames = F, include.rownames=F, 
      hline.after=seq(0, nrow(means.paper)-1, by = length(unique(res_both$algorithm))),
      sanitize.rownames.function=function(x){x}, sanitize.colnames.function = function(x){x},
      sanitize.text.function = function(x){x},
      NA.string = "", table.placement = "htp", add.to.row = header,
      caption.placement = "top", latex.environments = NULL,
      floating = TRUE, floating.environment = "sidewaystable",
      file="../../Projekte/MFC_blocks/paper/Revision2_Psychometrika/manuscript_Psychometrika_Revision2/textable_means_optimization.tex")

####---------------- plots both ---------------------####

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

res.b3.ordered.population <- res_both[res_both$blocksize == "3" & res_both$intercepts == "ordered" & res_both$target == "equal",]

plot.A.pop <- plot.algo("A", "Variances", res.b3.ordered.population)
plot.D.pop <- plot.algo("D", "Determinant", res.b3.ordered.population)
plot.Topt.pop <- plot.algo("T.opt", "Trace", res.b3.ordered.population)
plot.Frob.pop <- plot.algo("Frob", "Frobenius Norm Testinfo", res.b3.ordered.population)

res.b3.ordered.screening <- res_both[res_both$blocksize == "3" & res_both$intercepts == "ordered" & res_both$target == "single",]

plot.A.scr <- plot.algo("A", "Variances", res.b3.ordered.screening)
plot.D.scr <- plot.algo("D", "Determinant", res.b3.ordered.screening)
plot.Topt.scr <- plot.algo("T.opt", "Trace", res.b3.ordered.screening)
plot.Frob.scr <- plot.algo("Frob", "Frobenius Norm Testinfo", res.b3.ordered.screening)

ggsave("../../Projekte/MFC_blocks/paper/Revision2_Psychometrika/manuscript_Psychometrika_Revision2/plot_optimality.pdf",
       grid.arrange(plot.A.pop, plot.D.pop, plot.Topt.pop, plot.Frob.pop,
                    plot.A.scr, plot.D.scr, plot.Topt.scr, plot.Frob.scr,
                    nrow=2, ncol=4),
       width=26, height=8, units="in")
