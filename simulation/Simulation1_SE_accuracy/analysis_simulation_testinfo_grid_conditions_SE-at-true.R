####--------------- analysis simulation testinfo ---------------####
library(psych)
library(gridExtra)
library(ggplot2)
devtools::load_all("~/Dokumente/packages/MFCblockInfo/")
devtools::load_all("~/Dokumente/packages/DataAnalysisSimulation/")

####--------------- read in and check results ------------------####
setwd("~/Dokumente/Projekte/MFC_blocks/MFCblockInfo_simulations/")
res <- readRDS("results_simulation_testinfo_grid_dependent_conditions_SE-at-true.rds")
head(res)


table(res$n.err)
table(res$n.warn)
table(res$n.mess)
table(res$n.mess, res$likelihood) #all produced a message, but probably that is just telling how it converged

table(res$na.info) #total number of NAs
#0 NAs
table(res$na.info>0)/nrow(res) #percentage

# describe(res)

####---------------- reformat results --------------------------####

trait.seq <- seq(-2,2,.5) #trait levels that were varied
#add level of Trait 2
res$trait.level <- recode.df(res$p, 1:(length(trait.seq)*2), rep(trait.seq, 2))

#are all conditions included?
table(res$blocksize)
table(res$loadings)
table(res$length)
table(res$estimator)
table(res$likelihood)

# recode likelihood to new labels: genuine and independence
res$likelihood <- recode.df(res$likelihood, c("true", "dependent"), c("genuine", "independence"))

#change order of factor levels to be consistent with previous simulation
levels(res$estimator)
res$estimator <- factor(res$estimator, levels=c("ML","MAP"))
res$loadings <- factor(res$loadings, levels=c("low","high"))

#separate results for trait 2
res2 <- res[res$trait==2,]

#to long format: observed (at estimate) vs. expected (weighted with probability)
head(res)
# add empty columns rmse.exp, var.ratio.exp
res$rmse.exp <- res$var.ratio.exp <- rep(NA, nrow(res))
fixed <- c("likelihood","loadings","estimator","length","blocksize","p","trait.level","q","trait")
res.long <- rbind(res[,c(fixed,"mean.bias.obs","rmse.obs","mean.ratio.obs","var.ratio.obs")],
                  setNames(res[,c(fixed,"bias.exp","rmse.exp","ratio.exp","var.ratio.exp")], 
                           c(fixed,"mean.bias.obs","rmse.obs","mean.ratio.obs","var.ratio.obs")))
colnames(res.long) <- c(fixed,"mean.bias","rmse","mean.ratio","var.ratio")
#add indicator for method: expected vs. observed
res.long$method <- rep(c("observed","expected"), each=nrow(res))
#re-order columns
res.long <- res.long[,c(fixed,"method","mean.bias","rmse","mean.ratio","var.ratio")]
head(res.long)
tail(res.long)

#separate results for Trait 2
res.long.trait2 <- res.long[res.long$trait==2,]

####------------------------------- plots -----------------------------------------####

#length: short, long
#loading: high, low
#blocksize: 2,3,4
dvs <- c("se.true","mean.bias.est","rmse.est")
dv.labs <- c("SE","MB","RMSE")
names(dv.labs) <- dvs

#additional variable: length x loadings
res2$length.loadings <- paste(res2$length, res2$loadings, sep=" - ")
res2$length.loadings <- factor(res2$length.loadings, levels=c("short - low", "short - high", "long - low", "long - high"))

ggplot.means <- function(dv, ivs, ivx, res, dvlabs, size=16) {  
    form <- formula(paste(dv, "~", paste(c(ivx, ivs, "length.loadings","blocksize"), collapse=" + ")))
    
    means <- aggregate(form, FUN=mean, data=res2)
    sds <- aggregate(form, FUN=sd, data=res2)
    means$lower <- means[,dv] - sds[,dv]
    means$upper <- means[,dv] + sds[,dv]
    
    ggplot(data=means, aes(x=get(ivx), y=get(dv), ymin=lower, ymax=upper,
                           group=interaction(means[,c("likelihood","estimator")]), fill=estimator, linetype=likelihood)) +
      geom_line() +
      geom_ribbon(alpha=0.5) +
      xlab(as.expression(expression(theta))) +
      ylab(dvlabs[dv]) +
      theme(axis.text=element_text(size=size),
            axis.title=element_text(size=size),
            legend.text=element_text(size=size),
            legend.title=element_text(size=size)) +
      facet_grid(rows=vars(length.loadings), cols=vars(blocksize))
  }

plot.list <- lapply(dvs, ggplot.means, ivs=c("likelihood","estimator"), ivx="trait.level", res=res2, dvlabs=dv.labs, size=10)

ggsave("plot_recovery_dependent_conditions_SE-at-true.pdf",
       plot=lemon::grid_arrange_shared_legend(plot.list[[1]], plot.list[[2]], plot.list[[3]],
                                              nrow=1, ncol=3, position="right")
       ,width=30, height=16, units="cm")
file.copy(from="plot_recovery_dependent_conditions_SE-at-true.pdf",
          to="../paper/Revision3_Psychometrika/manuscript_Psychometrika_Revision3/plot_recovery.pdf",
          overwrite = TRUE)

#bias observed
plot.list.obs <- lapply(c("mean.bias.obs","rmse.obs", "mean.ratio.obs", "var.ratio.obs"), 
                        ggplot.means, ivs=c("likelihood","estimator"), ivx="trait.level", res=res2, 
                        dvlabs=c("mean.bias.obs"="MB", "rmse.obs"="RMSE", "mean.ratio.obs"="Mean Ratio" , "var.ratio.obs"="Var Ratio"), 
                        size=10)

ggsave("plot_SEs_observed_conditions_SE-at-true.pdf",
       plot=lemon::grid_arrange_shared_legend(plot.list.obs[[1]], plot.list.obs[[2]], plot.list.obs[[3]],
                                              nrow=1, ncol=3, position="right")
       ,width=30, height=16, units="cm")
file.copy(from="plot_SEs_observed_conditions_SE-at-true.pdf",
          to="../paper/Revision3_Psychometrika/manuscript_Psychometrika_Revision3/plot_SEs_observed.pdf",
          overwrite = TRUE)

#bias expected
plot.list.exp <- lapply(c("bias.exp", "ratio.exp"), ggplot.means, ivs=c("likelihood","estimator"), ivx="trait.level", res=res2, 
                        dvlabs=c("bias.exp"="MB", "ratio.exp"="Ratio"), size=10)

ggsave("plot_SEs_expected_conditions_SE-at-true.pdf",
       plot=lemon::grid_arrange_shared_legend(plot.list.exp[[1]], plot.list.exp[[2]],
                                              nrow=1, ncol=2, position="right")
       ,width=20, height=16, units="cm")
file.copy(from="plot_SEs_expected_conditions_SE-at-true.pdf",
          to="../paper/Revision3_Psychometrika/manuscript_Psychometrika_Revision3/plot_SEs_expected.pdf",
          overwrite = TRUE)

####----------------- ANOVA --------------------####

#observed vs. expected
lm.info.main <- calc.lms.main(dvs=c("mean.bias","mean.ratio"), 
                              ivs=c("loadings","estimator","length","method","likelihood","blocksize"), results=res.long.trait2)
var.expl.info <- rm.0rows(var.expl(lm.info.main))

#RMSE, var.ratio: only observed
lm.info.obs <- calc.lms.main(dvs=c("rmse","var.ratio"), 
                              ivs=c("loadings","estimator","length","likelihood","blocksize"), 
                             results=res.long.trait2[res.long.trait2$method == "observed", ])
var.expl.info.obs <- rm.0rows(var.expl(lm.info.obs))

# combine to one dataframe (including NA and 0 rows)
var.expl.info.both <- merge(var.expl.info, var.expl.info.obs[, "rmse"], by = 0, all = TRUE)
# bias and ratio: NA = 0
var.expl.info.both$mean.bias <- car::recode(var.expl.info.both$mean.bias, "NA = 0")
var.expl.info.both$mean.ratio <- car::recode(var.expl.info.both$mean.ratio, "NA = 0")
# to usual structure with rownames
rownames(var.expl.info.both) <- var.expl.info.both$Row.names
var.expl.info.both <- var.expl.info.both[, -c(1)]
var.expl.info.both <- rm.0rows(var.expl.info.both)

#estimates
lm.est.main <- calc.lms.main(dvs=c("se.true","mean.bias.est","rmse.est"), ivs=c("loadings","estimator","length","likelihood","blocksize"), results=res2)
var.expl.est <- rm.0rows(var.expl(lm.est.main))
#at 0
rm.0rows(var.expl(calc.lms.main(dvs=c("se.true","mean.bias.est","rmse.est"), ivs=c("loadings","estimator","length","likelihood","blocksize"), 
                                results=res2[res2$trait.level==0,])))
#at +2
rm.0rows(var.expl(calc.lms.main(dvs=c("se.true","mean.bias.est","rmse.est"), ivs=c("loadings","estimator","length","likelihood","blocksize"), 
                                results=res2[res2$trait.level==2,])))

#means and SDs by condition
#only for conditions that showed differences, probably in interaction
means.info <- mean.frame(dvs=c("mean.bias","mean.ratio", "rmse"), ivs=c("estimator","likelihood","blocksize","method"), results=res.long.trait2, na.rm=TRUE)
mean.frame(dvs=c("mean.bias","rmse"), ivs=c("estimator","likelihood","blocksize"), results=res.long.trait2[res.long.trait2$loadings=="high" & res.long.trait2$length=="short",], na.rm=TRUE)

means.info.obs <- mean.frame(dvs=c("rmse","var.ratio"), 
                             ivs=c("estimator","likelihood","blocksize", "length"), 
                             results=res.long.trait2[res.long.trait2$method == "observed", ], na.rm=TRUE)


#---------------- textable variance explanation ------------#

#format for latex
var.paper <- var.expl.info.both
rownames(var.paper) <- gsub(":","$\\\\times$",rownames(var.paper))
rownames(var.paper) <- gsub("residuals","Residuals",rownames(var.paper))
colnames(var.paper) <- c("MB","MR","RMSE")
var.paper[is.na(var.paper)] <- "--"

header <- list()
header$pos <- list(-1, nrow(var.paper))
header$command <- c("\\hline \n Factor & MB & MR & RMSE\\\\",
                    "\\hline \n \\multicolumn{4}{l}{\\small \\textit{Note.} MB = Mean Bias, MR = Mean Ratio, RMSE = Root Mean} \\\\
                    \\multicolumn{4}{l}{\\small Squared Error. For the expected \\textit{SE}s, the bias and the ratio are} \\\\
                    \\multicolumn{4}{l}{\\small shown. The RMSE was only computed for observed \\textit{SE}s.} \n")

print(xtable::xtable(var.paper, digits=0, align = "lrrr", 
                     caption="Variance in bias for information-based standard errors explained in \\% by the manipulated factors in simulation study 1 on standard error accuracy",
                     label="tb:var_sim_se"), include.colnames = F, include.rownames=T, hline.after=c(-1,nrow(var.paper)-1),
      sanitize.rownames.function=function(x){x}, sanitize.colnames.function = function(x){x},
      sanitize.text.function = function(x){x},
      NA.string = "", table.placement = "htp", add.to.row = header,
      caption.placement = "top", latex.environments = NULL,
      file="../paper/Revision3_Psychometrika/manuscript_Psychometrika_Revision3/textable_var_se.tex")

#----------------- textable means ----------------#
means.paper <- means.info
#re-order columns
factors <- 1:4
nlevels <- apply(means.paper[,factors], 2, function(cl) length(unique(cl)))

#only every first occurence of factor level
for(cl in factors) {
  if(cl > 1) {
    means.paper[,cl] <- as.character(means.paper[,cl])
    means.paper[-c(seq(1, nrow(means.paper), prod(nlevels[1:(cl-1)]))), cl] <- ""
  }
}
means.paper[,factors] <- means.paper[,rev(factors)]
means.paper <- add.brackets.SDs(means.paper)
means.paper[is.na(means.paper)] <- ""

header <- list()
posmidlines <- seq(prod(nlevels[1:2]), nrow(means.paper)-prod(nlevels[1:2]), prod(nlevels[1:2]))
posmidlines <- posmidlines[! posmidlines %in% (nrow(means.paper)/2)]
header$pos <- list(-1) 
header$pos <- append(header$pos, as.list(posmidlines))
header$pos <- append(header$pos, nrow(means.paper))
header$command <- c("\\hline \n Method & Blocksize & Likelihood & Estimator & \\multicolumn{2}{c}{MB} & \\multicolumn{2}{c}{MR} 
                    & \\multicolumn{2}{c}{RMSE}\\\\
                    \\cmidrule(rl){1-4} \\cmidrule(rl){5-6} \\cmidrule(rl){7-8} \\cmidrule(rl){9-10}",
                    rep("\\cmidrule(rl){2-10}", length(posmidlines)),
                    "\\hline \n \\multicolumn{10}{l}{\\small \\textit{Note.} MB = Mean Bias, MR = Mean Ratio, RMSE = Root Mean Squared Error, 
                    ML = Maximum} \\\\
                    \\multicolumn{10}{l}{\\small Likelihood, MAP = Maximum a posteriori. Standard deviations are given in parentheses.} \n")

print(xtable::xtable(means.paper, digits=2, align="lllllrlrlrl",
                     caption="Means of bias for information-based standard errors by condition in simulation study 1 on standard error accuracy",
                     label="tb:var_means_se"), include.colnames = F, include.rownames=F, hline.after=c(nrow(means.paper)/2),
      sanitize.rownames.function=function(x){x}, sanitize.colnames.function = function(x){x},
      sanitize.text.function = function(x){x},
      NA.string = "", table.placement = "htp", add.to.row = header,
      caption.placement = "top", latex.environments = NULL,
      file="../paper/Revision3_Psychometrika/manuscript_Psychometrika_Revision3/textable_means_se.tex")

####------------------------ descriptive in text -----------------------####

#- bias of information-based SEs

#- by likelihood
mean.frame(dvs=c("mean.bias"), ivs=c("likelihood"), results=res.long.trait2, na.rm=T)

mean.frame(dvs=c("mean.bias"), ivs=c("likelihood", "blocksize"), results=res.long.trait2, na.rm=T)

mean.frame(dvs=c("mean.bias"), ivs=c("likelihood", "estimator"), results=res.long.trait2, na.rm=T)

mean.frame(dvs=c("mean.bias"), ivs=c("likelihood", "estimator", "blocksize"), results=res.long.trait2, na.rm=T)

# range MR
mean.frame(dvs=c("mean.ratio"), ivs=c("likelihood", "estimator"), results=res.long.trait2, na.rm=T,
           fn1 = min, fn2 = max)

# at theta = 0
mean.frame(dvs=c("mean.ratio"), ivs=c("likelihood", "estimator"), 
           results=res.long.trait2[res.long.trait2$trait.level == 0,], na.rm=T,
           fn1 = min, fn2 = max)
mean.frame(dvs=c("mean.ratio"), ivs=c("likelihood"), 
           results=res.long.trait2[res.long.trait2$trait.level == 0,], na.rm=T,
           fn1 = min, fn2 = max)

# at theta = +/- 2
mean.frame(dvs=c("mean.ratio"), ivs=c("likelihood", "estimator"), 
           results=res.long.trait2[res.long.trait2$trait.level %in% c(-2,2),], na.rm=T,
           fn1 = min, fn2 = max)

# RMSE
mean.frame(dvs=c("rmse"), ivs=c("likelihood", "estimator"), results=res.long.trait2, na.rm=T)

#- expected versus observed
#range MB info
round(mean(res.long.trait2$mean.bias, na.rm=T), 2)
round(sd(res.long.trait2$mean.bias, na.rm=T), 2)

# - bias of trait estimates (SOM)

#range RMSE by length
tapply(res2$rmse.est, res2$length, function(vec) round(range(vec), 2))

#MB for +2 by estimator
mean.frame(dvs=c("mean.bias.est"), ivs=c("estimator"), results=res2[res2$trait.level==2,])

mean.frame(dvs=c("se.true"), ivs=c("loadings","estimator","length"), results=res2)

#mean RMSE by blocksize
mean.frame(dvs=c("rmse.est"), ivs=c("blocksize"), results=res2, na.rm=T)

# - empirical SEs
#how often did the box constraint activate, for +-2, by estimator
table(res2[res2$trait.level %in% c(-2,2),]$n.box.2, res2[res2$trait.level %in% c(-2,2),]$estimator)
t(as.matrix(table(res2[res2$trait.level %in% c(-2,2),]$n.box.2 > 0, res2[res2$trait.level %in% c(-2,2),]$estimator)))/
  as.vector(table(res2[res2$trait.level %in% c(-2,2),]$estimator))
mean.frame(dvs=c("n.box.2"), ivs=c("estimator"), results=res2[res2$trait.level %in% c(-2,2),], na.rm=T, rn=0)


####------------------------ additional analyses -----------------------####

