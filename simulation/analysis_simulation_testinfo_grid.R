####--------------- analysis simulation testinfo ---------------####
library(psych)
library(gridExtra)
library(ggplot2)
devtools::load_all()
devtools::load_all("../../DataAnalysisSimulation/")

#function to plot means by condition
ggplot.means <- function(res, dv, ivs, ivx, se.type="sd", main="Trait", N.obs=NULL, y_lab, size=16,...) {

  form <- formula(paste(dv, "~", paste(c(ivx, ivs), collapse=" + ")))

  means <- aggregate(form, FUN=mean, data=res)
  sds <- aggregate(form, FUN=sd, data=res)

  if(se.type=="confidence") {
    means$lower <- means[,dv] - qnorm(.975)*sds[,dv]/N.obs
    means$upper <- means[,dv] + qnorm(.975)*sds[,dv]/N.obs
  } else if(se.type=="sd") {
    means$lower <- means[,dv] - sds[,dv]
    means$upper <- means[,dv] + sds[,dv]
  } else if(se.type=="quantile") {
    means$lower <- aggregate(form, FUN=function(r) quantile(r, .05), data=res)[,length(c(ivx,ivs))+1]
    means$upper <- aggregate(form, FUN=function(r) quantile(r, .95), data=res)[,length(c(ivx,ivs))+1]
  }

  ggplot(data=means, aes(x=get(ivx), y=get(dv), ymin=lower, ymax=upper,
                         group=interaction(means[,ivs]), fill=estimator, linetype=loadings, ...)) +
    geom_line() +
    geom_ribbon(alpha=0.5) +
    ggtitle(main) +
    xlab(as.expression(expression(theta))) +
    ylab(y_lab) +
    theme(axis.text=element_text(size=size),
          axis.title=element_text(size=size),
          legend.text=element_text(size=size),
          legend.title=element_text(size=size))
}

####--------------- read in and check results ------------------####
res <- readRDS("simulation/results_simulation_testinfo_grid.rds")
head(res)

#check NAs = info that could not be estimated
anyNA(res)
table(res$na.info) #total number of NAs
table(res$na.info>0)/nrow(res) #percentage
#! much more than before
describe(res)

####---------------- reformat results --------------------------####
#recode loadings
res$loadings <- car::recode(res$loadings, " 'good'='high'; 'bad'='low' ")

trait.seq <- seq(-2,2,.5) #trait levels that were varied
#add indicator for type of person: all other traits 0 or equal
res$equal <- ifelse(res$p %in% 1:length(trait.seq), 0, 1)
#add level of Trait 2 / all traits
res$trait.level <- recode.df(res$p, 1:(length(trait.seq)*2), rep(trait.seq, 2))

#remove blocksize (was not varied)
table(res$blocksize)
res <- res[,grep("blocksize",colnames(res),invert=T)]

#separate persons with all other traits 0 from all other traits equal
res2 <- res[res$equal==0 & res$trait==2,]

#to long format: observed (at estimate) vs. expected (weighted with probability)
fixed <- c("loadings","estimator","length","p","trait.level","equal","q","trait")
res.long <- rbind(res[,c(fixed,"mean.bias.obs","rmse.obs")],
                  setNames(res[,c(fixed,"mean.bias.exp","rmse.exp")], c(fixed,"mean.bias.obs","rmse.obs")))
colnames(res.long) <- c(fixed,"mean.bias","rmse")
#add indicator for method: expected vs. observed
res.long$method <- rep(c("observed","expected"), each=nrow(res))
#re-order columns
res.long <- res.long[,c(fixed,"method","mean.bias","rmse")]
head(res.long)
tail(res.long)

#separate persons with all other traits 0 from all other traits equal
res.long.trait2 <- res.long[res.long$equal==0 & res.long$trait==2,]

####------------------------------- plots -----------------------------------------####

plot.ses.short <- ggplot.means(res2[res2$length=="short",], dv="se.true", ivx="trait.level", ivs=c("loadings","estimator"), main="", y_lab="SE")
plot.mab.short <- ggplot.means(res2[res2$length=="short",], dv="mean.bias.est", ivx="trait.level", ivs=c("loadings","estimator"), main="", y_lab="MB")
plot.rmse.short <- ggplot.means(res2[res2$length=="short",], dv="rmse.est", ivx="trait.level", ivs=c("loadings","estimator"), main="", y_lab="RMSE")
plot.ses.long <- ggplot.means(res2[res2$length=="long",], dv="se.true", ivx="trait.level", ivs=c("loadings","estimator"), main="", y_lab="SE")
plot.mab.long <- ggplot.means(res2[res2$length=="long",], dv="mean.bias.est", ivx="trait.level", ivs=c("loadings","estimator"), main="", y_lab="MB")
plot.rmse.long <- ggplot.means(res2[res2$length=="long",], dv="rmse.est", ivx="trait.level", ivs=c("loadings","estimator"), main="", y_lab="RMSE")
ggsave("simulation/plot_recovery.pdf",
       plot=lemon::grid_arrange_shared_legend(plot.ses.short, plot.mab.short, plot.rmse.short,
                                              plot.ses.long, plot.mab.long, plot.rmse.long,
                                              nrow=2, ncol=3, position="right"),
       width=10, height=8, units="in")
ggsave("../paper/plot_recovery.pdf",
       plot=lemon::grid_arrange_shared_legend(plot.ses.short, plot.mab.short, plot.rmse.short,
                                              plot.ses.long, plot.mab.long, plot.rmse.long,
                                              nrow=2, ncol=3, position="right"),
       width=10, height=8, units="in")


#bias observed
plot.mab <- ggplot.means(res2[res2$length=="short",], dv="mean.bias.obs", ivx="trait.level", ivs=c("loadings","estimator"), main="", y_lab="MB")
plot.rmse <- ggplot.means(res2[res2$length=="short",], dv="rmse.obs", ivx="trait.level", ivs=c("loadings","estimator"), main="", y_lab="RMSE")
plot.mab.long <- ggplot.means(res2[res2$length=="long",], dv="mean.bias.obs", ivx="trait.level", ivs=c("loadings","estimator"), main="", y_lab="MB")
plot.rmse.long <- ggplot.means(res2[res2$length=="long",], dv="rmse.obs", ivx="trait.level", ivs=c("loadings","estimator"), main="", y_lab="RMSE")
ggsave("simulation/plot_SEs_observed.pdf",
       plot=lemon::grid_arrange_shared_legend(plot.mab, plot.rmse, plot.mab.long, plot.rmse.long,
                                              nrow=2, ncol=2, position="right"),
       width=10, height=8, units="in")
ggsave("../paper/plot_SEs_observed.pdf",
       plot=lemon::grid_arrange_shared_legend(plot.mab, plot.rmse, plot.mab.long, plot.rmse.long,
                                              nrow=2, ncol=2, position="right"),
       width=10, height=8, units="in")

#bias expected
plot.mab <- ggplot.means(res[res$equal==0 & res$trait==2 & res$length=="short",], dv="mean.bias.exp", ivx="trait.level", ivs=c("loadings","estimator"), main="", y_lab="MB")
plot.rmse <- ggplot.means(res[res$equal==0 & res$trait==2 & res$length=="short",], dv="rmse.exp", ivx="trait.level", ivs=c("loadings","estimator"), main="", y_lab="RMSE")
plot.mab.long <- ggplot.means(res[res$equal==0 & res$trait==2 & res$length=="long",], dv="mean.bias.exp", ivx="trait.level", ivs=c("loadings","estimator"), main="", y_lab="MB")
plot.rmse.long <- ggplot.means(res[res$equal==0 & res$trait==2 & res$length=="long",], dv="rmse.exp", ivx="trait.level", ivs=c("loadings","estimator"), main="", y_lab="RMSE")
ggsave("simulation/plot_SEs_expected.pdf",
       plot=lemon::grid_arrange_shared_legend(plot.mab, plot.rmse, plot.mab.long, plot.rmse.long,
                                              nrow=2, ncol=2, position="right"),
       width=10, height=8, units="in")
ggsave("../paper/plot_SEs_expected.pdf",
       plot=lemon::grid_arrange_shared_legend(plot.mab, plot.rmse, plot.mab.long, plot.rmse.long,
                                              nrow=2, ncol=2, position="right"),
       width=10, height=8, units="in")

####----------------- ANOVA --------------------####

#observed vs. expected
lm.info.main <- calc.lms.main(dvs=c("mean.bias","rmse"), ivs=c("loadings","estimator","length","method"), results=res.long.trait2)
var.expl.info <- rm.0rows(var.expl(lm.info.main))

#means and SDs by condition
means.info <- mean.frame(dvs=c("mean.bias","rmse"), ivs=c("loadings","estimator","length","method"), results=res.long.trait2, na.rm=TRUE)

#---------------- textable variance explanation ------------#

#format for latex
var.paper <- var.expl.info
rownames(var.paper) <- gsub(":","$\\\\times$",rownames(var.paper))
rownames(var.paper) <- gsub("residuals","Residuals",rownames(var.paper))
colnames(var.paper) <- c("MB","RMSE")

header <- list()
header$pos <- list(-1, nrow(var.paper))
header$command <- c("\\hline \n Factor & MB & RMSE\\\\",
                    "\\hline \n \\multicolumn{3}{l}{\\small \\textit{Note.} MB = Mean Bias, RMSE = Root Mean} \\\\
                    \\multicolumn{3}{l}{\\small Squared Error, ML = Maximum Likelihood,} \\\\
                    \\multicolumn{3}{l}{\\small MAP = Maximum a posteriori. Expected vs. } \\\\
                    \\multicolumn{3}{l}{\\small observed explained less  than 1\\% of variance.} \n")

print(xtable::xtable(var.paper, digits=0,
                     caption="Variance in bias for information-based standard errors explained in \\% by the manipulated factors in simulation study 1 on standard error accuracy",
                     label="tb:var_sim_se"), include.colnames = F, include.rownames=T, hline.after=c(-1,nrow(var.paper)-1),
      sanitize.rownames.function=function(x){x}, sanitize.colnames.function = function(x){x},
      sanitize.text.function = function(x){x},
      NA.string = "", table.placement = "htp", add.to.row = header,
      caption.placement = "top", latex.environments = NULL,
      file="../paper/textable_var_se.tex")

#----------------- textable means ----------------#
means.paper <- means.info
#re-order columns
means.paper[,1:4] <- means.paper[,4:1]
#only every first occurence of factor level
for(cl in 1:4) {
  means.paper[,cl] <- as.character(means.paper[,cl])
  means.paper[-c(seq(1, nrow(means.paper), nrow(means.paper)/2^cl)), cl] <- ""
}
means.paper <- add.brackets.SDs(means.paper)

header <- list()
header$pos <- list(-1, nrow(means.paper)/4, 3*nrow(means.paper)/4, nrow(means.paper))
header$command <- c("\\hline \n Method & Length & Estimator & Loadings & \\multicolumn{2}{c}{MB} & \\multicolumn{2}{c}{RMSE} \\\\
                    \\cmidrule(rl){1-4} \\cmidrule(rl){5-6} \\cmidrule(rl){7-8}",
                    "\\cmidrule(rl){2-8}",
                    "\\cmidrule(rl){2-8}",
                    "\\hline \n \\multicolumn{8}{l}{\\small \\textit{Note.} MB = Mean Bias, RMSE = Root Mean Squared Error, ML = Maximum} \\\\
                    \\multicolumn{8}{l}{\\small Likelihood, MAP = Maximum a posteriori. Standard deviations are given in} \\\\
                    \\multicolumn{8}{l}{\\small brackets.} \n")

print(xtable::xtable(means.paper, digits=2, align="lllllrlrl",
                     caption="Means in bias for information-based standard errors by condition in simulation study 1 on standard error accuracy",
                     label="tb:var_means_se"), include.colnames = F, include.rownames=F, hline.after=c(nrow(means.paper)/2),
      sanitize.rownames.function=function(x){x}, sanitize.colnames.function = function(x){x},
      sanitize.text.function = function(x){x},
      NA.string = "", table.placement = "htp", add.to.row = header,
      caption.placement = "top", latex.environments = NULL,
      file="../paper/textable_means_se.tex")

####------------------------ descriptive in text -----------------------####

#range RMSE by length
tapply(res2$rmse.est, res2$length, function(vec) round(range(vec), 2))

#MB for +2 by estimator
mean.frame(dvs=c("mean.bias.est"), ivs=c("estimator"), results=res2[res2$trait.level==2,])

mean.frame(dvs=c("se.true"), ivs=c("loadings","estimator","length"), results=res2)

#range MB info
round(mean(res.long.trait2$mean.bias, na.rm=T), 2)
round(sd(res.long.trait2$mean.bias, na.rm=T), 2)

#mean MB by estimator
mean.frame(dvs=c("mean.bias"), ivs=c("estimator"), results=res.long.trait2, na.rm=T)

####------------------------ additional analyses -----------------------####

#info-based SEs
res$mean.se.obs <- res$mean.bias.obs + res$se.true
res$mean.se.exp <- res$mean.bias.exp + res$se.true
plot.se.obs <- ggplot.means(res[res$equal==0 & res$trait==2 & res$length=="short",], dv="mean.se.obs", ivx="trait.level",
                            ivs=c("loadings","estimator"), main="", y_lab="Mean SE")
plot.se.exp <- ggplot.means(res[res$equal==0 & res$trait==2 & res$length=="short",], dv="mean.se.exp", ivx="trait.level",
                            ivs=c("loadings","estimator"), main="", y_lab="Mean SE")
plot.combined <- lemon::grid_arrange_shared_legend(plot.se.obs, plot.se.exp,
                                                   nrow=1, ncol=2, position="right")

#all traits equal
plot.ses.equal.short <- ggplot.means(res[res$equal==1 & res$length=="short",], main="", dv="se.true", ivx="trait.level", ivs=c("loadings","estimator"), y_lab="SE")
plot.mb.equal.short <- ggplot.means(res[res$equal==1 & res$length=="short",], main="", dv="mean.bias.est", ivx="trait.level", ivs=c("loadings","estimator"), y_lab="SE")
plot.rmse.equal.short <- ggplot.means(res[res$equal==1 & res$length=="short",], main="", dv="rmse.est", ivx="trait.level", ivs=c("loadings","estimator"), y_lab="SE")
plot.ses.equal.long <- ggplot.means(res[res$equal==1 & res$length=="long",], main="", dv="se.true", ivx="trait.level", ivs=c("loadings","estimator"), y_lab="SE")
plot.mb.equal.long <- ggplot.means(res[res$equal==1 & res$length=="long",], main="", dv="mean.bias.est", ivx="trait.level", ivs=c("loadings","estimator"), y_lab="SE")
plot.rmse.equal.long <- ggplot.means(res[res$equal==1 & res$length=="long",], main="", dv="rmse.est", ivx="trait.level", ivs=c("loadings","estimator"), y_lab="SE")
plot.combined.equal <- lemon::grid_arrange_shared_legend(plot.ses.equal.short, plot.mb.equal.short, plot.rmse.equal.short,
                                                         plot.ses.equal.long, plot.mb.equal.long, plot.rmse.equal.long,
                                                         nrow=2, ncol=3, position="right")
ggsave("simulation/plot_recovery_equal.pdf", plot=plot.combined.equal,
       width=10, height=8, units="in")

#ANOVA bias SEs, separately for trait levels
#theta = 2
lm.info.main <- calc.lms.main(dvs=c("mean.bias","rmse"), ivs=c("loadings","estimator","length","method"),
                              results=res.long.trait2[res.long.trait2$trait.level==2,])
var.expl(lm.info.main)
#theta = 0
lm.info.main <- calc.lms.main(dvs=c("mean.bias","rmse"), ivs=c("loadings","estimator","length","method"),
                              results=res.long.trait2[res.long.trait2$trait.level==0,])
var.expl(lm.info.main)

#contrasts
#! currently: only rows where information was estimated, adapt calc.mm instead?
lm.info <- calc.lms(dvs=c("mean.bias","rmse"), ivs=c("loadings","estimator","length","method"), results=res.long.trait2[is.na(res.long.trait2$mean.bias)==FALSE,])
var.expl(lm.info)

#means and SDs across observed vs. expected
mean.frame(dvs=c("mean.bias","rmse"), ivs=c("loadings","estimator"), results=res.long.trait2, na.rm=TRUE)

#ANOVA SE true, MAB, RMSE
lm.est.main <- calc.lms.main(dvs=c("se.true","mean.bias.est","rmse.est"), ivs=c("loadings","estimator","length"), results=res[res$equal==0 & res$trait==2,])
var.expl(lm.est.main)
lm.est <- calc.lms(dvs=c("se.true","mean.bias.est","rmse.est"), ivs=c("loadings","estimator","length"), results=res[res$equal==0 & res$trait==2,])
var.expl(lm.est)
#separately for trait levels
#theta = 2
lm.est.main <- calc.lms.main(dvs=c("se.true","mean.bias.est","rmse.est"), ivs=c("loadings","estimator","length"), results=res[res$equal==0 & res$trait==2
                                                                                                                              & res$trait.level==2,])
var.expl(lm.est.main)
#theta = 0
lm.est.main <- calc.lms.main(dvs=c("se.true","mean.bias.est","rmse.est"), ivs=c("loadings","estimator","length"), results=res[res$equal==0 & res$trait==2
                                                                                                                              & res$trait.level==0,])
var.expl(lm.est.main)
#means and SDs
mean.frame(dvs=c("se.true","mean.bias.est","rmse.est"), ivs=c("loadings","estimator","length"), results=res[res$equal==0 & res$trait==2,])


####------------------------ base r plot function -----------------####

#plot by trait levels
plot.means <- function(res, ...) {
  ses.agg <- aggregate(se.true ~ p + estimator + loadings, FUN=mean, data=res)
  ses.agg$tr.levels <- trait.seq
  plot(trait.seq, ses.agg$se.true[ses.agg$estimator=="MAP" & ses.agg$loadings=="good"], type="l", xlab=expression(theta), ylab="SE",
       ylim=range(ses.agg$se.true), ...)
  lines(trait.seq, ses.agg$se.true[ses.agg$estimator=="MAP" & ses.agg$loadings=="bad"], lty="dashed")
  lines(trait.seq, ses.agg$se.true[ses.agg$estimator=="ML" & ses.agg$loadings=="good"], lty="solid", col="darkgrey")
  lines(trait.seq, ses.agg$se.true[ses.agg$estimator=="ML" & ses.agg$loadings=="bad"], lty="dashed", col="darkgrey")
  legend(x=1, y=.45, title="estimator", legend=c("ML","MAP"), fill=c("black","darkgrey"))
  legend(x=1, y=.55, title="loadings", legend=c("bad","good"), lty=c("solid","dashed"))
}
