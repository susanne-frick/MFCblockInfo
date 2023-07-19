####------------------- simulation testinfo -----------------------####
#for selected values of 5 traits
#compute SEs at true theta

# devtools::load_all()
# library(doParallel)

library(doMPI)
library(MFCblockInfo)
# source("~/Dokumente/Projekte/MFC_blocks/MFCblockInfo_simulations/functions_rel_emp_Diss.R")
source("functions_rel_emp_Diss.R")

####------------------- simulation design -------------------------####

# design.load.all <- readRDS("simulation/design_load_all_234.rds")
design.load.all <- readRDS("design_load_all_234.rds")

factor.blocksize <- c(2:4)
factor.loadings <- c("high","low")
factor.estimator <- c("MAP","ML")
factor.length <- c("short","long")
factor.likelihood <- c("true","dependent")

#trait levels (fixed)
# trait.levels <- create.grid(seq(-1,1,1), 5)
#fix 2nd trait (extraversion)
trait.levels <- matrix(0, length(seq(-2,2,.5)), 5)
trait.levels[,2] <- seq(-2,2,.5)
# trait.levels <- trait.levels[5:9,] # reduce for testing

Q <- 200 #500 #Number of questionnaires per condition

design.sim <- expand.grid("theta" = "true", 
                          "blocksize"=factor.blocksize, "loadings"=factor.loadings, "length"=factor.length,
                          "estimator"=factor.estimator, "likelihood"=factor.likelihood,
                          "p"=1:nrow(trait.levels), "q"=1:Q)

####-------------------- fixed conditions -------------------------####

#trait correlations (Big 5 from meta-analysis van der Linden et al.)
trait.cov <- matrix(c(1,-.36,-.17,-.36,-.43,
                      -.36,1,.43,.26,.29,
                      -.17,.43,1,.21,.20,
                      -.36,.26,.21,1,.43,
                      -.43,.29,.20,.43,1),
                    nrow=5,ncol=5)
load.range <- list("low"=c(.45,.75), "high"=c(.65,.95))

R <- 500 #number of responses per participant and questionnaire

####------------------ start simulation -------------------####

# cl <- makeCluster(10)
# registerDoParallel(cl)

cl <- startMPIcluster()
registerDoMPI(cl)

sinkWorkerOutput(paste0("worker_iter.out"))

# define chunkSize so that each cluster worker gets a single task chunk
chunkSize <- ceiling(nrow(design.sim)/getDoParWorkers())
mpiopts <- list(chunkSize=chunkSize)


res <- foreach (d=1:nrow(design.sim), .combine=rbind, .verbose=T, .packages=c("mvtnorm","numDeriv","devtools"),
                .inorder=F, .errorhandling="remove", .options.mpi=mpiopts) %dopar% {
                    set.seed(0805+d)

    design.load <- design.load.all[[as.character(design.sim[d,"blocksize"])]][["12"]]
    if(design.sim[d,"length"]=="long") design.load <- rbind(design.load, design.load)

    nb <- design.sim[d,"blocksize"]
    K <- nrow(design.load)/nb

    ####------------------ simulate item parameters -------------------####
    items <- sim.items(design.load=design.load, K=K, nb=nb,
                       load.range=load.range[[as.character(design.sim[d,"loadings"])]],
                       int.range=c(-1,1))

    ####-------------------------------- traits and responses --------------------------####

    traits <- matrix(rep(as.numeric(trait.levels[design.sim[d,"p"],]), R), R, ncol(design.load), byrow=TRUE)
    
    responses <- sim.responses(traits, items, design.load, K, nb, return.index=F)

    ####------------------- estimate traits ---------------------------####

    load.mat <- items$loads * design.load
    design.mat <- create.design.mat(K=K, nb=nb)
    gamma.true <- design.mat %*% items$u.mean

    if(design.sim[d,"likelihood"]=="true") {
      estimates.start <- est.MAP.old(FUN=posterior.tirt, responses=responses$outcomes, 
                                     design.mat=design.mat, comp.mat.a=design.mat %*% load.mat, int=items$u.mean, uni=items$uni,
                                     m.prior=rep(0,ncol(design.load)), s.prior=trait.cov, est=as.character(design.sim[d,"estimator"]))
      estimates <- est.MAP(FUN=lhb.mplus, responses=responses$rankindices,
                           int=gamma.true, loads=load.mat, uni=diag(items$uni),
                           perms=permute(1:nb), nb=nb,
                           m.prior=rep(0,ncol(design.load)), s.prior=trait.cov, est=as.character(design.sim[d,"estimator"]),
                           starts=estimates.start$traits)
      ses.obs <- est.hessian(FUN=lhb.mplus, responses=responses$rankindices, traits = traits,
                           int=gamma.true, loads=load.mat, uni=diag(items$uni),
                           perms=permute(1:nb), nb=nb,
                           m.prior=rep(0,ncol(design.load)), s.prior=trait.cov, est=as.character(design.sim[d,"estimator"]))$ses
    } else if(design.sim[d,"likelihood"]=="dependent") {
        estimates <- est.MAP.old(FUN=posterior.tirt, responses=responses$outcomes, 
                             design.mat=design.mat, comp.mat.a=design.mat %*% load.mat, int=items$u.mean, uni=items$uni,
                             m.prior=rep(0,ncol(design.load)), s.prior=trait.cov, est=as.character(design.sim[d,"estimator"]))
        ses.obs <- est.hessian.old(FUN=posterior.tirt, traits = traits, responses=responses$outcomes, 
                                 design.mat=design.mat, comp.mat.a=design.mat %*% load.mat, int=items$u.mean, uni=items$uni,
                                 m.prior=rep(0,ncol(design.load)), s.prior=trait.cov, est=as.character(design.sim[d,"estimator"]))$ses
    }
    n.err <- length(estimates$errors)
    n.warn <- length(estimates$warns)
    n.mess <- length(estimates$messages)
    
    #how often does the box constraint activate?
    #for the 2nd trait, out of R responses
    n.box.2 <- sum(estimates$traits[,2] %in% c(-3,3))
    #for the other traits (summed across traits -> 4*R responses)
    n.box.others <- sum(estimates$traits[,-c(2)] %in% c(-3,3))
    
    se.true <- apply(estimates$traits, 2, sd)
    
    #bias for estimates
    mean.bias.est <- colMeans( estimates$traits - traits)
    rmse.est <- sqrt(colMeans( (estimates$traits - traits)^2 ))

    ####------------------- SEs based on info --------------------------------####
    
    s.prior.d <- switch(as.character(design.sim[d,"estimator"]), "MAP"=trait.cov, "ML"=NULL)
    
    if(as.character(design.sim[d,"likelihood"])=="true") {
        
      info.blocks <- calc.info.block(lhb.mplus, traits=traits, 
                                     int=gamma.true, loads=load.mat, uni=diag(items$uni),
                                     K=K, nb=nb, responses = responses$rankindices)
      ses.info <- info2se(info.blocks$expected, summed=T, prior=s.prior.d)[1,, drop = FALSE] #here, the traits are all the same
    } else {
      
      info.pairs <- calc.info.pairs(traits=traits, responses = responses$outcomes, 
                                    int=gamma.true, loads=load.mat, uni=items$uni, design.mat=design.mat)
      ses.info <- info2se(info.pairs$expected, summed=T, prior=s.prior.d)[1,, drop = FALSE]
    }
    
    #number of persons with at least one NA
    na.info <- sum(is.na(ses.info))
    
    mean.bias.obs <- rowMeans( t(ses.obs) - se.true)
    rmse.obs <- sqrt(rowMeans( (t(ses.obs) - se.true)^2 ))
    mean.ratio.obs <- rowMeans( t(ses.obs) / se.true)
    var.ratio.obs <- apply(t(ses.obs) / se.true, 1, var)
    
    bias.exp <- as.vector(ses.info - se.true) # RMSE does not make sense here
    ratio.exp <- as.vector(ses.info/se.true)
    
    ####-------------------- save results -------------------------------------####
    res.r <- data.frame(design.sim[d,], "trait"=1:ncol(design.load), n.err, n.warn, n.mess, n.box.2, n.box.others, na.info,
                        se.true, 
                        mean.bias.obs, rmse.obs, mean.ratio.obs, var.ratio.obs, 
                        bias.exp, ratio.exp, 
                        mean.bias.est, rmse.est)
    saveRDS(res.r, file=paste0("results_dependent_SE-at-true/results_simulation_dependent_SE-at-true_d",d,".rds"))
    res.r
    #res[d,] <- res.r
}

# saveRDS(res, file="simulation/results_simulation_testinfo_grid.rds")
saveRDS(res, file="results_simulation_testinfo_grid_dependent_conditions_SE-at-true.rds")

# stopCluster(cl)
# colnames(res) <- resnames
closeCluster(cl)
mpi.quit()

