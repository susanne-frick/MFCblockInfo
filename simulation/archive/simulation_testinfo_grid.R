####------------------- simulation testinfo -----------------------####
#for selected values of 5 traits

# devtools::load_all()
# library(doParallel)

library(doMPI)
library(MFCblockInfo)

####------------------- simulation design -------------------------####

# design.load.all <- readRDS("simulation/design_load_all_234.rds")
design.load.all <- readRDS("design_load_all_234.rds")

factor.blocksize <- c(2,3,4)
factor.loadings <- c("good","bad")
factor.estimator <- c("ML","MAP")
factor.length <- c("short","long")

#trait levels (fixed)
# trait.levels <- create.grid(seq(-1,1,1), 5)
#fix 2nd trait (extraversion)
trait.levels <- matrix(0, length(seq(-2,2,.5))*2, 5)
trait.levels[,2] <- seq(-2,2,.5)

Q <- 500 #Number of questionnaires per condition

design.sim <- expand.grid("blocksize"=factor.blocksize, "loadings"=factor.loadings, "length"=factor.length, "estimator"=factor.estimator,
                          "p"=1:nrow(trait.levels), "q"=1:Q)

####-------------------- fixed conditions -------------------------####

#trait correlations (Big 5 from meta-analysis van der Linden et al.)
trait.cov <- matrix(c(1,-.36,-.17,-.36,-.43,
                      -.36,1,.43,.26,.29,
                      -.17,.43,1,.21,.20,
                      -.36,.26,.21,1,.43,
                      -.43,.29,.20,.43,1),
                    nrow=5,ncol=5)
load.range <- list("bad"=c(.45,.75), "good"=c(.65,.95))

R <- 500 #number of responses per participant and questionnaire

####------------------ start simulation -------------------####

# set.seed(2407)
#res <- matrix(NA, nrow(design.sim), length(resnames))
#colnames(res) <- resnames

# cl <- makeCluster(10)
# registerDoParallel(cl)

cl <- startMPIcluster()
registerDoMPI(cl)

sinkWorkerOutput(paste0("worker_iter.out"))

# define chunkSize so that each cluster worker gets a single task chunk
chunkSize <- ceiling(R/getDoParWorkers())
mpiopts <- list(chunkSize=chunkSize)

# res <- foreach(d=1:nrow(design.sim), .packages=c("mvtnorm","numDeriv","devtools"), .combine=rbind) %dopar% {
#

res <- foreach (d=1:nrow(design.sim), .combine=rbind, .verbose=T, .packages=c("mvtnorm","numDeriv","devtools"),
                .inorder=F, .errorhandling="remove", .options.mpi=mpiopts) %dopar% {
                    #for (r in 2:R) {
                    set.seed(1204+d)

    # load_all()
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
    #traits <- mvtnorm::rmvnorm(n=J, mean=rep(0,ncol(design.load)), sigma = trait.cov, method="chol")

    responses <- sim.responses(traits, items, design.load, K, nb, return.index=F)

    ####------------------- estimate traits ---------------------------####

    load.mat <- items$loads * design.load
    gamma.true <- create.design.mat(K=K, nb=nb) %*% items$u.mean

    estimates <- est.MAP(FUN=lhb.mplus, responses=responses$rankindices,
                              int=gamma.true, loads=load.mat, uni=diag(items$uni),
                              perms=permute(1:nb), nb=nb,
                              m.prior=rep(0,ncol(design.load)), s.prior=trait.cov, est=as.character(design.sim[d,"estimator"]))

    se.true <- apply(estimates$traits, 2, sd)
    mean.bias.obs <- rowMeans( t(estimates$ses) - se.true)
    rmse.obs <- sqrt(rowMeans( (t(estimates$ses) - se.true)^2 ))

    #bias for estimates
    mean.bias.est <- colMeans( estimates$traits - traits)
    rmse.est <- sqrt(colMeans( (estimates$traits - traits)^2 ))

    ####------------------- SEs based on info --------------------------------####

    s.prior.d <- switch(as.character(design.sim[d,"estimator"]), "MAP"=trait.cov, "ML"=NULL)
    ses.info <- info2se(calc.info.block(lhb.mplus, traits=estimates$traits, int=gamma.true, loads=load.mat, uni=diag(items$uni),
                        K=K, nb=nb), summed=T, s.prior=s.prior.d)
    #number of persons with at least one NA
    na.info <- sum(rowSums(is.na(ses.info))>0)

    mean.bias.exp <- rowMeans( t(ses.info) - se.true)
    rmse.exp <- sqrt(rowMeans( (t(ses.info) - se.true)^2 ))

    ####-------------------- save results -------------------------------------####
    res.r <- data.frame(design.sim[d,], "trait"=1:ncol(design.load), na.info, se.true, mean.bias.obs, mean.bias.exp, rmse.obs, rmse.exp, mean.bias.est, rmse.est)
    res.r
    #res[d,] <- res.r
}

# saveRDS(res, file="simulation/results_simulation_testinfo_grid.rds")
saveRDS(res, file="results_simulation_testinfo_grid.rds")

# stopCluster(cl)
# colnames(res) <- resnames
closeCluster(cl)
mpi.quit()

