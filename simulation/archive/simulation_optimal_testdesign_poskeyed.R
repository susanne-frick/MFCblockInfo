####------------------- simulation testinfo -----------------------####
# devtools::load_all()
# library(doParallel)

library(doMPI)
library(MFCblockInfo)

####------------------- simulation design -------------------------####

# design.load.all <- readRDS("simulation/design_load_all_5-15_equal.rds")
design.load.all <- readRDS("design_load_all_5-15_equal.rds")

factor.blocksize <- 2:4
factor.keying <- "0" # c("0","12","23")
factor.int <- "large"
factor.load <- "acceptable"
factor.length <- "long"
factor.algorithm <- c("opt","r2","loads","random")
factor.constraints <- c("unconstrained","constrained")
factor.target <- c("weighted","equal")
factor.ntraits <- "5" #c("5","15")

#number of replications
R <- 200#500

design.sim <- expand.grid("blocksize"=factor.blocksize, "keying"=factor.keying, "length"=factor.length, "intercepts"=factor.int,
                          "loads"=factor.load, 
                          "constraints"=factor.constraints, "target"=factor.target, "ntraits"=factor.ntraits,
                          "rep"=1:R)

####-------------------- fixed conditions -------------------------####

#trait correlations (Big 5 from meta-analysis van der Linden et al.)
trait.cov.5 <- matrix(c(1,-.36,-.17,-.36,-.43,
                        -.36,1,.43,.26,.29,
                        -.17,.43,1,.21,.20,
                        -.36,.26,.21,1,.43,
                        -.43,.29,.20,.43,1),
                      nrow=5,ncol=5)

#trait.cov.m.15: randomly drawn such that 59% are small, 40% medium and 1% negligible
#covariance matrix for inverse Wishart: all covariances set to .3
cov15 <- diag(15)
cov15[cov15==0] <- .3
set.seed(515)
#draw from inverse Wishart with df=100
trait.cov.15 <- cov2cor(SimDesign::rinvWishart(1, 100, cov15))
#reverse traits 1, 6, 11
trait.cov.15[,1][trait.cov.15[,1]!=1] <- trait.cov.15[,1][trait.cov.15[,1]!=1]*-1
trait.cov.15[1,][trait.cov.15[1,]!=1] <- trait.cov.15[1,][trait.cov.15[1,]!=1]*-1
trait.cov.15[,6][trait.cov.15[,6]!=1] <- trait.cov.15[,6][trait.cov.15[,6]!=1]*-1
trait.cov.15[6,][trait.cov.15[6,]!=1] <- trait.cov.15[6,][trait.cov.15[6,]!=1]*-1
trait.cov.15[,11][trait.cov.15[,11]!=1] <- trait.cov.15[,11][trait.cov.15[,11]!=1]*-1
trait.cov.15[11,][trait.cov.15[11,]!=1] <- trait.cov.15[11,][trait.cov.15[11,]!=1]*-1

#combine both trait covariance matrices to a list
trait.cov.list <- list("5"=trait.cov.5, "15"=trait.cov.15)

int.range <- list("small"=c(-1,1), "large"=c(-2,2))
load.range <- list("good"=c(.65,.95), "acceptable"=c(.45,.95))

#create grid of traits if traits are not given
tr.levels <- c(-1,0,1)
tr.list <- vector("list", ncol(trait.cov.5))
for(tr in 1:length(tr.list)) tr.list[[tr]] <- tr.levels
traits.grid.5 <- expand.grid(tr.list)
#for 15 traits this is not feasible -> 3^15 = 4348907 trait levels
#draw random sample of 500 instead
set.seed(515)
traits.grid.15 <- mvtnorm::rmvnorm(n=500, sigma=diag(15)) # uncorrelated
#add all 0 as reference point
if(isFALSE(any(rowSums(traits.grid.15==0)==15))) traits.grid.15 <- rbind(traits.grid.15, rep(0,15))
traits.grid.list <- list("5"=traits.grid.5, "15"=traits.grid.15)

#reduce for testing
#traits.grid.list <- lapply(traits.grid.list, function(tl) rbind(tl[1:5,], rep(0, ncol(tl))))

J <-  500 #Number of participants

####------------------ start simulation -------------------####

# cl <- makeCluster(4)
# registerDoParallel(cl)

cl <- startMPIcluster()
registerDoMPI(cl)

sinkWorkerOutput(paste0("worker_iter_opt.out"))

# define chunkSize so that each cluster worker gets a single task chunk
chunkSize <- ceiling(R/getDoParWorkers())
mpiopts <- list(chunkSize=chunkSize)

# res <- foreach(d=1:nrow(design.sim), .packages=c("mvtnorm","numDeriv","devtools"), .combine=rbind) %dopar% {
#

res <- foreach (d=1:nrow(design.sim), .combine=rbind, .verbose=T, .packages=c("mvtnorm","numDeriv","devtools","lpSolveAPI","MFCblockInfo"),
                .inorder=F, .errorhandling="remove", .options.mpi=mpiopts) %dopar% {
                  set.seed(1204+d)
  
                  ####------------------- test design -----------------------####
                  #select design load according to test design condition
                  design.load <- design.load.all[[as.character(design.sim[d,"ntraits"])]][[as.character(design.sim[d,"blocksize"])]][["12"]]
                  #quadruple
                  design.load <- rbind(design.load, design.load, design.load, design.load)
                                  
                  nb <- design.sim[d,"blocksize"]
                  K <- nrow(design.load)/nb
                  blocks <- create.block.ind(K, nb)
                  
                  #specifications for lp model
                  #final test length
                  K.final <- K/4
                  
                  ####------------------ simulate item parameters -------------------####
  
                  items <- sim.items(design.load=design.load, K=K, nb=nb,
                                     load.range=load.range[[as.character(design.sim[d,"loads"])]],
                                     int.range=int.range[[as.character(design.sim[d,"intercepts"])]])
                  
                  ####-------------------------------- traits and responses --------------------------####
                  trait.cov <- trait.cov.list[[as.character(design.sim[d,"ntraits"])]]
                  traits.grid <- traits.grid.list[[as.character(design.sim[d,"ntraits"])]]
                  
                  if(as.character(design.sim[d,"target"])=="equal") {
                    #grid as traits -> to evaluate equal weights
                    traits <- rbind(traits.grid, traits.grid)
                    #weights equal
                    weights.grid <- a.all <- rep(1,nrow(traits.grid))
                  } else {
                    traits <- rmvnorm(n=J, mean=rep(0,ncol(design.load)), sigma = trait.cov, method="chol")
                    #weights for opt
                    weights.grid <- NULL #(implemented in function)
                    }
                  
                  responses <- sim.responses(traits, items, design.load, K, nb, return.index=F)
                  
                  ####------------------- constraints -----------------------------------####
                  
                  if(as.character(design.sim[d,"constraints"])=="constrained") {
                    #constraints on items per trait
                    traits.blocks <- create.traits.blocks(loads=design.load, which.blocks=1:K, nb=nb)
                    traits.blocks.ind <- do.call(cbind, (lapply(1:ncol(design.load), function(f, tb) apply(tb, 1, function(rw) ifelse(f %in% rw, 1, 0)), tb=traits.blocks)))
                    n.traits <- rep(K.final/ncol(design.load)*nb, ncol(design.load))
                    #constraints on item keying (comparisons between opposite-keyed items)
                    #at least 3/4 of blocks are mixed keyed
                    loads.blocks <- t(apply(blocks, 1, function(b, dl) colSums(dl[b,]), dl=design.load))
                    block.mixed <- ifelse(rowSums(loads.blocks)==nb, 0, 1)
                    n.mixed <- 3/4*K.final
                    #at least 1 negatively keyed item per trait
                    traits.neg.ind <- apply(loads.blocks, 2, function(rw) ifelse(rw==-1, 1, 0))
                    n.neg <- rep(1, ncol(design.load))
                    
                    #combine to constraint.list
                    constraint.list <- list("left"=cbind(traits.blocks.ind, block.mixed, traits.neg.ind),
                                            "operator"=c(rep("=",ncol(traits.blocks.ind)), ">=", rep(">=",ncol(traits.neg.ind))),
                                            "right"=c(n.traits, n.mixed, n.neg))
                    
                  } else {
                    constraint.list <- NULL
                  }
                  
                  ####------------------- T-optimality --------------------------------####
                  
                  load.mat <- items$loads * design.load
                  gamma.true <- create.design.mat(K=K, nb=nb) %*% items$u.mean
                  
                  infos <- calc.info.block(lhb.mplus, traits=as.matrix(traits.grid), int=gamma.true, loads=load.mat, uni=diag(items$uni),
                                           K=K, nb=nb)
                  #trace for each block (and grid point)
                  info.trace <- do.call(rbind, lapply(infos, function(ip) apply(ip, 1, function(i) sum(diag(i)))))
                  
                  ####------------ MIP algorithm with T-optimality -----------------------####
                  results.opt <- select.optimal(info.sum=info.trace, traits.grid=traits.grid, K=K, K.final=K.final,
                                                constraint.list=constraint.list, weights.grid=weights.grid)
                  
                  ####---------------- block selection based on R^2 ----------------------####
                  
                  if(as.character(design.sim[d,"target"])=="weighted") {
                    #weights for R^2
                    a.all <- rowSums(info2se(infos, var.out=T), na.rm=T)
                    a.all <- a.all/a.all[which(rowSums(traits.grid==0)==ncol(traits.grid))]
                  }
                  
                  #R^2 mean across persons and traits, weighted
                  means.r2 <- do.call(c, lapply(1:K, function(k,i) mean(rowMeans(calc.info.block.r2(i, wo.blocks=k))*a.all), i=infos))
                  means.r2 <- matrix(means.r2, 1, K)
                  results.r2 <- select.optimal(info.sum=means.r2, traits.grid=matrix(0,1,1), K=K, K.final=K.final,
                                               constraint.list=constraint.list)
                  #ind.r2 <- order(means.r2, decreasing = T)[1:K.final] #old version without constraints
                  
                  ####---------------- item selection based on loadings ------------------####
                  loads.blocks <- t(apply(blocks, 1, function(b, dl) colSums(dl[b,]), dl=load.mat))
                  means.loads <- matrix(rowMeans(abs(loads.blocks)), 1, K)
                  results.loads <- select.optimal(info.sum=means.loads, traits.grid=matrix(0,1,1), K=K, K.final=K.final,
                                                  constraint.list=constraint.list)
                  
                  #ind.loads <- order(means.loads, decreasing = T)[1:K.final] #old version without constraints
                  
                  ####---------------- random item selection ------------------####
                  results.rand <- select.optimal(info.sum=matrix(runif(K, 0, 1), 1, K), traits.grid=matrix(0,1,1), K=K, K.final=K.final,
                                                 constraint.list=constraint.list)
                  
                  ####---------------- trait estimation and summary measures -------------####
                  results.list <- list("opt"=results.opt, "r2"=results.r2, "loads"=results.loads, "random"=results.rand)
                  res.r <- NULL
                  for (a in factor.algorithm) {
                    if(results.list[[as.character(a)]]$solved==0) {
                      blocks.ind <- c(t(blocks[results.list[[as.character(a)]]$ind.opt,]))
                      gamma.true <- create.design.mat(K=K.final, nb=nb) %*% items$u.mean[blocks.ind]
                      
                      #estimate traits based on new questionnaire
                      estimates <- est.MAP(FUN=lhb.mplus, responses=responses$rankindices[,results.list[[as.character(a)]]$ind.opt],
                                           int=gamma.true, loads=load.mat[blocks.ind,], uni=diag(items$uni)[blocks.ind,blocks.ind],
                                           perms=permute(1:nb), nb=nb,
                                           m.prior=rep(0,ncol(design.load)), s.prior=trait.cov, SE=FALSE)
                      
                      rec <- diag(cor(estimates$traits, traits))
                      RMSE <- colMeans((estimates$traits - traits)^2)
                      MAB <- colMeans(abs(estimates$traits - traits))
                      
                      res.r <- rbind(res.r,
                                     data.frame(design.sim[d,], "trait"=1:ncol(design.load), "algorithm"=a, rec, RMSE, MAB))
                      
                    } else {
                      res.r <- data.frame(design.sim[d,], "trait"=1:ncol(design.load), "algorithm"="opt", rec=NA, RMSE=NA, MAB=NA)
                    }
                  }
                  saveRDS(res.r, file=paste0("results_opt_poskeyed/results_simulation_opt_poskeyed_d",d,".rds"))
                  res.r
                }

saveRDS(res, file="results_simulation_opt_poskeyed.rds")

# stopCluster(cl)
closeCluster(cl)
mpi.quit()

