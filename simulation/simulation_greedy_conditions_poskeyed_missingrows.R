####------------------- simulation testinfo -----------------------####
# devtools::load_all()
# library(doParallel)

library(doMPI)
library(MFCblockInfo)
library(mvtnorm)
library(lpSolveAPI)

####------------------- simulation design -------------------------####

# design.load.all <- readRDS("simulation/design_load_all_5-15_equal.rds")
design.load.all <- readRDS("design_load_all_5-15_equal.rds")

factor.blocksize <- c(2,3,4)
factor.keying <- "0" # c("0","12","23")
factor.int <- "large" #c("small","large")
factor.load <- "acceptable" #c("good","acceptable")
factor.algorithm <- c("greedy-a","greedy-d","opt-t","r2","loads","random")
factor.target <- c("weighted","equal")
factor.ntraits <- c("5") #c("5,","15)

#number of replications
R <-  200 #500

design.sim <- expand.grid("blocksize"=factor.blocksize, "keying"=factor.keying, "intercepts"=factor.int,
                          "loads"=factor.load, "target"=factor.target, "ntraits"=factor.ntraits, "rep"=1:R)

#read in which rows are missing
id.miss <- readRDS("missing_rows_greedy_conditions_poskeyed.rds")
id.design <- paste(design.sim$blocksize, design.sim$target, design.sim$constraints, design.sim$rep, sep="_")
#compare and extract which rows (d) in design.sim are missing
d.miss <- which(id.design %in% setdiff(id.design, id.miss))
design.sim <- design.sim[d.miss,]

####-------------------- fixed conditions -------------------------####

#trait correlations (Big 5 from meta-analysis van der Linden et al.)
trait.cov <- matrix(c(1,-.36,-.17,-.36,-.43,
                      -.36,1,.43,.26,.29,
                      -.17,.43,1,.21,.20,
                      -.36,.26,.21,1,.43,
                      -.43,.29,.20,.43,1),
                    nrow=5,ncol=5)

int.range <- list("small"=c(-1,1), "large"=c(-2,2))
load.range <- list("good"=c(.65,.95), "acceptable"=c(.45,.95))

#create grid of traits if traits are not given
tr.levels <- c(-1,0,1)
tr.list <- vector("list", ncol(trait.cov))
for(tr in 1:length(tr.list)) tr.list[[tr]] <- tr.levels
traits.grid <- expand.grid(tr.list)

#reduce for testing
# traits.grid <- traits.grid[1:5,]
# traits.grid <- rbind(traits.grid, rep(0,5))

J <- 500
####------------------ start simulation -------------------####

# cl <- makeCluster(4)
# registerDoParallel(cl)
# mpiopts <- NULL

cl <- startMPIcluster()
registerDoMPI(cl)

sinkWorkerOutput(paste0("worker_iter_greedy.out"))

# define chunkSize so that each cluster worker gets a single task chunk
chunkSize <- ceiling(R/getDoParWorkers())
mpiopts <- list(chunkSize=chunkSize)

res <- foreach (d=1:nrow(design.sim), .combine=rbind, .verbose=T, .packages=c("mvtnorm","numDeriv","devtools","lpSolveAPI","MFCblockInfo"),
                .inorder=F, .errorhandling="remove", .options.mpi=mpiopts) %dopar% {
                  
                  set.seed(1204+d)
                  
                  #select design load according to test design condition
                  design.load <- design.load.all[[as.character(design.sim[d,"ntraits"])]][[as.character(design.sim[d,"blocksize"])]][[as.character(design.sim[d,"keying"])]]
                  #quadruple
                  design.load <- rbind(design.load, design.load, design.load, design.load)
                  
                  
                  nb <- design.sim[d,"blocksize"]
                  K <- nrow(design.load)/nb
                  blocks <- create.block.ind(K, nb)
                  
                  #12 pairwise comparisons, constant across block sizes
                  K.start <- 12/choose(nb,2)
                  K.final <- K/4
                  
                  ####------------------ simulate item parameters -------------------####
                  
                  #select blocks.start so that all traits are included (fixed)
                  if (design.sim[d,"blocksize"]==2) {
                    blocks.start <- seq(1, 144, 12)
                  } else {
                    blocks.start <- 1:K.start
                  }
                  blocks.part2 <- (1:K)[! (1:K) %in% blocks.start]
                  #re-order design.load
                  design.load <- design.load[c(t(blocks[c(blocks.start, blocks.part2),])),]
                    
                  #separately for first and second half
                  items1 <- sim.items(design.load=design.load[c(t(blocks[1:K.start,])),], K=K.start, nb=nb,
                                      load.range=c(.65,.95),
                                      int.range=c(-1,1))
                  items2 <- sim.items(design.load=design.load[c(t(blocks[(K.start+1):K,])),], K=K-K.start, nb=nb,
                                      load.range=load.range[[as.character(design.sim[d,"loads"])]],
                                      int.range=int.range[[as.character(design.sim[d,"intercepts"])]])
                  items <- rbind(items1, items2)
                  
                  ####-------------------------------- traits and responses --------------------------####
                  
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
                  
                  ####------------------- calculate info --------------------------------####
                  load.mat <- items$loads * design.load
                  gamma.true <- create.design.mat(K=K, nb=nb) %*% items$u.mean
                  
                  infos <- calc.info.block(lhb.mplus, traits=as.matrix(traits.grid), int=gamma.true, loads=load.mat, uni=diag(items$uni),
                                           K=K, nb=nb)
                  
                  ####-------------- greedy algorithm with A-optimality and D-optimality ----------------####
                  ind.a.opt <- select.greedy(infos, calc.a.optimality, traits.grid, K, K.start, K.final, maximize=FALSE, weights.grid=weights.grid)
                  ind.d.opt <- select.greedy(infos, calc.d.optimality, traits.grid, K, K.start, K.final, maximize=TRUE, weights.grid=weights.grid)
                  
                  ####---------------- MIP T-optimality ----------------------------------####
                  #trace for each block (and grid point)
                  info.trace <- calc.info.trace(infos)
                  
                  #across blocks (for weights on grid points)
                  info.trace.pool <- rowSums(info.trace)
                  info.trace.pool.start <- rowSums(info.trace[,1:K.start])
                  
                  #without start
                  info.trace <- info.trace[,-c(1:K.start)]
                  results.opt <- select.optimal(info.trace, traits.grid, info.trace.pool.start, K-K.start, K.final-K.start, weights.grid=weights.grid)
                  ind.t.opt <- c(1:K.start, K.start + results.opt$ind.opt)
                  
                  ####---------------- block selection based on R^2 ----------------------####
                  #R^2 mean across persons and traits, weighted
                  if(as.character(design.sim[d,"target"])=="weighted") {
                    #weights for R^2
                    a.all <- calc.a.optimality(infos)
                    a.all <- a.all[which(rowSums(traits.grid==0)==ncol(traits.grid))]/a.all
                  }
                  
                  means.r2 <- do.call(c, lapply(1:K, function(k,i,a) mean(rowMeans(calc.info.block.r2(i, wo.blocks=k))*a), i=infos, a=a.all))
                  ind.r2 <- c(1:K.start, K.start + order(means.r2[(K.start+1):K], decreasing = T)[1:(K.final - K.start)])
                  
                  ####---------------- item selection based on loadings ------------------####
                  loads.blocks <- t(apply(blocks, 1, function(b, dl) colSums(dl[b,]), dl=load.mat))
                  means.loads <- rowMeans(abs(loads.blocks))
                  ind.loads <- c(1:K.start, K.start + order(means.loads[(K.start+1):K], decreasing = T)[1:(K.final - K.start)])
                  
                  ####---------------- random item selection ------------------####
                  ind.rand <- c(1:K.start, sample((K.start + 1):K, K.final-K.start))
                  
                  ####---------------- trait estimation and summary measures -------------####
                  
                  ind.list <- list("greedy-a"=ind.a.opt, "greedy-d"=ind.d.opt, "opt-t"=ind.t.opt,
                                   "r2"=ind.r2, "loads"=ind.loads, "random"=ind.rand)
                  res.r <- NULL
                  for (a in factor.algorithm) {
                    
                    if(anyNA(ind.list[[as.character(a)]])) {
                      res.r <- rbind(res.r,
                                     data.frame(design.sim[d,], "trait"=1:ncol(design.load), "algorithm"=a, rec=NA, RMSE=NA, MAB=NA))
                    } else {
                      blocks.ind <- c(t(blocks[ind.list[[as.character(a)]],]))
                      gamma.true <- create.design.mat(K=K.final, nb=nb) %*% items$u.mean[blocks.ind]
                      
                      #estimate traits based on new questionnaire
                      estimates <- est.MAP(FUN=lhb.mplus, responses=responses$rankindices[,ind.list[[as.character(a)]]],
                                           int=gamma.true, loads=load.mat[blocks.ind,], uni=diag(items$uni)[blocks.ind,blocks.ind],
                                           perms=permute(1:nb), nb=nb,
                                           m.prior=rep(0,ncol(design.load)), s.prior=trait.cov, SE=FALSE)
                      
                      rec <- diag(cor(estimates$traits, traits))
                      RMSE <- colMeans((estimates$traits - traits)^2)
                      MAB <- colMeans(abs(estimates$traits - traits))
                      
                      res.r <- rbind(res.r,
                                     data.frame(design.sim[d,], "trait"=1:ncol(design.load), "algorithm"=a, rec, RMSE, MAB))
                    }
                  }
                  saveRDS(res.r, file=paste0("results_greedy_conditions_poskeyed_missingrows/results_greedy_d",d,".rds"))
                  res.r
                }

saveRDS(res, file="results_simulation_greedy_conditions_poskeyed_missingrows.rds")

# stopCluster(cl)
closeCluster(cl)
mpi.quit()

