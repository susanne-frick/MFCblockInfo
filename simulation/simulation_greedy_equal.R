####------------------- simulation testinfo -----------------------####
# devtools::load_all()
# library(doParallel)

library(doMPI)
library(MFCblockInfo)
library(mvtnorm)
library(lpSolveAPI)

####------------------- simulation design -------------------------####

# design.load.all <- readRDS("simulation/design_load_all_234.rds")
design.load.all <- readRDS("design_load_all_234.rds")
design.load.1 <- design.load.all[["3"]][["12"]]
design.load <- rbind(design.load.1, design.load.1, design.load.1, design.load.1)


factor.blocksize <- 3 # c(2,3,4)
factor.keying <- "12" # c("0","12","23")
factor.int <- "large" #c("small","large")
factor.load <- "acceptable" #c("good","acceptable")
factor.algorithm <- c("greedy-a","greedy-d","opt-t","r2","loads","random")

#number of replications
R <- 500

design.sim <- expand.grid("blocksize"=factor.blocksize, "keying"=factor.keying, "intercepts"=factor.int,
                          "loads"=factor.load, "rep"=1:R)
#first reduced design: only blocksize 3, 1/2 mixed comparisons

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

                    nb <- design.sim[d,"blocksize"]
                    K <- nrow(design.load)/nb
                    blocks <- create.block.ind(K, nb)

                    K.start <- 3
                    K.final <- 20

                    ####------------------ simulate item parameters -------------------####

                    #separately for first and second half
                    items1 <- sim.items(design.load=design.load[c(t(blocks[1:(K.start),])),], K=K.start, nb=nb,
                                        load.range=c(.65,.95),
                                        int.range=c(-1,1))
                    items2 <- sim.items(design.load=design.load[c(t(blocks[(K.start+1):K,])),], K=K-K.start, nb=nb,
                                        load.range=load.range[[as.character(design.sim[d,"loads"])]],
                                        int.range=int.range[[as.character(design.sim[d,"intercepts"])]])
                    items <- rbind(items1, items2)

                    ####-------------------------------- traits and responses --------------------------####

                    #grid as traits -> to evaluate equal weights
                    traits <- rbind(traits.grid, traits.grid)

                    responses <- sim.responses(traits, items, design.load, K, nb, return.index=F)

                    ####------------------- calculate info --------------------------------####
                    load.mat <- items$loads * design.load
                    gamma.true <- create.design.mat(K=K, nb=nb) %*% items$u.mean

                    infos <- calc.info.block(lhb.mplus, traits=as.matrix(traits.grid), int=gamma.true, loads=load.mat, uni=diag(items$uni),
                                             K=K, nb=nb)

                    ####-------------- greedy algorithm with A-optimality and D-optimality ----------------####
                    ind.a.opt <- select.greedy(infos, calc.a.optimality, traits.grid,
                                               K, K.start, K.final, maximize=FALSE, weights.grid=rep(1, nrow(traits.grid)))
                    ind.d.opt <- select.greedy(infos, calc.d.optimality, traits.grid,
                                               K, K.start, K.final, maximize=TRUE, weights.grid=rep(1, nrow(traits.grid)))

                    ####---------------- MIP T-optimality ----------------------------------####
                    #trace for each block (and grid point)
                    info.trace <- do.call(rbind, lapply(infos, function(ip) apply(ip, 1, function(i) sum(diag(i)))))

                    #across blocks (for weights on grid points)
                    info.trace.pool <- rowSums(info.trace)
                    info.trace.pool.start <- rowSums(info.trace[,1:K.start])

                    #without start
                    info.trace <- info.trace[,-c(1:K.start)]
                    results.opt <- select.optimal(info.sum=info.trace, traits.grid=traits.grid,
                                                  info.start=info.trace.pool.start, K=K-K.start, K.final=K.final-K.start,
                                                  weights.grid=rep(1,nrow(traits.grid)))
                    ind.t.opt <- c(1:K.start, K.start + results.opt$ind.opt)

                    ####---------------- block selection based on R^2 ----------------------####
                    #R^2 mean across persons and traits, equal weights
                    means.r2 <- do.call(c, lapply(1:K, function(k,i) mean(rowMeans(calc.info.block.r2(i, wo.blocks=k))), i=infos))
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
                    saveRDS(res.r, file=paste0("results_greedy_equal/results_greedy_equal_d",d,".rds"))
                    res.r
                }

saveRDS(res, file="results_simulation_greedy_equal.rds")

# stopCluster(cl)
closeCluster(cl)
mpi.quit()

