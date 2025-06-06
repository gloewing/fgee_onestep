library(lme4) ## mixed models
library(refund) ## fpca.face 
library(dplyr) ## organize lapply results
library(mgcv) ## smoothing in step 2
library(mvtnorm) ## joint CI
library(parallel) ## mcapply
library(Rfast)
library(caret)
library(data.table)
library(sanic)
library(SuperGauss)

simParams <- expand.grid(c(25, 50, 100),
                         c("fastkfold"), #c(TRUE, "kfold", "fastCV"), 
                         c(TRUE), 
                         c(FALSE),  # c(FALSE, TRUE, "fastBoot")
                         c(5, 25, 100), 
                         c(100), # L
                         c(0.25, 0.5, 0.75) 
)
colnames(simParams) <- c("n", "cv", "rho.smooth", "var.type", "Ni", "L", "rho")

cluserInd <- TRUE

# save and load
args = commandArgs(TRUE)
runNum <- as.integer( as.numeric(args[1]) )
iter <- as.integer( Sys.getenv('SLURM_ARRAY_TASK_ID') ) # seed index from array id

# paths
wd <- "/home/photometry_fglmm/code_utils/"
data_path <-"/home/photometry_fglmm/data" # path to original Matlab files for data
save.folder <- "/home/fgee/fgee_sim_gauss_ar" # photometry sims # 26 had fixed knots, 27 had varied
saveWD <- "/home/smtl/code/" 
source(paste0(saveWD, "SimFn.R"))
source(paste0(saveWD, "saveFn.R"))
# ********* must add working directories below especially wd4 that is called in function for robust sandwich stimator********
source(paste0(wd,"sandwich_normal_pen_fGEE.R"))
source(paste0(wd, "sim_func_data_Li.R"))
source(paste0(wd, "fun.gee1step_joint.R")) # function code

  
# simulation parameters
reps <- 300 # replicates of simulations
boots <- 5000 # for bootstrap permutation test
nknots <- 10 #
n_sim <- simParams[runNum, 1] 
N_i <- simParams[runNum, 5]
cv <- simParams[runNum, 2] 
rho.smooth <- simParams[runNum, 3] 
var.type <- simParams[runNum, 4] 
bs <- "ps"
comparison_methods <- TRUE 
L <- simParams[runNum, 6] 
p <- 3 # number of covariates
nknots_min <- 10
z11 <- 5
z12 <- 3 
sigma_sq <- 10 # 
n_delta <- "n.5" # variability in sample size
cov.type <- "ar1"
comb_pen <- TRUE
rho <- simParams[runNum, 7]
rho.vec <- rep(rho, L) # ar9 uses this
long_rho <- rho

# cv tuning grid
cv.grid <- vector(length = 3, "list")
# pre-grid (never include 0 in pre-grid!)
cv.grid[[1]] <- sort(unique( c(10^(-6:4), 
                               10^(-6:4) * 0.25,
                               10^(-6:4) * 0.5,
                               10^(-6:4) * 0.75,
                               10^(-6:4) * 0.9)))
cv.grid[[2]] <- c(0, 1e-2, 0.1, 1, 10, 100)
cv.grid[[3]] <- c(0, 0.5, 0.75, 1, 1.3, 2, 5) 

cv.grid_pre <- NULL # not necessary anymore

cnt.idx <- 12 + p # length of result vector for one method (before lambda)
c.up <- 13 #cnt.idx - 1 # for lambda indices

fileNm <- paste0("fgee_sims_",
                 "-", cov.type, 
                 "_comPen_", comb_pen,
                 "_rho_sm_", rho.smooth,
                 "_nkn_", nknots,
                 "_z11_", z11,
                 '_z12_', z12,
                 "_sSq_", sigma_sq,
                 "_nD_", n_delta,
                 "_bs_", bs,
                 "_var_", var.type,
                 "_cv_", cv,
                 "_L_", L,
                 "_rho_", rho,
                 "_Ni_", N_i,
                 "_n_trgt_", n_sim)

# results matrix
res_mat <- matrix(NA, ncol = (cnt.idx)*8, nrow = reps)
nm_vec <- c("pCI_int", "jCI_int", "jCI", "pCI", "rmse", 
            "jCI_slope", "pCI_slope", "slope_rmse", 
            "slope_bias", "pSig", "jSig", "time",
            paste0("lambda", seq(0,p-1)))

colnames(res_mat) <- c(paste0(c("fGEE_"), nm_vec),
                       paste0(c("pffr_"), nm_vec),
                       paste0(c("robust_"), nm_vec),
                       paste0(c("exact_"), nm_vec),
                       paste0(c("unwght_"), nm_vec),
                       paste0(c("ar1_"), nm_vec),
                       paste0(c("comb_"), nm_vec),
                       paste0(c("Li_"), nm_vec))
   

##################
# simulate data
##################
set.seed(iter)
beta_idx <- 2:3 # this corresponds to slope effects
n_delta <- round(N_i/2)

dat_sim <-  GenerateData(Nsubj=n_sim, 
                         numFunctPoints = L, 
                         min_visit=N_i,
                         max_visit=N_i,
                         numLongiPoints = N_i, 
                         sigma_sq = sigma_sq, 
                         sigma_z11 = z11, sigma_z12 = z12, 
                         sigma_z21 = 2, sigma_z22 = 1,
                         corstr = cov.type,
                         rho.func = long_rho, # within function correlation
                         rho.vec = rho.vec)

dat_sim$beta <- t(dat_sim$beta)
#--------------------------------------------------------------------------------------------------------------------
#########################################################################
# analyze data with various approaches and assess 95% CI coverage 
#########################################################################

#----------------------------------------------------------
# A) functional GEE
#----------------------------------------------------------
##############################
# fit model to simulated data
##############################
# fit model and get model CIs
timeStart <- Sys.time()
Y <- as.data.frame(dat_sim$data$Y)
L <- ncol(Y)
X <- model.matrix(~X1 + X2, data = data.frame(dat_sim$data$Cov))
d = data.frame(Y = I(Y), ID = dat_sim$data$subjID, X, time = dat_sim$time)

fit <- fun.gee1step(formula = Y ~ X1 + X2,
                    data = d,
                    cluster = "ID",
                    family="gaussian",
                    algorithm = "bam",
                    rho.smooth = rho.smooth,
                    weighted = TRUE, # use exchangeable weights or not
                    exact = FALSE,
                    time = "time",
                    index = "yindex.vec",
                    cv = cv,
                    cv.grid = cv.grid,
                    cov.type = cov.type,
                    knots = nknots_min, 
                    m.pffr = c(2, 2),
                    bs = bs,
                    sandwich = var.type,
                    parallel = FALSE, 
                    wd = wd)

fit1 <- fit
timeEnd <- Sys.time()

##############################
# generate 95% CIs
##############################
cnt <- 0
lower <- upper <- lower.joint <- upper.joint <- matrix(NA, nrow = p, ncol = L) # initialize CIs
joint_incl <- naive_incl <- lower # initialize inclusion indicator matrix
dd <- data.frame(X1 = 1, X2 = 1, yindex.vec = 1:L)
fit$model$Ve <- fit$model$Vc <- fit$model$Vp <- fit$Vp <- fit$vb
preds <- mgcv::predict.gam(fit$model, dd, type="iterms", se.fit=TRUE)
preds$fit[,1]=preds$fit[,1] + fit$beta[1]
fit$betaHat <- t(preds$fit)
preds$fit <- lapply(1:nrow(fit$betaHat), function(iii) matrix(preds$fit[,iii], nrow = 1) )
preds$se <- lapply(1:nrow(fit$betaHat), function(iii) matrix(preds$se[,iii], nrow = 1) )

# intercept
d0 <- data.frame(X1 = 0, X2 = 0, yindex.vec = 1:L)
preds0 <- mgcv::predict.gam(fit$model, d0, type="iterms", se.fit=TRUE)
preds$se[[1]][1,] <- preds0$se.fit[,1]

# iterate across regression coefficients, calculate PIs and determine inclusion across fn. domain
for(r in 1:p){
  # joint CIs
  lower.joint[r,] = preds$fit[[r]][1,] - fit$qn[[r]] * preds$se[[r]][1,] #fit$betaHat[r,] - fit$qn[r]*sqrt(diag(fit$betaHat.var[,,r]))
  upper.joint[r,] = preds$fit[[r]][1,] + fit$qn[[r]] * preds$se[[r]][1,] #fit$betaHat[r,] + fit$qn[r]*sqrt(diag(fit$betaHat.var[,,r]))
  
  # check whether estimate in CI
  joint_incl[r,] <- dat_sim$beta[r,] <= upper.joint[r,] & dat_sim$beta[r,] >= lower.joint[r,]
  
  # naive CIs
  lower[r,] = preds$fit[[r]][1,] - 1.96 * preds$se[[r]][1,] 
  upper[r,] = preds$fit[[r]][1,] + 1.96 * preds$se[[r]][1,] 
  
  # check whether estimate in CI
  naive_incl[r,] <- dat_sim$beta[r,] <= upper[r,] & dat_sim$beta[r,] >= lower[r,]
}

# average of CI inclusion intercept
res_mat[iter, cnt+1] <- mean( naive_incl[1, ] )  # whether pointwise CI contains estimate
res_mat[iter, cnt+2] <- mean( joint_incl[1, ] )  # whether joint CI contains estimate

# calculate regression coefficient inclusion in 95% CI (across functional domain)
joint_incl <- rowSums(joint_incl) / ncol(joint_incl)
naive_incl <- rowSums(naive_incl) / ncol(naive_incl)

# average CI coverage across regression coefficients
res_mat[iter,cnt+3] <- mean(joint_incl)
res_mat[iter,cnt+4] <- mean(naive_incl)
res_mat[iter,cnt+5] <- sqrt( mean( (dat_sim$beta - fit$betaHat)^2 ) ) # rmse of model fit

# evaluate performance w/o intercept (so only cue effect here since one covariate)
res_mat[iter,cnt+6] <- mean(joint_incl[-1])
res_mat[iter,cnt+7] <- mean(naive_incl[-1])
res_mat[iter,cnt+8] <- sqrt( mean( (dat_sim$beta[-1,] - fit$betaHat[-1,])^2 ) ) # rmse of model fit
res_mat[iter,cnt+9] <- mean(dat_sim$beta[beta_idx,] - fit$betaHat[beta_idx,]) # bias of cue term

# average of significance
res_mat[iter, cnt+10] <- mean( 1 * I(upper[beta_idx, ] * lower[beta_idx, ] > 0) )
res_mat[iter, cnt+11] <- mean( 1 * I(upper.joint[beta_idx, ] * lower.joint[beta_idx, ] > 0) )
res_mat[iter, cnt+12] <- as.numeric(difftime(timeEnd, timeStart, units = "secs"))

# penalty hyperparameter value
if(length(cv.grid) > 1){
  res_mat[iter,seq(cnt+c.up, cnt+c.up + p - 1)] <- as.numeric(fit$lambda)
}else{
  res_mat[iter,seq(cnt+c.up, cnt+c.up + p - 1)] <- rep(cv.grid, 2)
}
#--------------------------------------------------------------------------------------------------------------------

if(comparison_methods){
  #----------------------------------------------------------
  # B) pffr with robust variance estimator
  #----------------------------------------------------------
  
  ###############################################
  # use original pffr model's coefficients to calculate variance
  ###############################################
  timeStart <- Sys.time()
  pen_mat <- pen_mat_construct(fit$pffr_original, 
                               nonSmooths=0)
  
  # fit for timing
  fit_pffr <- refund::pffr(formula = Y ~ X1 + X2,
                           algorithm = "bam",
                           family = "gaussian",
                           discrete = TRUE,
                           bs.yindex = list(bs = bs, 
                                            k = nknots_min+1,
                                            m = c(2, 2)),
                           data = d)
  
  fit$model$dat <- fit$data
  x <- as.data.frame(fit$data)
  if(cov.type == "ar1"){
    rm.vec <- c("Y", "cluster", "yindex.vec", "time")
  }else{
    rm.vec <- c("Y", "cluster", "yindex.vec")
  }
  
  var_array <- sandwich_normal_pen(data =  as.data.frame(fit$data), 
                                   model = fit$model, 
                                   pen_mat = pen_mat,
                                   X = x[!colnames(x) %in% rm.vec],
                                   Y = fit$data$Y,
                                   w= rep(1, nrow(fit$data)),
                                   id_name = "cluster",
                                   beta_hat = fit$pffr_original$coefficients,
                                   nonSmooths = 0, wd = wd)$cov # our estimator
  
  fit$Vp <- fit$fit$vb <- fit$pffr_original$Vp <- var_array
  dd <- data.frame(X1 = 1, X2 = 1, yindex.vec = 1:L)
  fit$model$Ve <- fit$model$Vc <- fit$model$Vp <- fit$Vp <- fit$vb
  preds <- mgcv::predict.gam(fit$pffr_original, dd, type="iterms", se.fit=TRUE)
  preds$fit[,1]=preds$fit[,1] + fit$pffr_original$coefficients[1]
  fit$betaHat <- t(preds$fit)
  preds$fit <- lapply(1:nrow(fit$betaHat), function(iii) matrix(preds$fit[,iii], nrow = 1) )
  preds$se <- lapply(1:nrow(fit$betaHat), function(iii) matrix(preds$se[,iii], nrow = 1) )
  
  # intercept
  d0 <- data.frame(X1 = 0, X2 = 0, yindex.vec = 1:L)
  preds0 <- mgcv::predict.gam(fit$pffr_original, d0, type="iterms", se.fit=TRUE)
  preds$se[[1]][1,] <- preds0$se.fit[,1]

  qn <- joint.qn(glmfit = fit$pffr_original, NN.sim = 10000)
  timeEnd <- Sys.time()
  
  lower <- upper <- lower.joint <- upper.joint <- matrix(NA, nrow = p, ncol = L) # initialize CIs
  joint_incl <- naive_incl <- lower # initialize inclusion indicator matrix
  lower.avg <- upper.avg <- lower # CIs for averaging
  n <- length(unique(dat_sim$data$ID))
  
  # iterate across regression coefficients, calculate PIs and determine inclusion across fn. domain
  for(r in 1:p){
    # joint CIs
    lower.joint[r,] = preds$fit[[r]][1,] - qn[[r]] * preds$se[[r]][1,] #fit$betaHat[r,] - fit$qn[r]*sqrt(diag(fit$betaHat.var[,,r]))
    upper.joint[r,] = preds$fit[[r]][1,] + qn[[r]] * preds$se[[r]][1,] #fit$betaHat[r,] + fit$qn[r]*sqrt(diag(fit$betaHat.var[,,r]))
    
    # check whether estimate in CI
    joint_incl[r,] <- dat_sim$beta[r,] <= upper.joint[r,] & dat_sim$beta[r,] >= lower.joint[r,]
    
    # naive CIs
    lower[r,] = preds$fit[[r]][1,] - 1.96 * preds$se[[r]][1,] #fit$betaHat[r,] - 1.96*sqrt(diag(fit$betaHat.var[,,r]))
    upper[r,] = preds$fit[[r]][1,] + 1.96 * preds$se[[r]][1,] #fit$betaHat[r,] + 1.96*sqrt(diag(fit$betaHat.var[,,r]))
    
    # check whether estimate in CI
    naive_incl[r,] <- dat_sim$beta[r,] <= upper[r,] & dat_sim$beta[r,] >= lower[r,]
  }
  
  cnt <- cnt + cnt.idx
  # average of CI inclusion over intercept
  res_mat[iter, cnt+1] <- mean( naive_incl[1, ] )  # whether pointwise CI contains estimate
  res_mat[iter, cnt+2] <- mean( joint_incl[1, ] )  # whether joint CI contains estimate
  
  # calculate regression coefficient inclusion in 95% CI (across functional domain)
  joint_incl <- rowSums(joint_incl) / ncol(joint_incl)
  naive_incl <- rowSums(naive_incl) / ncol(naive_incl)
  
  # average CI coverage across regression coefficients
  res_mat[iter,cnt+3] <- mean(joint_incl)
  res_mat[iter,cnt+4] <- mean(naive_incl)
  res_mat[iter,cnt+5] <- sqrt( mean( (dat_sim$beta - fit$betaHat)^2 ) ) # rmse of model fit
  
  # evaluate performance w/o intercept (so only cue effect here since one covariate)
  res_mat[iter,cnt+6] <- mean(joint_incl[-1])
  res_mat[iter,cnt+7] <- mean(naive_incl[-1])
  res_mat[iter,cnt+8] <- sqrt( mean( (dat_sim$beta[-1,] - fit$betaHat[-1,])^2 ) ) # rmse of model fit
  res_mat[iter,cnt+9] <- mean(dat_sim$beta[beta_idx,] - fit$betaHat[beta_idx,]) # bias of cue term
  
  # average of significance
  res_mat[iter, cnt+10] <- mean( 1 * I(upper[beta_idx, ] * lower[beta_idx, ] > 0) )
  res_mat[iter, cnt+11] <- mean( 1 * I(upper.joint[beta_idx, ] * lower.joint[beta_idx, ] > 0) )
  res_mat[iter, cnt+12] <- as.numeric(difftime(timeEnd, timeStart, units = "secs"))
  
  # penalty hyperparameter value
  if(length(cv.grid) > 1){
    res_mat[iter,seq(cnt+c.up, cnt+c.up + p - 1)] <- as.numeric(fit$lambda)
  }else{
    res_mat[iter,seq(cnt+c.up, cnt+c.up + p - 1)] <- rep(cv.grid, 2)
  }
}
#-------------------------------------------------------------------
# Robust Sandwich Estimator: working independence
#-------------------------------------------------------------------
timeStart <- Sys.time()

pen_mat <- pen_mat_construct(fit$model, 
                             nonSmooths=0)
fit$model$dat <- fit$data
if(cov.type == "ar1"){
  rm.vec <- c("Y", "cluster", "yindex.vec", "time")
}else{
  rm.vec <- c("Y", "cluster", "yindex.vec")
}
var_array <- sandwich_normal_pen(data =  fit$data, 
                                 model = fit$model, 
                                 pen_mat = pen_mat,
                                 X = fit$data[!colnames(fit$data) %in% rm.vec],
                                 Y = fit$data$Y,
                                 w= rep(1, nrow(fit$data)),
                                 id_name = "cluster",
                                 beta_hat = fit$model$coefficients,
                                 nonSmooths = 0, wd = wd)$cov # our estimator

fit$Vp <- fit$fit$vb <- fit$pffr_original$Vp <- var_array
dd <- data.frame(X1 = 1, X2 = 1, yindex.vec = 1:L)
fit$model$Ve <- fit$model$Vc <- fit$model$Vp <- fit$Vp <- fit$vb
preds <- mgcv::predict.gam(fit$model, dd, type="iterms", se.fit=TRUE)
preds$fit[,1]=preds$fit[,1] + fit$beta[1]
fit$betaHat <- t(preds$fit)
preds$fit <- lapply(1:nrow(fit$betaHat), function(iii) matrix(preds$fit[,iii], nrow = 1) )
preds$se <- lapply(1:nrow(fit$betaHat), function(iii) matrix(preds$se[,iii], nrow = 1) )

# intercept
d0 <- data.frame(X1 = 0, X2 = 0, yindex.vec = 1:L)
preds0 <- mgcv::predict.gam(fit$model, d0, type="iterms", se.fit=TRUE)
preds$se[[1]][1,] <- preds0$se.fit[,1]

qn <- joint.qn(glmfit = fit$pffr_original, NN.sim = 10000) # joint CIs
timeEnd <- Sys.time()

##############################
# generate 95% CIs
##############################

lower <- upper <- lower.joint <- upper.joint <- matrix(NA, nrow = nrow(dat_sim$beta), ncol = ncol(dat_sim$beta)) # initialize CIs
joint_incl <- naive_incl <- lower # initialize inclusion indicator matrix

# iterate across regression coefficients, calculate PIs and determine inclusion across fn. domain
for(r in 1:p){
  lower.joint[r,] = preds$fit[[r]][1,] - qn[[r]] * preds$se[[r]][1,]
  upper.joint[r,] = preds$fit[[r]][1,] + 1.96 * preds$se[[r]][1,] 
  
  # check whether estimate in CI
  joint_incl[r,] <- dat_sim$beta[r,] <= upper.joint[r,] & dat_sim$beta[r,] >= lower.joint[r,]
  
  # naive CIs
  lower[r,] = preds$fit[[r]][1,] - 1.96 * preds$se[[r]][1,] 
  upper[r,] = preds$fit[[r]][1,] + 1.96 * preds$se[[r]][1,] 
  
  # check whether estimate in CI
  naive_incl[r,] <- dat_sim$beta[r,] <= upper[r,] & dat_sim$beta[r,] >= lower[r,]
}

cnt <- cnt + cnt.idx
# average of CI inclusion intercept
res_mat[iter, cnt+1] <- mean( naive_incl[1, ] )  # whether pointwise CI contains estimate
res_mat[iter, cnt+2] <- mean( joint_incl[1, ] )  # whether joint CI contains estimate

# calculate regression coefficient inclusion in 95% CI (across functional domain)
joint_incl <- rowSums(joint_incl) / ncol(joint_incl)
naive_incl <- rowSums(naive_incl) / ncol(naive_incl)

# average CI coverage across regression coefficients
res_mat[iter,cnt+3] <- mean(joint_incl)
res_mat[iter,cnt+4] <- mean(naive_incl)
res_mat[iter,cnt+5] <- sqrt( mean( (dat_sim$beta - fit$betaHat)^2 ) ) # rmse of model fit

# evaluate performance w/o intercept (so only cue effect here since one covariate)
res_mat[iter,cnt+6] <- mean(joint_incl[-1])
res_mat[iter,cnt+7] <- mean(naive_incl[-1])
res_mat[iter,cnt+8] <- sqrt( mean( (dat_sim$beta[-1,] - fit$betaHat[-1,])^2 ) ) # rmse of model fit
res_mat[iter,cnt+9] <- mean(dat_sim$beta[beta_idx,] - fit$betaHat[beta_idx,]) # bias of cue term

# average of significance
res_mat[iter, cnt+10] <- mean( 1 * I(upper[beta_idx, ] * lower[beta_idx, ] > 0) )
res_mat[iter, cnt+11] <- mean( 1 * I(upper.joint[beta_idx, ] * lower.joint[beta_idx, ] > 0) )
res_mat[iter, cnt+12] <- as.numeric(difftime(timeEnd, timeStart, units = "secs"))

# penalty hyperparameter value
if(length(cv.grid) > 1){
  res_mat[iter,seq(cnt+c.up, cnt+c.up + p - 1)] <- as.numeric(fit$lambda)
}else{
  res_mat[iter,seq(cnt+c.up, cnt+c.up + p - 1)] <- rep(cv.grid, 2)
}
#}

#----------------------------------------------------------
# C) Exact Weighted functional GEE
#----------------------------------------------------------
rm(fit)
##############################
# fit model to simulated data
##############################
# fit model and get model CIs
timeStart <- Sys.time()
Y <- as.data.frame(dat_sim$data$Y)
L <- ncol(Y)
X <- model.matrix(~X1 + X2, data = data.frame(dat_sim$data$Cov))
d = data.frame(Y = I(Y), ID = dat_sim$data$subjID, X, time = dat_sim$time)

fit <- fun.gee1step(formula = Y ~ X1 + X2,
                    data = d,
                    cluster = "ID",
                    family="gaussian",
                    algorithm = "bam",
                    rho.smooth = rho.smooth,
                    time = "time",
                    weighted = TRUE, # use exchangeable weights or not
                    exact = TRUE,
                    cv = cv,
                    cv.grid = cv.grid,
                    cov.type = cov.type,
                    knots = nknots, 
                    m.pffr = c(2, 2),
                    bs = bs,
                    sandwich = var.type,
                    parallel = FALSE, 
                    wd = wd)

##############################
# generate 95% CIs
##############################
lower <- upper <- lower.joint <- upper.joint <- matrix(NA, nrow = p, ncol = L) # initialize CIs
joint_incl <- naive_incl <- lower # initialize inclusion indicator matrix

dd <- data.frame(X1 = 1, X2 = 1, yindex.vec = 1:L)
fit$model$Ve <- fit$model$Vc <- fit$model$Vp <- fit$Vp <- fit$vb
preds <- mgcv::predict.gam(fit$model, dd, type="iterms", se.fit=TRUE)
preds$fit[,1]=preds$fit[,1] + fit$beta[1]
fit$betaHat <- t(preds$fit)
preds$fit <- lapply(1:nrow(fit$betaHat), function(iii) matrix(preds$fit[,iii], nrow = 1) )
preds$se <- lapply(1:nrow(fit$betaHat), function(iii) matrix(preds$se[,iii], nrow = 1) )

# intercept
d0 <- data.frame(X1 = 0, X2 = 0, yindex.vec = 1:L)
preds0 <- mgcv::predict.gam(fit$model, d0, type="iterms", se.fit=TRUE)
preds$se[[1]][1,] <- preds0$se.fit[,1]

fit$pffr_original$Ve <- fit$pffr_original$Vc <- fit$pffr_original$Vp <- fit$pffr_original$Vp <- fit$vb
fit$qn <- joint.qn(glmfit = fit$pffr_original, NN.sim = 10000)
timeEnd <- Sys.time()

# iterate across regression coefficients, calculate PIs and determine inclusion across fn. domain
for(r in 1:p){
  # joint CIs
  lower.joint[r,] = preds$fit[[r]][1,] - fit$qn[[r]] * preds$se[[r]][1,] #fit$betaHat[r,] - fit$qn[r]*sqrt(diag(fit$betaHat.var[,,r]))
  upper.joint[r,] = preds$fit[[r]][1,] + fit$qn[[r]] * preds$se[[r]][1,] #fit$betaHat[r,] + fit$qn[r]*sqrt(diag(fit$betaHat.var[,,r]))
  
  # check whether estimate in CI
  joint_incl[r,] <- dat_sim$beta[r,] <= upper.joint[r,] & dat_sim$beta[r,] >= lower.joint[r,]
  
  # naive CIs
  lower[r,] = preds$fit[[r]][1,] - 1.96 * preds$se[[r]][1,] 
  upper[r,] = preds$fit[[r]][1,] + 1.96 * preds$se[[r]][1,] 
  
  # check whether estimate in CI
  naive_incl[r,] <- dat_sim$beta[r,] <= upper[r,] & dat_sim$beta[r,] >= lower[r,]
}

cnt <- cnt + cnt.idx
# average of CI inclusion over reward period
res_mat[iter, cnt+1] <- mean( naive_incl[1, ] )  # whether pointwise CI contains estimate
res_mat[iter, cnt+2] <- mean( joint_incl[1, ] )  # whether joint CI contains estimate

# calculate regression coefficient inclusion in 95% CI (across functional domain)
joint_incl <- rowSums(joint_incl) / ncol(joint_incl)
naive_incl <- rowSums(naive_incl) / ncol(naive_incl)

# average CI coverage across regression coefficients
res_mat[iter,cnt+3] <- mean(joint_incl)
res_mat[iter,cnt+4] <- mean(naive_incl)
res_mat[iter,cnt+5] <- sqrt( mean( (dat_sim$beta - fit$betaHat)^2 ) ) # rmse of model fit

# evaluate performance w/o intercept (so only cue effect here since one covariate)
res_mat[iter,cnt+6] <- mean(joint_incl[-1])
res_mat[iter,cnt+7] <- mean(naive_incl[-1])
res_mat[iter,cnt+8] <- sqrt( mean( (dat_sim$beta[-1,] - fit$betaHat[-1,])^2 ) ) # rmse of model fit
res_mat[iter,cnt+9] <- mean(dat_sim$beta[beta_idx,] - fit$betaHat[beta_idx,]) # bias of cue term

# average of significance
res_mat[iter, cnt+10] <- mean( 1 * I(upper[beta_idx, ] * lower[beta_idx, ] > 0) )
res_mat[iter, cnt+11] <- mean( 1 * I(upper.joint[beta_idx, ] * lower.joint[beta_idx, ] > 0) )
res_mat[iter, cnt+12] <- as.numeric(difftime(timeEnd, timeStart, units = "secs"))

# penalty hyperparameter value
if(length(cv.grid) > 1){
  res_mat[iter,seq(cnt+c.up, cnt+c.up + p - 1)] <- as.numeric(fit$lambda)
}else{
  res_mat[iter,seq(cnt+c.up, cnt+c.up + p - 1)] <- rep(cv.grid, 2)
}
# sample size of sims (i.e., number of trials analyzed)

#--------------------------------------------------------------------------------------------------------------------

#------------------------------------------------------------
# D) Exact *Unweighted* Functional GEE - Working Independence
#------------------------------------------------------------
rm(fit)
##############################
# fit model to simulated data
##############################
# fit model and get model CIs
timeStart <- Sys.time()
Y <- as.data.frame(dat_sim$data$Y)
L <- ncol(Y)
X <- model.matrix(~X1 + X2, data = data.frame(dat_sim$data$Cov))
d = data.frame(Y = I(Y), ID = dat_sim$data$subjID, X, time = dat_sim$time)

fit <- fun.gee1step(formula = Y ~ X1 + X2,
                    data = d,
                    cluster = "ID",
                    family="gaussian",
                    algorithm = "bam",
                    rho.smooth = rho.smooth,
                    weighted = FALSE, # working independence covariance
                    exact = TRUE,
                    cv = cv,
                    cv.grid = cv.grid,
                    cov.type = "exchangeable", # keep exchangeable regardless b/c its unweighted
                    knots = nknots, 
                    m.pffr = c(2, 2),
                    bs = bs,
                    sandwich = var.type,
                    parallel = FALSE, 
                    wd = wd)

##############################
# generate 95% CIs
##############################
lower <- upper <- lower.joint <- upper.joint <- matrix(NA, nrow = p, ncol = L) # initialize CIs
joint_incl <- naive_incl <- lower # initialize inclusion indicator matrix
dd <- data.frame(X1 = 1, X2 = 1, yindex.vec = 1:L)
fit$model$Ve <- fit$model$Vc <- fit$model$Vp <- fit$Vp <- fit$vb
preds <- mgcv::predict.gam(fit$model, dd, type="iterms", se.fit=TRUE)
preds$fit[,1]=preds$fit[,1] + fit$beta[1]
fit$betaHat <- t(preds$fit)
preds$fit <- lapply(1:nrow(fit$betaHat), function(iii) matrix(preds$fit[,iii], nrow = 1) )
preds$se <- lapply(1:nrow(fit$betaHat), function(iii) matrix(preds$se[,iii], nrow = 1) )

# intercept
d0 <- data.frame(X1 = 0, X2 = 0, yindex.vec = 1:L)
preds0 <- mgcv::predict.gam(fit$model, d0, type="iterms", se.fit=TRUE)
preds$se[[1]][1,] <- preds0$se.fit[,1]

fit$pffr_original$Ve <- fit$pffr_original$Vc <- fit$pffr_original$Vp <- fit$pffr_original$Vp <- fit$vb
fit$qn <- joint.qn(glmfit = fit$pffr_original, NN.sim = 10000)
timeEnd <- Sys.time()

# iterate across regression coefficients, calculate PIs and determine inclusion across fn. domain
for(r in 1:p){
  # joint CIs
  lower.joint[r,] = preds$fit[[r]][1,] - fit$qn[[r]] * preds$se[[r]][1,] 
  upper.joint[r,] = preds$fit[[r]][1,] + fit$qn[[r]] * preds$se[[r]][1,] 
  
  # check whether estimate in CI
  joint_incl[r,] <- dat_sim$beta[r,] <= upper.joint[r,] & dat_sim$beta[r,] >= lower.joint[r,]
  
  # naive CIs
  lower[r,] = preds$fit[[r]][1,] - 1.96 * preds$se[[r]][1,] 
  upper[r,] = preds$fit[[r]][1,] + 1.96 * preds$se[[r]][1,] 
  
  # check whether estimate in CI
  naive_incl[r,] <- dat_sim$beta[r,] <= upper[r,] & dat_sim$beta[r,] >= lower[r,]
}

cnt <- cnt + cnt.idx
# average of CI inclusion intercept
res_mat[iter, cnt+1] <- mean( naive_incl[1, ] )  # whether pointwise CI contains estimate
res_mat[iter, cnt+2] <- mean( joint_incl[1, ] )  # whether joint CI contains estimate

# calculate regression coefficient inclusion in 95% CI (across functional domain)
joint_incl <- rowSums(joint_incl) / ncol(joint_incl)
naive_incl <- rowSums(naive_incl) / ncol(naive_incl)

# average CI coverage across regression coefficients
res_mat[iter,cnt+3] <- mean(joint_incl)
res_mat[iter,cnt+4] <- mean(naive_incl)
res_mat[iter,cnt+5] <- sqrt( mean( (dat_sim$beta - fit$betaHat)^2 ) ) # rmse of model fit

# evaluate performance w/o intercept (so only cue effect here since one covariate)
res_mat[iter,cnt+6] <- mean(joint_incl[-1])
res_mat[iter,cnt+7] <- mean(naive_incl[-1])
res_mat[iter,cnt+8] <- sqrt( mean( (dat_sim$beta[-1,] - fit$betaHat[-1,])^2 ) ) # rmse of model fit
res_mat[iter,cnt+9] <- mean(dat_sim$beta[beta_idx,] - fit$betaHat[beta_idx,]) # bias of cue term

# average of significance
res_mat[iter, cnt+10] <- mean( 1 * I(upper[beta_idx, ] * lower[beta_idx, ] > 0) )
res_mat[iter, cnt+11] <- mean( 1 * I(upper.joint[beta_idx, ] * lower.joint[beta_idx, ] > 0) )
res_mat[iter, cnt+12] <- as.numeric(difftime(timeEnd, timeStart, units = "secs"))

# penalty hyperparameter value
if(length(cv.grid) > 1){
  res_mat[iter,seq(cnt+c.up, cnt+c.up + p - 1)] <- as.numeric(fit$lambda)
}else{
  res_mat[iter,seq(cnt+c.up, cnt+c.up + p - 1)] <- rep(cv.grid, 2)
}
#--------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------
# E) Functional Direction Only (AR1)
#----------------------------------------------------------
##############################
# fit model to simulated data
##############################
# fit model and get model CIs
timeStart <- Sys.time()
Y <- as.data.frame(dat_sim$data$Y)
L <- ncol(Y)
X <- model.matrix(~X1 + X2, data = data.frame(dat_sim$data$Cov))
d = data.frame(Y = I(Y), ID = dat_sim$data$subjID, X, time = dat_sim$time)

fit <- fun.gee1step(formula = Y ~ X1 + X2,
                    data = d,
                    cluster = "ID",
                    family="gaussian",
                    algorithm = "bam",
                    rho.smooth = rho.smooth,
                    weighted = TRUE, # use exchangeable weights or not
                    exact = FALSE,
                    time = "time",
                    index = "time",
                    cv = cv,
                    cv.grid = cv.grid,
                    cov.type = cov.type,
                    # cv.grid_pre = cv.grid_pre,
                    knots = nknots_min, #nknots_min,
                    m.pffr = c(2, 2),
                    bs = bs,
                    sandwich = var.type,
                    parallel = FALSE, #!cluserInd,
                    wd = wd)

fit2 <- fit # save for combined below

timeEnd <- Sys.time()

##############################
# generate 95% CIs
##############################
cnt <- cnt + cnt.idx
lower <- upper <- lower.joint <- upper.joint <- matrix(NA, nrow = p, ncol = L) # initialize CIs
joint_incl <- naive_incl <- lower # initialize inclusion indicator matrix

dd <- data.frame(X1 = 1, X2 = 1, yindex.vec = 1:L)
fit$model$Ve <- fit$model$Vc <- fit$model$Vp <- fit$Vp <- fit$vb
preds <- mgcv::predict.gam(fit$model, dd, type="iterms", se.fit=TRUE)
preds$fit[,1]=preds$fit[,1] + fit$beta[1]
fit$betaHat <- t(preds$fit)
preds$fit <- lapply(1:nrow(fit$betaHat), function(iii) matrix(preds$fit[,iii], nrow = 1) )
preds$se <- lapply(1:nrow(fit$betaHat), function(iii) matrix(preds$se[,iii], nrow = 1) )

# intercept
d0 <- data.frame(X1 = 0, X2 = 0, yindex.vec = 1:L)
preds0 <- mgcv::predict.gam(fit$model, d0, type="iterms", se.fit=TRUE)
preds$se[[1]][1,] <- preds0$se.fit[,1]

# iterate across regression coefficients, calculate PIs and determine inclusion across fn. domain
for(r in 1:p){
  # joint CIs
  lower.joint[r,] = preds$fit[[r]][1,] - fit$qn[[r]] * preds$se[[r]][1,] #fit$betaHat[r,] - fit$qn[r]*sqrt(diag(fit$betaHat.var[,,r]))
  upper.joint[r,] = preds$fit[[r]][1,] + fit$qn[[r]] * preds$se[[r]][1,] #fit$betaHat[r,] + fit$qn[r]*sqrt(diag(fit$betaHat.var[,,r]))
  
  # check whether estimate in CI
  joint_incl[r,] <- dat_sim$beta[r,] <= upper.joint[r,] & dat_sim$beta[r,] >= lower.joint[r,]
  
  # naive CIs
  lower[r,] = preds$fit[[r]][1,] - 1.96 * preds$se[[r]][1,] #fit$betaHat[r,] - 1.96*sqrt(diag(fit$betaHat.var[,,r]))
  upper[r,] = preds$fit[[r]][1,] + 1.96 * preds$se[[r]][1,] #fit$betaHat[r,] + 1.96*sqrt(diag(fit$betaHat.var[,,r]))
  
  # check whether estimate in CI
  naive_incl[r,] <- dat_sim$beta[r,] <= upper[r,] & dat_sim$beta[r,] >= lower[r,]
}

# average of CI inclusion intercept
res_mat[iter, cnt+1] <- mean( naive_incl[1, ] )  # whether pointwise CI contains estimate
res_mat[iter, cnt+2] <- mean( joint_incl[1, ] )  # whether joint CI contains estimate

# calculate regression coefficient inclusion in 95% CI (across functional domain)
joint_incl <- rowSums(joint_incl) / ncol(joint_incl)
naive_incl <- rowSums(naive_incl) / ncol(naive_incl)

# average CI coverage across regression coefficients
res_mat[iter,cnt+3] <- mean(joint_incl)
res_mat[iter,cnt+4] <- mean(naive_incl)
res_mat[iter,cnt+5] <- sqrt( mean( (dat_sim$beta - fit$betaHat)^2 ) ) # rmse of model fit

# evaluate performance w/o intercept (so only cue effect here since one covariate)
res_mat[iter,cnt+6] <- mean(joint_incl[-1])
res_mat[iter,cnt+7] <- mean(naive_incl[-1])
res_mat[iter,cnt+8] <- sqrt( mean( (dat_sim$beta[-1,] - fit$betaHat[-1,])^2 ) ) # rmse of model fit
res_mat[iter,cnt+9] <- mean(dat_sim$beta[beta_idx,] - fit$betaHat[beta_idx,]) # bias of cue term

# average of significance
res_mat[iter, cnt+10] <- mean( 1 * I(upper[beta_idx, ] * lower[beta_idx, ] > 0) )
res_mat[iter, cnt+11] <- mean( 1 * I(upper.joint[beta_idx, ] * lower.joint[beta_idx, ] > 0) )
res_mat[iter, cnt+12] <- as.numeric(difftime(timeEnd, timeStart, units = "secs"))

# penalty hyperparameter value
if(length(cv.grid) > 1){
  res_mat[iter,seq(cnt+c.up, cnt+c.up + p - 1)] <- as.numeric(fit$lambda)
}else{
  res_mat[iter,seq(cnt+c.up, cnt+c.up + p - 1)] <- rep(cv.grid, 2)
}
#--------------------------------------------------------------------------------------------------------------------

####################################
# save to CSV on cluster
####################################
saveFn_new(file = res_mat,
           fileNm = fileNm,
           iterNum = iter,
           iters = reps,
           save.folder = save.folder)

####################################
# save to CSV on cluster
####################################
saveFn_new(file = res_mat,
           fileNm = fileNm,
           iterNum = iter,
           iters = reps,
           save.folder = save.folder)

#----------------------------------------------------------
####################################
# fGEE - exchangeable correlation
####################################
if(cov.type == "exchangeable"){
  timeStart <- Sys.time()
  # estimate mean functions with exchangeable covariance
  exch.estim <- fgee(formula = Y ~ X1 + X2, 
                     Y=dat_sim$data$Y, 
                     Cov=dat_sim$data$Cov,
                     s=dat_sim$data$funcArg, 
                     subjID=d$ID, 
                     Tij=dat_sim$data$Tij, 
                     numLongiPoints=L,
                     control=list(FPCA1.pve=0.95, FPCA1.knots=20),
                     corstr = "exchangeable")
  
  timeEnd <- Sys.time()
  
  fit <- exch.estim$final.fit
  fit$model <- exch.estim$final.fit
  ##############################
  # generate 95% CIs
  ##############################
  lower <- upper <- lower.joint <- upper.joint <- matrix(NA, nrow = p, ncol = L) # initialize CIs
  joint_incl <- naive_incl <- lower # initialize inclusion indicator matrix
  fit$betaHat <- t(exch.estim$beta)
  fit$qn <- list() # there is no joint CIs with this method
  fit$qn[1:p] <- 1.96
  
  # iterate across regression coefficients, calculate PIs and determine inclusion across fn. domain
  for(r in 1:p){
    # joint CIs
    lower.joint[r,] = fit$betaHat[r,] - fit$qn[[r]] * exch.estim$beta.se[,r] #fit$betaHat[r,] - fit$qn[r]*sqrt(diag(fit$betaHat.var[,,r]))
    upper.joint[r,] = fit$betaHat[r,] + fit$qn[[r]] * exch.estim$beta.se[,r] #fit$betaHat[r,] + fit$qn[r]*sqrt(diag(fit$betaHat.var[,,r]))
    
    # check whether estimate in CI
    joint_incl[r,] <- dat_sim$beta[r,] <= upper.joint[r,] & dat_sim$beta[r,] >= lower.joint[r,]
    
    # naive CIs
    lower[r,] = fit$betaHat[r,] - 1.96 * exch.estim$beta.se[,r] #fit$betaHat[r,] - 1.96*sqrt(diag(fit$betaHat.var[,,r]))
    upper[r,] = fit$betaHat[r,] + 1.96 * exch.estim$beta.se[,r] #fit$betaHat[r,] + 1.96*sqrt(diag(fit$betaHat.var[,,r]))
    
    # check whether estimate in CI
    naive_incl[r,] <- dat_sim$beta[r,] <= upper[r,] & dat_sim$beta[r,] >= lower[r,]
  }
  
  cnt <- cnt + 2*cnt.idx
  # average of CI inclusion over reward period
  res_mat[iter, cnt+1] <- mean( naive_incl[1, ] )  # whether pointwise CI contains estimate
  res_mat[iter, cnt+2] <- mean( joint_incl[1, ] )  # whether joint CI contains estimate
  
  # calculate regression coefficient inclusion in 95% CI (across functional domain)
  joint_incl <- rowSums(joint_incl) / ncol(joint_incl)
  naive_incl <- rowSums(naive_incl) / ncol(naive_incl)
  
  # average CI coverage across regression coefficients
  res_mat[iter,cnt+3] <- mean(joint_incl)
  res_mat[iter,cnt+4] <- mean(naive_incl)
  res_mat[iter,cnt+5] <- sqrt( mean( (dat_sim$beta - fit$betaHat)^2 ) ) # rmse of model fit
  
  # evaluate performance w/o intercept (so only cue effect here since one covariate)
  res_mat[iter,cnt+6] <- mean(joint_incl[-1])
  res_mat[iter,cnt+7] <- mean(naive_incl[-1])
  res_mat[iter,cnt+8] <- sqrt( mean( (dat_sim$beta[-1,] - fit$betaHat[-1,])^2 ) ) # rmse of model fit
  res_mat[iter,cnt+9] <- mean(dat_sim$beta[beta_idx,] - fit$betaHat[beta_idx,]) # bias of cue term
  
  # average of significance
  res_mat[iter, cnt+10] <- mean( 1 * I(upper[beta_idx, ] * lower[beta_idx, ] > 0) )
  res_mat[iter, cnt+11] <- mean( 1 * I(upper.joint[beta_idx, ] * lower.joint[beta_idx, ] > 0) )
  res_mat[iter, cnt+12] <- as.numeric(difftime(timeEnd, timeStart, units = "secs"))
  
  ####################################
  # save to CSV on cluster
  ####################################
  saveFn_new(file = res_mat,
             fileNm = fileNm,
             iterNum = iter,
             iters = reps,
             save.folder = save.folder)
}
