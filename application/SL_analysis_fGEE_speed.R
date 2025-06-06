#----------------------------
# SL data analysis with FGEE
#----------------------------
library(data.table)
library(ggplot2)
library(gridExtra)
# source files
source("~/Desktop/functional_GEE/gee1step-main/R/fun.gee1step.plot.R")
source("/Users/Desktop/Research/functional_GEE/gee1step-main/R/fun.gee1step_joint.R")
wd3 <- "/Users/Desktop/Research/Causal/msm_hr/Pen_Spline/Sims/Code/"
source(paste0(wd3,"sandwich_normal_pen_fGEE.R"))

# read data
wd <- "/Users/Desktop/Research/functional_GEE/data/SL/Chronic/matlab/pooled/pre_process"
save.wd <- "/Users/Desktop/Research/functional_GEE/Figures/Lee_data"
setwd(wd)
dat <- data.table::fread("SL_s1_speed.csv")

dat$nID <- paste0(dat$id, "_", dat$neuron) # unique combo of ID and neuron #

# downsample number of neurons and timepoints
set.seed(2)
ntrg <- 500#500
col_targ <- 150
neurons <- unique(dat$nID)
n_incl <- sample(neurons, ntrg)
dat <- dat[nID %in% n_incl,]
nm <- colnames(dat)
fn.len <- sum(grepl("whisk.w", nm, fixed = TRUE))
hz <- 1 / 0.03373076

# evenly sample across fn. domain
sub.idx <- 1:col_targ # do not downsample

# find column numbers of different covariates and down-sample
whisk.idx <- which(grepl("whisk.w", nm, fixed = TRUE))[sub.idx]
spd.idx <- which(grepl("speed.spd", nm, fixed = TRUE))[sub.idx]
spk.idx <- which(grepl("spikes.s", nm, fixed = TRUE))[sub.idx]

# construct data frame for pffr()
setDF(dat)
d <- data.frame(ID = as.numeric(as.factor(dat$nID)),
                trial = dat$trial,
                neuron = dat$neuron,
                speed = I(dat[,spd.idx]),
                avg_spd = rowMeans(dat[,spd.idx]),
                Y = I(dat[,spk.idx]),
                whisk = I(dat[,whisk.idx]))

# center and scale whisker stim distibution
d$whisk <- scale(d$whisk)
d$speed <- scale(d$speed)

# cv grid
cv.grid <- vector(length = 3, "list")
# pre-grid (never include 0 in pre-grid!)
cv.grid[[1]] <- sort(unique( c(10^(-6:4), 
                               10^(-6:4) * 0.25,
                               10^(-6:4) * 0.5,
                               10^(-6:4) * 0.75,
                               10^(-6:4) * 0.9)))
cv.grid[[2]] <- c(0, 1e-2, 0.1, 1, 10, 100)
cv.grid[[3]] <- c(0, 0.5, 0.75, 1, 1.3, 2, 5) 

#----------------------------------------------------
# fit model
timeStart <- Sys.time()
L <- ncol(d$Y)

wd <- "/Users/Desktop/Research/functional_GEE/gee1step-main/sources/"
cov.type <- "ar1"
fit.ar <- fun.gee1step(formula = Y ~ speed,
                    data = d,
                    cluster = "ID",
                    family="binomial",
                    algorithm = "bam",
                    rho.smooth = TRUE,
                    weighted = TRUE, # use exchangeable weights or not
                    exact = FALSE,
                    time = "trial",
                    index = "yindex.vec",
                    cv = "fastkfold",
                    cv.grid = cv.grid,
                    cov.type = cov.type,
                    knots = 20, 
                    bs = "bs",
                    sandwich = FALSE,
                    parallel = FALSE, 
                    wd = wd)

timeEnd <- Sys.time()
difftime(timeEnd, timeStart, units = "mins")

plot.ar1 <- fgee.plot(fit.ar, 
                     xlab = "Time (sec)", 
                     title_names = c("AR1: Intercept", "AR1: Speed"), 
                     align_x = NULL, 
                     x_rescale = hz, 
                     y_val_lim = 1,
                     return = FALSE) 

# -------------
# exchangeable 
# -------------
timeStart <- Sys.time()
cov.type <- "exchangeable"
fit.ex <- fun.gee1step(formula = Y ~ speed,
                    data = d,
                    cluster = "ID",
                    family="binomial",
                    algorithm = "bam",
                    rho.smooth = TRUE,
                    weighted = TRUE, 
                    exact = FALSE,
                    time = "trial",
                    index = "yindex.vec",
                    cv = "fastkfold",
                    cv.grid = cv.grid,
                    cov.type = cov.type,
                    knots = 20, 
                    bs = "bs",
                    sandwich = FALSE,
                    parallel = FALSE, 
                    wd = wd)

timeEnd <- Sys.time()
difftime(timeEnd, timeStart, units = "mins")

plot.ex <- fgee.plot(fit.ex, 
                   xlab = "Time (sec)", 
                   title_names = c("Exch: Intercept", "Exch: Speed"), 
                   align_x = NULL, 
                   x_rescale = hz, 
                   y_val_lim = 1,
                   return = FALSE) 

# standard sandwich - independent
fit <- fit.ar
pen_mat <- pen_mat_construct(fit$pffr_original, 
                             nonSmooths=0)
fit$model$dat <- fit$data
x <- as.data.frame(fit$data)
rm.vec <- c("Y", "cluster", "yindex.vec", "time")

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
fit$model$coefficients <- fit$beta <- fit$pffr_original$coefficients
fit$model$Vp <- fit$model$Ve <- fit$model$Vc <- var_array

# sandwich on original fit
plot.sand <- fgee.plot(fit, 
                       xlab = "Time (sec)", 
                       title_names = c("Ind: Intercept", "Ind: Speed"), 
                       align_x = NULL, 
                       x_rescale = hz, 
                       y_val_lim = 1,
                       return = FALSE) 


# combine and save
setwd(save.wd)
plot.list <- list(plot.sand, plot.ex, plot.ar1)
plt_comb <- do.call("grid.arrange", c(plot.list, nrow = 3))

setwd(save.wd)
ggsave( "whisk-activit-comb.pdf",
        plot = plt_comb,
        width = 5,
        height = 7)

saveRDS(plot.list, file = "model.list_whisk_activity.Rdata")

# to plot saved objects use the following
mod.list <- readRDS(file="model.list_whisk_activity.Rdata")
grid::grid.draw(mod.list[[1]]) # first plot in 3-element list
