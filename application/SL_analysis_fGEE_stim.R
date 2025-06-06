#----------------------------
# SL data analysis with FGEE
#----------------------------
library(data.table)
library(ggplot2)
library(gridExtra)
library(dplyr)
# source files
source("~/Desktop/functional_GEE/gee1step-main/R/fun.gee1step.plot.R")
source("/Users/Desktop/Research/functional_GEE/gee1step-main/R/fun.gee1step_joint.R")
wd3 <- "/Users/Desktop/Research/Causal/msm_hr/Pen_Spline/Sims/Code/"
source(paste0(wd3,"sandwich_normal_pen_fGEE.R"))

# read data
wd <- "/Users/Desktop/Research/functional_GEE/data/SL/Chronic/matlab/pooled/pre_process"
save.wd <- "/Users/Desktop/Research/functional_GEE/Figures/Lee_data"
setwd(wd)
dat <- data.table::fread("SL_s1_stim.csv")

dat$nID <- paste0(dat$id, "_", dat$neuron) # unique combo of ID and neuron #

# downsample number of neurons and timepoints
set.seed(2)
ntrg <- 500
col_targ <- 120 # number of trial timepoints
neurons <- unique(dat$nID)
n_incl <- sample(neurons, ntrg)
epoch_incl <- 300:600
dat <- dat[nID %in% n_incl,]
nm <- colnames(dat)
fn.len <- sum(grepl("whisk.w", nm, fixed = TRUE))
hz <- 1 / 0.03373076
pre.stim <- round(hz) # 1 second pre stim

# evenly sample across fn. domain
sub.idx <- 1:col_targ # do not downsample

# find column numbers of different covariates and down-sample
whisk.idx <- which(grepl("whisk.", nm, fixed = TRUE))[sub.idx]
spd.idx <- which(grepl("sp", nm, fixed = TRUE))[sub.idx]
spk.idx <- which(grepl("spikes.", nm, fixed = TRUE))[sub.idx]

# create reward number
n.tab <- dat %>% 
  as_tibble() %>% 
  select(epoch, id, trial) %>%
  group_by(id, trial) %>% 
  summarise(n = max(epoch), .groups = 'drop') %>%
  group_by(id) %>% 
  mutate(N = cumsum(n)) %>%
  ungroup() %>%
  bind_rows(data.frame(id = unique(dat$id), trial = 0, N = 0, n = 0)) %>% 
  mutate(trial = trial+1) %>%
  arrange(id, trial) #%>% print(n = "Inf")

dat <- left_join(dat, n.tab, by = c("id", "trial"))
dat$epoch.num <- dat$epoch + dat$N # culmulative trial
dat <- dat[dat$epoch.num %in% epoch_incl,]
dat$epoch.num <- dat$epoch.num - min(epoch_incl) + 1

# construct data frame for pffr()
setDF(dat)
d <- data.frame(ID = as.numeric(as.factor(dat$nID)),
                trial = dat$trial,
                neuron = dat$neuron,
                epoch = dat$epoch.num,
                stim = dat$stim,
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
fit.ar <- fun.gee1step(formula = Y ~ stim,
                    data = d,
                    cluster = "ID",
                    family="binomial",
                    algorithm = "bam",
                    rho.smooth = TRUE,
                    weighted = TRUE, # use exchangeable weights or not
                    exact = FALSE,
                    time = "epoch",
                    index = "yindex.vec",
                    cv = "fastkfold",
                    cv.grid = cv.grid,
                    cov.type = cov.type,
                    knots = 10, 
                    bs = "bs",
                    sandwich = FALSE,
                    parallel = FALSE, 
                    wd = wd)

plot.ar1 <- fgee.plot(fit.ar, 
                     xlab = "Time (sec)", 
                     title_names = c("AR1: Intercept", "AR1: Stimulation"), 
                     align_x = 1, 
                     x_rescale = hz, 
                     y_val_lim = 1,
                     return = FALSE) 

# -------------
# exchangeable 
# -------------
cov.type <- "exchangeable"
fit.ex <- fun.gee1step(formula = Y ~ speed,
                    data = d,
                    cluster = "ID",
                    family="binomial",
                    algorithm = "bam",
                    rho.smooth = FALSE,
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

# standard sandwich - independent
fit <- fit.ar
pen_mat <- pen_mat_construct(fit$pffr_original, 
                             nonSmooths=0)
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
fit$model$coefficients <- fit$beta <- fit$pffr_original$coefficients
fit$model$Vp <- fit$model$Ve <- fit$model$Vc <- var_array

# sandwich on original fit
plot.sand <- fgee.plot(fit, 
                       xlab = "Time (sec)", 
                       title_names = c("Ind: Intercept", "Ind: Stimulation"), 
                       align_x = 1, 
                       x_rescale = hz, 
                       y_val_lim = 1,
                       return = FALSE) 


# exchangeable rho = 0 everywhere so equivalent to independence 
plot.ex <- fgee.plot(fit, 
                     xlab = "Time (sec)", 
                     title_names = c("Exch: Intercept", "Exch: Stimulation"), 
                     align_x = 1, 
                     x_rescale = hz, 
                     y_val_lim = 1,
                     return = FALSE) 

# save
setwd(save.wd)

# combine and save
plot.list <- list(plot.sand, plot.ex, plot.ar1)
plt_comb <- do.call("grid.arrange", c(plot.list, nrow = 3))

setwd(save.wd)
ggsave( "stim-comb.pdf",
        plot = plt_comb,
        width = 5,
        height = 7)

saveRDS(plot.list, file = "model.list_stim.Rdata")

# to plot saved objects use the following
mod.list <- readRDS(file="model.list.Rdata")
grid::grid.draw(mod.list[[1]]) # first plot in 3-element list
