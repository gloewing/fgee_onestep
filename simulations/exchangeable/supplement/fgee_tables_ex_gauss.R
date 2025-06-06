#######################################################
# functional gee simulations
#######################################################
wd <- "/Users/Desktop/Research/fgee_sim_gauss_ex"
setwd(wd)
library(dplyr)
library(latex2exp)
library(kableExtra)
library(tidyverse)
library(ggpubr)
process.wd <- "/Users/Desktop/Research/Code/" 
source(paste0(process.wd,"saveFn.R"))

process.ind <- FALSE # process files flag
#-----------------------
# simulation parameters
reps <- 300 # replicates of simulations
boots <- 5000 # for bootstrap permutation test
nknots <- 10 #
cv.vec <- c("fastkfold")
var.vec <- c(FALSE, TRUE, "fastBoot") # 
nVec <- c(25, 30, 50, 100, 500)
NiVec <- c(5, 10, 25, 50,100, 500)
rhoVec <- c( TRUE) # FALSE,
rhos <- c(0.25,0.5, 0.75)
comb_pen <- TRUE
bs <- "ps"
p <- 3 # number of covariates
z11 <- 5
z12 <- 3 # 1.5
sigma_vec <- c(0.5, 3,10)
n_delta <- "n.5" #10 
Lvec <- c(100,200)

# column names
nm_vec <- c("pCI_int", "jCI_int", "jCI", "pCI", "rmse", 
            "jCI_slope", "pCI_slope", "slope_rmse", "slope_bias", "pSig", "jSig", "time",
            paste0("lambda", seq(0,p-1)))
method_vec <- c("fGEE_", "pffr_", "robust_", "exact_", "unwght_", "ar1_", "comb_", "Li_")
ref.nm <- "pffr_"

colnm <- c(paste0(c("fGEE_"), nm_vec),
           paste0(c("pffr_"), nm_vec),
           paste0(c("robust_"), nm_vec),
           paste0(c("exact_"), nm_vec),
           paste0(c("unwght_"), nm_vec),
           paste0(c("ar1_"), nm_vec),
           paste0(c("comb_"), nm_vec),
           paste0(c("Li_"), nm_vec))

ls_p <- ls_rmse_rel <- ls_rmse_slp <- ls_pci_slp <- list()
ls_rmse <- ls_pci <- ls_jci <- ls_time <- list()
#-----------------------
itrs <- reps # number of simulation iterations
cnt <- 0
rmseMat <- matrix(nc = 4, nr = itrs)
for(r in rhos){
  for(l in Lvec){
    for(sigma_sq in sigma_vec){
      for(cv in cv.vec){
        for(var.type in var.vec){
          for(n_sim in nVec){
            for(N_i in NiVec){
              for(rho.smooth in rhoVec){
                flNm <-   paste0("fgee_sims_-exchangeable_",
                                 "comPen_", comb_pen,
                                 "_rho_sm_", rho.smooth,
                                 "_nkn_", nknots,
                                 "_z11_", z11,
                                 '_z12_', z12,
                                 "_sSq_", sigma_sq,
                                 "_nD_", n_delta,
                                 "_bs_", bs,
                                 "_var_", var.type,
                                 "_cv_", cv,
                                 "_L_", l,
                                 "_rho_", r,
                                 "_Ni_", N_i,
                                 "_n_trgt_", n_sim)
                if(process.ind){
                  filePrefix <- paste0(flNm, "_")
                  fls <- list.files(getwd()) 
                  fls_nms <- grep(paste0("^", filePrefix), fls) # files with names that match
                  
                  if(length(fls_nms) > 0 ){
                    print(flNm)
                    process_files(fileNm=flNm, save.folder = wd, itrs = itrs, colnm = colnm)
                  }   
                }
                
                flNm <- paste0(flNm, ".csv") # processed file is .csv
                
                if(file.exists(flNm)){
                  # check to see if file exists
                  cnt <- cnt + 1
                  d <- read.csv(flNm) 
                  
                  # remove NAs
                  ind <- apply(d, 1, function(x) all(is.na(x)))
                  d <- d[!ind,]
                  
                  # rmse
                  nm_vec <- paste0(method_vec, "rmse")
                  ls_rmse[[cnt]] <- cbind(gather(d[,nm_vec]), 
                                          cv, var.type, n_sim, N_i, rho.smooth, sigma_sq, l, r )
                  
                  # rmse relative
                  nm_vec <- paste0(method_vec[method_vec != ref.nm], "rmse")
                  ref.val <- paste0(ref.nm, "rmse")
                  ls_rmse_rel[[cnt]] <- cbind(gather(d[,nm_vec] / d[,ref.val]), 
                                              cv, var.type, n_sim, N_i, rho.smooth, sigma_sq, l, r )
                  
                  # pointwise coverage
                  nm_vec <- paste0(method_vec, "pCI")
                  ls_pci[[cnt]] <- cbind(gather(d[,nm_vec]), 
                                         cv, var.type, n_sim, N_i, rho.smooth, sigma_sq, l, r )
                  
                  # joint coverage
                  nm_vec <- paste0(method_vec, "jCI")
                  ls_jci[[cnt]] <- cbind(gather(d[,nm_vec]), 
                                         cv, var.type, n_sim, N_i, rho.smooth, sigma_sq, l, r )
                  
                  # power 
                  nm_vec <- paste0(method_vec, "pSig")
                  ls_p[[cnt]] <- cbind(gather(d[,nm_vec]), 
                                       cv, var.type, n_sim, N_i, rho.smooth, sigma_sq, l, r )
                  
                  # rmse
                  nm_vec <- paste0(method_vec, "slope_rmse")
                  ls_rmse_slp[[cnt]] <- cbind(gather(d[,nm_vec]), 
                                          cv, var.type, n_sim, N_i, rho.smooth, sigma_sq, l, r )
                  
                  # rmse
                  nm_vec <- paste0(method_vec, "pCI_slope")
                  ls_pci_slp[[cnt]] <- cbind(gather(d[,nm_vec]), 
                                              cv, var.type, n_sim, N_i, rho.smooth, sigma_sq, l, r )
                  
                  # time
                  nm_vec <- paste0(method_vec, "time")
                  ls_time[[cnt]] <- cbind(gather(d[,nm_vec]), 
                                             cv, var.type, n_sim, N_i, rho.smooth, sigma_sq, l, r )
                  
                  d1 <- d
                  rm(d)
                }
              }
            }
          }
        }  
      }
    }
  }
}


###########################
# RMSE Relative RMSE Table
###########################
# factors
dat <- do.call(rbind, ls_rmse_rel)
dat$cv <- as.factor(dat$cv)
dat$var.type <- as.factor(dat$var.type)
dat$n_sim <- as.factor(dat$n_sim)
dat$N_i <- as.factor(dat$N_i)
dat$rho.smooth <- as.factor(dat$rho.smooth)
dat$sigma_sq <- as.factor(dat$sigma_sq)
dat$r <- as.factor(dat$r)

library(tables)
n_included <- nVec
Ni <- 5
metric.nm <- "rmse"
tune.type <- "fastkfold"
sig <- 10
x=dat %>% tibble %>%  
  dplyr::mutate(key = sub(paste0("_", metric.nm), "", key)) %>% # remove rmse
  dplyr::filter(cv == tune.type,
                sigma_sq == sig,
                n_sim %in% n_included) %>%
  dplyr::mutate(key=str_replace_all(key, fixed("fGEE"), "One-Step"),
                key=str_replace_all(key, fixed("exact"), "GLS-Ex"),
                key=str_replace_all(key, fixed("unwght"), "GLS-Ind"),
                key=str_replace_all(key, fixed("Li"), "Marginal")
  ) %>%
  dplyr::mutate(key=fct_relevel(key, c("One-Step","GLS-Ex","GLS-Ind", "Marginal", "pffr"))) %>%
  as.data.frame()

# functions
stderr <- function(y) round(sd(y,na.rm=TRUE)/sqrt(length(y[!is.na(y)])), 2)
mm <- function(y) round(mean(y,na.rm=TRUE), 2)
x$key <- as.factor(x$key)
x$n_sim <- as.factor(x$n_sim)
x$var.type <- as.factor(x$var.type)
x$cv <- as.factor(x$cv)

stderr <- function(y) round(sd(y,na.rm=TRUE)/sqrt(length(y[!is.na(y)])), 2)
mm <- function(y) round(mean(y,na.rm=TRUE), 2)
x$key <- as.factor(x$key)
x$n_sim <- as.factor(x$n_sim)
x$var.type <- as.factor(x$var.type)
x$cv <- as.factor(x$cv)

x$cv <- droplevels(x$cv, c("kfold", "FALSE"))
x$key <- droplevels(x$key, c("robust", "comb", NA, "ar1", "NA")) # "unwght", "exact", 
x$var.type <- droplevels(x$var.type, c("TRUE", "fastBoot")) # "TRUE", 

# table
toLatex( tabular( n_sim * RowFactor(N_i,levelnames=unique(N_i))
                  ~ key*value*
                    PlusMinus(mm, stderr, digits=6, scipen = 999),
                  data = x))

x %>% dplyr::as_tibble() %>%
  group_by(N_i, n_sim, sigma_sq, key, r,cv) %>%
  summarize(mean(value, na.rm=TRUE)) %>% print(n = "Inf")


########################
# Power CI Table
########################
# factors 
dat <- do.call(rbind, ls_p)
dat$cv <- as.factor(dat$cv)
dat$var.type <- as.factor(dat$var.type)
dat$n_sim <- as.factor(dat$n_sim)
dat$N_i <- as.factor(dat$N_i)
dat$rho.smooth <- as.factor(dat$rho.smooth)
dat$r <- as.factor(dat$r)
dat$l <- as.factor(dat$l)

library(tables)
n_included <- nVec
Ni <- 5
metric.nm <- "pSig"
tune.type <- "fastkfold"
x=dat %>% tibble %>%  
  dplyr::mutate(key = sub(paste0("_", metric.nm), "", key)) %>% # remove rmse
  dplyr::filter(cv == tune.type,
                n_sim %in% n_included) %>%
  dplyr::mutate(key=str_replace_all(key, fixed("fGEE"), "One-Step"),
                key=str_replace_all(key, fixed("exact"), "GLS-Ex"),
                key=str_replace_all(key, fixed("unwght"), "GLS-Ind"),
                key=str_replace_all(key, fixed("Li"), "Marginal")
  ) %>%
  dplyr::mutate(key=fct_relevel(key, c("One-Step","GLS-Ex","GLS-Ind", "Marginal", "pffr"))) %>%
  as.data.frame()

# functions
stderr <- function(y) round(sd(y,na.rm=TRUE)/sqrt(length(y[!is.na(y)])), 2)
mm <- function(y) round(mean(y,na.rm=TRUE), 2)
x$key <- as.factor(x$key)
x$n_sim <- as.factor(x$n_sim)
x$var.type <- as.factor(x$var.type)
x$cv <- as.factor(x$cv)

x$key <- droplevels(x$key, c("robust", "comb", "unwght", "exact", NA, "ar1"))
x$cv <- droplevels(x$cv, c("FALSE", "kfold")) # "TRUE", 
x$var.type <- droplevels(x$var.type, c("TRUE", "fastBoot")) # "TRUE", 

# table
toLatex( tabular( n_sim * RowFactor(N_i,levelnames=unique(N_i))
                  ~ key*r*value*
                    PlusMinus(mm, stderr, digits=6, scipen = 999),
                  data = x))

x %>% dplyr::as_tibble() %>%
  group_by(N_i, n_sim, sigma_sq, key, cv, var.type, r) %>%
  summarize(mean(value, na.rm=TRUE)) %>% print(n = "Inf")

########################
# Joint CI Table
########################
# factors
dat <- do.call(rbind, ls_jci)
dat$cv <- as.factor(dat$cv)
dat$var.type <- as.factor(dat$var.type)
dat$n_sim <- as.factor(dat$n_sim)
dat$N_i <- as.factor(dat$N_i)
dat$rho.smooth <- as.factor(dat$rho.smooth)
dat$l <- as.factor(dat$l)

library(tables)
n_included <- nVec
Ni <- 5
metric.nm <- "jCI"
tune.type <- "fastkfold"
x=dat %>% tibble %>%  
  dplyr::mutate(key = sub(paste0("_", metric.nm), "", key)) %>% # remove rmse
  dplyr::filter(cv == tune.type,
                n_sim %in% n_included) %>%
  dplyr::mutate(key=str_replace_all(key, fixed("fGEE"), "One-Step"),
                key=str_replace_all(key, fixed("exact"), "GLS-Ex"),
                key=str_replace_all(key, fixed("unwght"), "GLS-Ind"),
                key=str_replace_all(key, fixed("Li"), "Marginal")
  ) %>%
  dplyr::mutate(key=fct_relevel(key, c("One-Step","GLS-Ex","GLS-Ind", "Marginal", "pffr"))) %>%
  as.data.frame()

# functions
stderr <- function(y) round(sd(y,na.rm=TRUE)/sqrt(length(y[!is.na(y)])), 2)
mm <- function(y) round(mean(y,na.rm=TRUE), 2)
x$key <- as.factor(x$key)
x$n_sim <- as.factor(x$n_sim)
x$var.type <- as.factor(x$var.type)
x$cv <- as.factor(x$cv)
x$r <- as.factor(x$r)

x$cv <- droplevels(x$cv, c("kfold","FALSE"))
x$key <- droplevels(x$key, c("robust", "comb", NA, "ar1", "NA", "Marginal")) # "unwght", "exact", 
x$var.type <- droplevels(x$var.type, c("TRUE", "fastBoot")) # "TRUE", 

# table
toLatex( tabular( n_sim * RowFactor(N_i,levelnames=unique(N_i))
                  ~ key*value*
                    PlusMinus(mm, stderr, digits=6, scipen = 999),
                  data = x))

x %>% dplyr::as_tibble() %>%
  group_by(N_i, n_sim, sigma_sq, key, cv, var.type, r) %>%
  summarize(mean(value, na.rm=TRUE)) %>% print(n = "Inf")

# ------------------------------------------
# works for multi-column
toLatex( tabular( key * Multicolumn(var.type, width=2,
                                         levelnames=unique(var.type))
                  ~ n_sim*cv*value*
                    PlusMinus(mm, stderr, digits=2, scipen = 999),
                  data = x))


########################
# Pointwise CI Table
########################
# factors
dat <- do.call(rbind, ls_pci)
dat$cv <- as.factor(dat$cv)
dat$var.type <- as.factor(dat$var.type)
dat$n_sim <- as.factor(dat$n_sim)
dat$N_i <- as.factor(dat$N_i)
dat$rho.smooth <- as.factor(dat$rho.smooth)
dat$l <- as.factor(dat$l)

library(tables)
n_included <- nVec
Ni <- 5
metric.nm <- "pCI"
tune.type <- "fastkfold"
x=dat %>% tibble %>%  
  dplyr::mutate(key = sub(paste0("_", metric.nm), "", key)) %>% # remove rmse
  dplyr::filter(cv == tune.type,
                sigma_sq == sig,
                n_sim %in% n_included) %>%
  dplyr::mutate(key=str_replace_all(key, fixed("fGEE"), "One-Step"),
                key=str_replace_all(key, fixed("exact"), "GLS-Ex"),
                key=str_replace_all(key, fixed("unwght"), "GLS-Ind"),
                key=str_replace_all(key, fixed("Li"), "Marginal")
  ) %>%
  dplyr::mutate(key=fct_relevel(key, c("One-Step","GLS-Ex","GLS-Ind", "Marginal", "pffr"))) %>%
  as.data.frame()

# functions
stderr <- function(y) round(sd(y,na.rm=TRUE)/sqrt(length(y[!is.na(y)])), 2)
mm <- function(y) round(mean(y,na.rm=TRUE), 2)
x$key <- as.factor(x$key)
x$n_sim <- as.factor(x$n_sim)
x$var.type <- as.factor(x$var.type)
x$cv <- as.factor(x$cv)
x$r <- as.factor(x$r)

x$cv <- droplevels(x$cv, c("kfold", "FALSE"))
x$key <- droplevels(x$key, c("robust", "comb", NA, "ar1", "NA")) # "unwght", "exact", 
x$var.type <- droplevels(x$var.type, c("TRUE", "fastBoot")) # "TRUE", 

# table
toLatex( tabular( n_sim * RowFactor(N_i,levelnames=unique(N_i))
                  ~ key*value*
                    PlusMinus(mm, stderr, digits=6, scipen = 999),
                  data = x))

x %>% dplyr::as_tibble() %>%
  group_by(N_i, n_sim, sigma_sq, key, cv, var.type, r) %>%
  summarize(mean(value, na.rm=TRUE)) %>% print(n = "Inf")

########################
# Time
########################
# factors
dat <- do.call(rbind, ls_time)
dat$cv <- as.factor(dat$cv)
dat$var.type <- as.factor(dat$var.type)
dat$n_sim <- as.factor(dat$n_sim)
dat$N_i <- as.factor(dat$N_i)
dat$rho.smooth <- as.factor(dat$rho.smooth)
dat$sigma_sq <- as.factor(dat$sigma_sq)
dat$r <- as.factor(dat$r)

library(tables)
n_included <- nVec
Ni <- 5
metric.nm <- "time"
tune.type <- "fastkfold"
sig <- 10
x=dat %>% tibble %>%  
  dplyr::mutate(key = sub(paste0("_", metric.nm), "", key)) %>% # remove rmse
  dplyr::filter(cv == tune.type,
                sigma_sq == sig,
                n_sim %in% n_included) %>%
  dplyr::mutate(key=str_replace_all(key, fixed("fGEE"), "One-Step"),
                key=str_replace_all(key, fixed("exact"), "GLS-Ex"),
                key=str_replace_all(key, fixed("unwght"), "GLS-Ind"),
                key=str_replace_all(key, fixed("Li"), "Marginal")
  ) %>%
  dplyr::mutate(key=fct_relevel(key, c("One-Step","GLS-Ex","GLS-Ind", "Marginal", "pffr"))) %>%
  as.data.frame()

# functions
stderr <- function(y) round(sd(y,na.rm=TRUE)/sqrt(length(y[!is.na(y)])), 2)
mm <- function(y) round(mean(y,na.rm=TRUE), 2)
x$key <- as.factor(x$key)
x$n_sim <- as.factor(x$n_sim)
x$var.type <- as.factor(x$var.type)
x$cv <- as.factor(x$cv)

stderr <- function(y) round(sd(y,na.rm=TRUE)/sqrt(length(y[!is.na(y)])), 2)
mm <- function(y) round(mean(y,na.rm=TRUE), 2)
x$key <- as.factor(x$key)
x$n_sim <- as.factor(x$n_sim)
x$var.type <- as.factor(x$var.type)
x$cv <- as.factor(x$cv)

x$key <- droplevels(x$key, c("robust", "comb", "unwght", "exact", NA, "ar1"))
x$cv <- droplevels(x$cv, c("FALSE", "kfold")) # "TRUE", 
x$var.type <- droplevels(x$var.type, c("TRUE", "fastBoot")) # "TRUE", 

# table
toLatex( tabular( n_sim * RowFactor(N_i,levelnames=unique(N_i))
                  ~ key*value*
                    PlusMinus(mm, stderr, digits=6, scipen = 999),
                  data = x))

x %>% dplyr::as_tibble() %>%
  group_by(N_i, n_sim, sigma_sq, key, r,cv) %>%
  summarize(mean(value, na.rm=TRUE)) %>% print(n = "Inf")

# ------------------------------------------
# works for multi-column
toLatex( tabular( key * Multicolumn(var.type, width=2,
                                    levelnames=unique(var.type))
                  ~ n_sim*cv*value*
                    PlusMinus(mm, stderr, digits=2, scipen = 999),
                  data = x))

#################################################################################
# **** extra tables not included in paper ****
#################################################################################
###########################
# Pointwise Slope CI Table
###########################
# factors
dat <- do.call(rbind, ls_pci_slp)
dat$cv <- as.factor(dat$cv)
dat$var.type <- as.factor(dat$var.type)
dat$n_sim <- as.factor(dat$n_sim)
dat$N_i <- as.factor(dat$N_i)
dat$rho.smooth <- as.factor(dat$rho.smooth)
dat$l <- as.factor(dat$l)

library(tables)
n_included <- nVec
Ni <- 5
metric.nm <- "pci"
tune.type <- "fastkfold"
x=dat %>% tibble %>%  
  dplyr::mutate(key = sub(paste0("_", metric.nm), "", key)) %>% # remove rmse
  dplyr::filter(cv == tune.type,
                var.type %in% var.vec,
                n_sim %in% n_included) %>%
  as.data.frame()

# functions
stderr <- function(y) round(sd(y,na.rm=TRUE)/sqrt(length(y[!is.na(y)])), 2)
mm <- function(y) round(mean(y,na.rm=TRUE), 2)
x$key <- as.factor(x$key)
x$n_sim <- as.factor(x$n_sim)
x$var.type <- as.factor(x$var.type)
x$cv <- as.factor(x$cv)


# functions
stderr <- function(y) round(sd(y,na.rm=TRUE)/sqrt(length(y[!is.na(y)])), 2)
mm <- function(y) round(mean(y,na.rm=TRUE), 2)
x$key <- as.factor(x$key)
x$n_sim <- as.factor(x$n_sim)
x$var.type <- as.factor(x$var.type)
x$cv <- as.factor(x$cv)
x$r <- as.factor(x$r)

x$cv <- droplevels(x$cv, "kfold")
x$key <- droplevels(x$key, c("robust", "comb", "unwght", "exact", NA, "ar1"))
x$cv <- droplevels(x$cv, c("FALSE")) # "TRUE", 
x$var.type <- droplevels(x$var.type, c("TRUE", "fastBoot")) # "TRUE", 

# table
toLatex( tabular( n_sim * RowFactor(N_i,levelnames=unique(N_i))
                  ~ key*r*value*
                    PlusMinus(mm, stderr, digits=6, scipen = 999),
                  data = x))

x %>% dplyr::as_tibble() %>%
  group_by(N_i, n_sim, sigma_sq, key, cv, var.type, r) %>%
  summarize(mean(value, na.rm=TRUE)) %>% print(n = "Inf")

# ------------------------------------------
# works for multi-column
toLatex( tabular( key * Multicolumn(var.type, width=2,
                                    levelnames=unique(var.type))
                  ~ n_sim*cv*value*
                    PlusMinus(mm, stderr, digits=2, scipen = 999),
                  data = x))


########################
# RMSE CI Table
########################
# factors
dat <- do.call(rbind, ls_rmse)
dat$cv <- as.factor(dat$cv)
dat$var.type <- as.factor(dat$var.type)
dat$n_sim <- as.factor(dat$n_sim)
dat$N_i <- as.factor(dat$N_i)
dat$rho.smooth <- as.factor(dat$rho.smooth)
dat$sigma_sq <- as.factor(dat$sigma_sq)
dat$l <- as.factor(dat$l)
dat$r <- as.factor(dat$r)


library(tables)
n_included <- nVec
Ni <- 5
metric.nm <- "rmse"
tune.type <- "fastkfold"
sig <- 10
x=dat %>% tibble %>%  
  dplyr::mutate(key = sub(paste0("_", metric.nm), "", key)) %>% # remove rmse
  dplyr::filter(cv == tune.type,
                sigma_sq == sig,
                n_sim %in% n_included) %>%
  as.data.frame()

# functions
stderr <- function(y) round(sd(y,na.rm=TRUE)/sqrt(length(y[!is.na(y)])), 2)
mm <- function(y) round(mean(y,na.rm=TRUE), 2)
x$key <- as.factor(x$key)
x$n_sim <- as.factor(x$n_sim)
x$var.type <- as.factor(x$var.type)
x$cv <- as.factor(x$cv)

x$cv <- droplevels(x$cv, "kfold")
x$key <- droplevels(x$key, c("robust", "comb", "unwght", "exact", NA))
x$cv <- droplevels(x$cv, c("FALSE")) # "TRUE", 
x$var.type <- droplevels(x$var.type, c("TRUE", "fastBoot")) # "TRUE", 


# table
toLatex( tabular( n_sim * RowFactor(N_i,levelnames=unique(N_i))
                  ~ key*cv*value*
                    PlusMinus(mm, stderr, digits=6, scipen = 999),
                  data = x))

x %>% dplyr::as_tibble() %>%
  group_by(N_i, n_sim, sigma_sq, cv, l, r, key) %>%
  filter(var.type == FALSE) %>%
  summarize(mean(value, na.rm=TRUE)) %>% print(n = "Inf")

########################
# RMSE Table SLOPE
########################
# factors
dat <- do.call(rbind, ls_rmse_slp)
dat$cv <- as.factor(dat$cv)
dat$var.type <- as.factor(dat$var.type)
dat$n_sim <- as.factor(dat$n_sim)
dat$N_i <- as.factor(dat$N_i)
dat$rho.smooth <- as.factor(dat$rho.smooth)
dat$sigma_sq <- as.factor(dat$sigma_sq)
dat$l <- as.factor(dat$l)
dat$r <- as.factor(dat$r)

library(tables)
n_included <- nVec
Ni <- 5
metric.nm <- "rmse"
tune.type <- "fastkfold"
sig <- 10
x=dat %>% tibble %>%  
  dplyr::mutate(key = sub(paste0("_", metric.nm), "", key)) %>% # remove rmse
  dplyr::filter(cv == tune.type,
                #var.type %in% var.vec,
                sigma_sq == sig,
                n_sim %in% n_included) %>%
  as.data.frame()

# functions
stderr <- function(y) round(sd(y,na.rm=TRUE)/sqrt(length(y[!is.na(y)])), 2)
mm <- function(y) round(mean(y,na.rm=TRUE), 2)
x$key <- as.factor(x$key)
x$n_sim <- as.factor(x$n_sim)
x$var.type <- as.factor(x$var.type)
x$cv <- as.factor(x$cv)

x$cv <- droplevels(x$cv, "kfold")
x$key <- droplevels(x$key, "robust")
x$cv <- droplevels(x$cv, c("TRUE", "FALSE"))

# table
toLatex( tabular( n_sim * RowFactor(N_i,levelnames=unique(N_i))
                  ~ key*cv*value*
                    PlusMinus(mm, stderr, digits=4, scipen = 999),
                  data = x))

x %>% dplyr::as_tibble() %>%
  group_by(N_i, n_sim, sigma_sq, key, cv, l, r) %>%
  summarize(mean(value, na.rm=TRUE)) %>% print(n = "Inf")
