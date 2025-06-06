# pre-process SL data
library(R.matlab)
library(tidyverse)
library(ggplot2)
library(refund)
library(dplyr)

# read data and extract variable names
wd <- "/Users/Desktop/Research/functional_GEE/data/SL/Chronic/matlab/pooled/v6"
setwd(wd)
files <- list.files()
n <- length(files) # number of animals

# raw data
hz <- 1 / 0.03373076
min.len <- 4495 # min length of bins that all trials from all animals have

# target
target.sec <- 10 # number of sec of data
target.len <- 150 #ceiling(target.sec*hz)

# names
matdat <- R.matlab::readMat(files[1])
trial <- 2 # example trial number
nm <- rownames(matdat$fTrial[trial,,][[1]])

extract_number <- function(input_string) {
  # Use regex to find the number between "Mouse" and "_"
  match <- str_extract(input_string, "(?<=Mouse)\\d+(?=_)")
  return(as.numeric(match))
}

animal.list <- vector(length = n, "list")
for(j in 1:n){
  print(paste(j, "of", n))
  
  # read data
  id <- extract_number(files[j])
  matdat <- R.matlab::readMat(files[j])
  
  num.trials <- dim(matdat$fTrial)[1]
  cell.ind <- I(matdat$iscell[,1] == 1) # channel is a cell
  num.trials <- dim(matdat$fTrial)[1]
  
  # iterate over trials
  t.list <- vector(length = num.trials, "list")
  for(tt in 1:num.trials){
    
    # speed to find time
    trg <- "fSpeed"
    r <- which(nm == trg)
    sp <- matdat$fTrial[tt,,][[1]][r,,][[1]]
    speed.sd <- sd(sp)
    speed.idx <- which(sp > 2 * speed.sd)[1]
    if(speed.idx + target.len > min.len){
      message("over by ": speed.idx + target.len - min.len)
      speed.idx <- min.len - target.len
    }
   
    
    # spikes
    trg <- "fSpks"
    r <- which(nm == trg)
    spikes <- matdat$fTrial[tt,,][[1]][r,,][[1]]
    spikes <- 1*(spikes[speed.idx:(speed.idx + target.len-1), cell.ind]>0) # rows are length of time, columns filter non-cells
    spikes <- t(spikes)
    colnames(spikes) <- paste0("s_", 1:target.len)
    #spikes <- cbind(1:sum(cell.ind), spikes) # add neuron number
    
    # whiskers
    trg <- "fWhisk1"
    r <- which(nm == trg)
    whisk <- matdat$fTrial[tt,,][[1]][r,,][[1]]
    whisk <- whisk[speed.idx:(speed.idx + target.len-1), ] # rows are length of time
    whisk <- matrix(whisk, nrow = 1)
    colnames(whisk) <- paste0("w_", 1:target.len)
    
    # speed
    trg <- "fSpeed"
    r <- which(nm == trg)
    spd <- matdat$fTrial[tt,,][[1]][r,,][[1]]
    spd <- spd[speed.idx:(speed.idx + target.len-1), ] # rows are length of time
    spd <- matrix(spd, nrow = 1)
    colnames(spd) <- paste0("spd_", 1:target.len)
    
    t.list[[tt]] <- data.frame(id = id,
                               trial = tt,
                               neuron = 1:sum(cell.ind),
                               spikes = spikes,
                               speed = spd,
                               whisk = whisk)
  }
  
  animal.list[[j]] <- do.call(rbind, t.list) # pool trials for animal j
  
}


# write data
df <- do.call(rbind, animal.list) # pool across animals
rm(animal.list, t.list, matdat)

setwd("/Users/Desktop/Research/functional_GEE/data/SL/Chronic/matlab/pooled/pre_process")
data.table::fwrite(df, "SL_s1_speed.csv", row.names = FALSE)
