library(MASS)
GenerateData <- function(Nsubj=100, numFunctPoints = 101, min_visit=8, max_visit=12,
                         numLongiPoints = 41, sigma_sq = 1.5, 
                         sigma_z11 = 3, sigma_z12 = 1.5, 
                         sigma_z21 = 2, sigma_z22 = 1, rho = 0.75,
                         rho.func = 0.9,
                         divFactor = 10, # divide betas by this for binomial to avoid all 1s
                         corstr = c("unspecified", "exchangeable", "independent", "ar1"),
                         family = "gaussian",
                         true.cov = FALSE,
                         rho.vec){
  # Arguments:
  # Nsubj -- number of subjects, default to 100
  # numFunctPoints - number of points at which functions are observed
  # min_visit = 8 #(min # visits is 8) - minimum number of visits for each subject
  # max_visit = 12 #(max # visits is 12) - maximum number of visits for each subject
  # numLongiPoints -- number of points to evaluate inner eigenfunctions (for visits) on 
  # inner (visits) eigenfunction scores
  # sigma_sq -- white noise variance defaults to 1.5 for SNR = 5 (see Ana-Maria paper)
  # sigma_z11 <- 3
  # sigma_z12 <- 1.5
  # sigma_z21 <- 2
  # sigma_z22 <- 1
  ############################################################
  true.Cov <- NULL
  psi = NULL 
  xi = NULL
  zeta = NULL
  sigma = NULL
  time.vec <- NULL
  # time points at which functional responses are collected
  s <- seq(0, 0.99, length = numFunctPoints)
  
  ########################################
  # Select time points for visits j
  ########################################
  # vector of visit argvalues T of length numLongiPoints
  T <- seq(0, 1, len = numLongiPoints) 
  
  # select number of visits per subject from uniform distribution
  runifdisc <- function(n, min=0, max=1) {sample(min:max, n, replace=TRUE)}
  if(min_visit == max_visit){
    m_i <- rep(min_visit, Nsubj)
  }else{
    m_i <- runifdisc(n=Nsubj, min = min_visit, max = max_visit)
  }
  
  # vector of subject visits
  subjID <- rep(1:Nsubj, m_i)
  n <- length(subjID)
  
  # select for each subject i from all possible argvalues T, subject visits of length m_i
  replace.ind <- I(max(m_i) > Nsubj) # GCL added: use to draw Tij with replacement
  Tij <- lapply(m_i, function(x){sort(sample(T, size = x, replace = replace.ind))}) # GCL changed replace = TRUE
  # keep indicators for the visit number selected
  visitID <- lapply(Tij, function(x) which(T %in% x))
  # convert into a vector
  Tij <- unlist(Tij)
  visitID <- unlist(visitID)
  # create column id names 
  times <- paste("t_", 1:length(s), sep = "")
  # create row id names
  visits <- paste("S_", subjID, "v_", visitID , sep="")
  
  ########################################
  # Define k=2 outer basis functions 
  ########################################
  b_per <- 2
  phi_1 = function(s) {rep(1, length(s))}
  phi_2 = function(s, b_per=2) {sqrt(2) * sin (b_per*pi*s)}
  # phi values
  phi <- list()
  phi$phi_1 <- phi_1(s)
  phi$phi_2 <- phi_2(s)
  
  ########################################
  # Define covariates
  ########################################
  # covariate 1
  Cov1.pred <- rnorm(Nsubj)

  Cov1 <- rep(Cov1.pred, m_i)
  a <- 1
  pho <- 0.7
  
  # covariate 2
  ar.sim <- unlist(lapply(m_i,function(b){arima.sim(model=list(ar=pho), n=b)}))
  Cov2 <- a*Tij + ar.sim
  Cov <- cbind(Cov1,Cov2)
  rownames(Cov) <- visits
  colnames(Cov) <- c("X1", "X2")
  
  ########################################
  # Define mean response functions
  ########################################
  #  exchangeable
  if(corstr == "exchangeable" & family == "gaussian"){
    beta0 <- 2 + sin(pi*s) + sqrt(2)*cos(3*pi*s) + 1
    beta1 <- 2 + cos(2*pi*s) + sqrt(2)*cos(3*pi*s) + 1
    beta2 <- 1/60*dnorm( (s-0.2) / 0.1^2 ) +
      1/200*dnorm( (s-0.35) / 0.1^2 ) -
      1/250*dnorm( (s-0.65) / 0.06^2 ) +
      1/60*dnorm( (s-1) / 0.07^2 )
  }else if(family == "binomial" | corstr == "ar1"){
    beta0 <- 2 + sin(pi*s) + sqrt(2)*cos(3*pi*s) + 1
    beta1 <- 2 + cos(2*pi*s) + sqrt(2)*cos(3*pi*s) + 1
    beta2 <- 5*(dnorm( (s-0.35) / 0.1 ) -
                 dnorm( (s-0.65) / 0.2 ))

    if(family == "binomial"){
      beta0 <- beta0 / divFactor
      beta1 <- beta1 / divFactor 
      beta2 <- beta2 / divFactor
    }
  }

  ############################################################
  # Define mean function at each time (s) and visit (T) point 
  #############################################################
  mean <- Cov%*%t(cbind(beta1,beta2)) + matrix(rep(beta0,n), nr=n, byrow = TRUE) 
  rownames(mean) <- visits
  colnames(mean) <- times
  if(corstr == "unspecified"){
    ########################################################
    # Define time varying components
    # xi_ik(T) = zeta_ik1 psi_k1(T) + zeta_ik2 psi_k2(T)
    ########################################################
    
    #define psi_k1 and psi_k2
    psi_11 <- function(T) {sqrt(2) * cos (2*pi*T)}
    psi_12 <- function(T) {sqrt(2) * sin (2*pi*T)}
    psi_21 <- function(T) {sqrt(2) * cos (4*pi*T)}
    psi_22 <- function(T) {sqrt(2) * sin (4*pi*T)}
    
    # calculate eigenfunctions values at each point
    psi <- list()
    psi$psi_11 <- psi_11(T)
    psi$psi_12 <- psi_12(T) 
    psi$psi_21 <- psi_21(T)
    psi$psi_22 <- psi_22(T) 
    
    # define zeta_k1 and zeta_k2
    zeta_11 <- rnorm(n = Nsubj, mean = 0, sd = sqrt(sigma_z11-1))
    zeta_12 <- rnorm(n = Nsubj, mean = 0, sd = sqrt(sigma_z12))
    zeta_21 <- rnorm(n = Nsubj, mean = 0, sd = sqrt(sigma_z21-0.5))
    zeta_22 <- rnorm(n = Nsubj, mean = 0, sd = sqrt(sigma_z22))
    # zeta values that will be used for evaluation, inner eigenfunctions scores
    zeta <- list()
    zeta$zeta_11 <- zeta_11
    zeta$zeta_12 <- zeta_12
    zeta$zeta_21 <- zeta_21
    zeta$zeta_22 <- zeta_22
    
    # create xi's
    xi_1 <- rep(zeta_11, m_i)*psi$psi_11[visitID] + rep(zeta_12,m_i)*psi$psi_12[visitID]
    names(xi_1) <- visits
    xi_2 <- rep(zeta_21, m_i)*psi$psi_21[visitID] + rep(zeta_22,m_i)*psi$psi_22[visitID]
    names(xi_2) <- visits
    # xi values that will be used for evaluation
    xi <-list()
    xi$xi_1 <- xi_1
    xi$xi_2 <- xi_2
    
    ########################################
    # define X_i(s, T_ij)
    # X_i(s, T_ij) = xi_1*phi_1 + xi_2*phi_2
    ########################################
    X <- xi_1%*%matrix(phi_1(s), nrow=1) + xi_2%*%matrix(phi_2(s), nrow=1)
    rownames(X) <- visits
    colnames(X) <- times
    
    ########################################
    # Generate random error terms
    # there are (i=length(s))*(j=sum(m_i)) = #time points * #visits for all subjects
    ########################################
    # noise/signal = 1:5
    inner.noise_1 <- rnorm(n = n, mean = 0, sd = sqrt(1))
    inner.noise_2 <- rnorm(n = n, mean = 0, sd = sqrt(0.5))
    epsilon.inner <- (inner.noise_1)%*%matrix(phi_1(s), nrow=1) 
    + (inner.noise_2)%*%matrix(phi_2(s), nrow=1)
    rownames(epsilon.inner) <- visits
    colnames(epsilon.inner) <- times
    
    # outer noise
    epsilon <- matrix(rnorm(n*numFunctPoints, mean = 0, sd = sqrt(sigma_sq)), 
                      nrow = n, ncol = numFunctPoints)
    # assign row and column names to error terms matrix epsilon_ij
    rownames(epsilon) <- visits
    colnames(epsilon) <- times
    
    ########################################
    # get data Y_ij(s) = mu(s,T_ij) + X_i(s, T_ij) + epsilon.inner_ij(s) + epsilon_ij(s)
    ########################################
    Y = mean + X + epsilon.inner + epsilon 
    Y.star = mean + X 
  }#end if corstr == unspecified
  
  else if (corstr == "exchangeable"){
    if(family == "gaussian"){
      # k=1
      # subject specific random effects
      xi_i1_coef <- rnorm(n = Nsubj, mean = 0, sd = sqrt(sigma_z11))
      names(xi_i1_coef) <- paste("S_", 1:Nsubj, sep="")
      xi_i1 <- rep(xi_i1_coef, m_i)
      names(xi_i1) <- visits
      # subject-visit specific effects
      xi_ij1 <- rnorm(sum(m_i), mean =0, sd = sqrt(sigma_z12))
      names(xi_ij1) <- visits
      # k=2
      # define subject-visit specific effects
      xi_i2_coef <- rnorm(n = Nsubj, mean = 0, sd = sqrt(sigma_z21))
      names(xi_i2_coef) <- paste("S_", 1:Nsubj, sep="")
      xi_i2 <- rep(xi_i2_coef, m_i)
      names(xi_i2) <- visits
      # subject-visit specific effects
      xi_ij2 <- rnorm(sum(m_i), mean =0, sd = sqrt(sigma_z22))  
      names(xi_ij2) <- visits
      
      # xi values that will be used for evaluation
      # subject scores
      xi <-list()
      xi$xi_i1 <- xi_i1_coef
      xi$xi_i2 <- xi_i2_coef
      # subject-visit specific
      xi$xi_ij1 <- xi_ij1
      xi$xi_ij2 <- xi_ij2
      
      # define X_i(s, T_ij)
      # subject specific random process U_i
      U_i = t(matrix(phi_1(s))%*%xi_i1) + t(matrix(phi_2(s))%*%xi_i2)
      rownames(U_i) <- visits
      colnames(U_i) <- times
      browser()
      # subject-visit specific random process V_ij
      V_ij <-  t(matrix(phi_1(s))%*%xi_ij1) + t(matrix(phi_2(s))%*%xi_ij2)
      rownames(V_ij) <- visits
      colnames(V_ij) <- times
      # response X = U_i + V_ij 
      X = U_i + V_ij 
      ########################################
      # Generate random error terms
      # there are (i=length(s))*(j=sum(m_i)) = #time points * #visits for all subjects
      ########################################
      epsilon <- matrix(rnorm(n*numFunctPoints, mean = 0, sd = sqrt(sigma_sq)), 
                        nrow = n, ncol = numFunctPoints)
      # assign row and column names to error terms matrix epsilon_ij
      rownames(epsilon) <- visits
      colnames(epsilon) <- times
      
      ########################################
      # combine to get data Y_ij(s) = mu(s,T_ij) + X_i(s, T_ij) + epsilon_ij(s)
      ########################################
      Y = mean + X + epsilon #response function used for evaluation (matrix form)
      Y.star = mean + X #used for accuracy evaluation
    #end if corstr == exchangeable
  }else{
    # binomial exchangeable
    
    ###################################
    # ar1 error longitudinal direction
    ###################################
    message("binary exchangeable")
    beta.mat <- t(cbind(beta0, beta1,beta2))
    XX <- cbind(1, Cov)
    m_i2 <- c(0, m_i)

    ar.error <- lapply(1:Nsubj, function(ii)
      lapply(1:L, function(s){
        X.mat <- XX[(sum(m_i2[1:ii])+1):sum(m_i[1:ii]),];
        rbin(clsize = m_i[ii], 
             intercepts = 0,
             betas = beta.mat[,s], 
             xformula = ~X.mat, 
             xdata = X.mat,
             cor.matrix = stats::toeplitz(c(1, rep( rho.vec[s], m_i[ii]-1)) ) ,
             link = "logit")$simdata[,c("y")] # "time",
      } 
      )
    )
    
    ar.error <- lapply(ar.error, function(xx) do.call(cbind, xx))
    ar.error <- lapply(ar.error, function(xx) cbind(1:nrow(xx), xx)) # add in time
    ar.error <- do.call(rbind, ar.error)
    time.vec <- ar.error[,1]
    ar.error <- ar.error[,-1] # remove time
    
    Y <- Y.star <- ar.error 
    X <- NULL
    
    if(true.cov){
      #---------------------------------
      # estimate true covariance matrix
      #---------------------------------
      rownames(Y) <- visits
      colnames(Y) <- times
      
      num.extract <- function(string){ as.numeric(str_extract(string, "(?<=S_)\\d+(?=v)"))}
      ids <- sapply(rownames(Y), function(xx) num.extract(xx))
      rows <- lapply(unique(ids), function(ii) which(ids == ii))
      subj.X <- lapply(rows, function(ii) as.vector((Y[ii,]) )) # X[ii, s.idx] vectorize subject-specific matrices of X
      XX <- do.call(rbind, subj.X)  
      true.Cov <- cor(XX)
    }else{
      true.Cov <- NULL
    }
  } # end binomial exchangeable
  }else if (corstr == "independent"){
    # k=1
    # subject-visit specific effects
    xi_ij1 <- rnorm(sum(m_i), mean =0, sd = sqrt(sigma_z11+sigma_z12))
    names(xi_ij1) <- visits
    # k=2
    # subject-visit specific effects
    xi_ij2 <- rnorm(sum(m_i), mean =0, sd = sqrt(sigma_z21+sigma_z22))  
    names(xi_ij2) <- visits
    
    # xi values that will be used for evaluation
    # subject-visit specific scores
    xi <-list()
    xi$xi_ij1 <- xi_ij1
    xi$xi_ij2 <- xi_ij2
    
    # define X_i(s, T_ij)
    # subject-visit specific random process V_ij
    V_ij <-  t(matrix(phi_1(s))%*%xi_ij1) + t(matrix(phi_2(s))%*%xi_ij2)
    rownames(V_ij) <- visits
    colnames(V_ij) <- times
    # response X = V_ij 
    X = V_ij 
    ########################################
    # Generate random error terms
    # there are (i=length(s))*(j=sum(m_i)) = #time points * #visits for all subjects
    ########################################
    epsilon <- matrix(rnorm(n*numFunctPoints, mean = 0, sd = sqrt(sigma_sq)), 
                      nrow = n, ncol = numFunctPoints)
    # assign row and column names to error terms matrix epsilon_ij
    rownames(epsilon) <- visits
    colnames(epsilon) <- times
    
    ########################################
    # combine to get data Y_ij(s) = mu(s,T_ij) + X_i(s, T_ij) + epsilon_ij(s)
    ########################################
    Y = mean + X + epsilon #response function used for evaluation (matrix form)
    Y.star = mean + X #used for accuracy evaluation
  }#end if corstr == independent
  
  else if(corstr == "ar1"){
    L <- numFunctPoints
    # Radial basis function kernel and i.i.d. multivariate Normal draws
    if(family == "gaussian"){
      rbfkernel <- function(x, y, sigma = sigma_sq, l = 20) {
        sigma * exp(-(x - y)^2 / (2 * l^2))
      }
      
      fun_cov_sig <- matrix(0, nrow = L, ncol = L)
      
      for (i in 1:L) {
        for (j in i:L) {
          fun_cov_sig[i, j] <- fun_cov_sig[j, i] <- rbfkernel(i, j)
        }
      }
      
      # functional domain
      X <- MASS::mvrnorm(n = sum(m_i), mu = rep(0,L), Sigma = fun_cov_sig)
      
      ###################################
      # ar1 error longitudinal direction
      ###################################
      ar.error <- lapply(m_i, function(ii) 
        lapply(1:L, function(s)
          SuperGauss::rnormtz(1,
                              acf = sigma_sq*(rho.vec[s])^(0:(ii-1))))

      ) 
      ar.error <- lapply(ar.error, function(xx) do.call(cbind, xx))
      ar.error <- lapply(ar.error, function(xx) cbind(1:nrow(xx), xx)) # add in time
      ar.error <- do.call(rbind, ar.error)
      time.vec <- ar.error[,1]
      ar.error <- ar.error[,-1] # remove time
      
      #################################
      # ar1 error functional direction
      #################################
      ar.error.fun <- SuperGauss::rnormtz(sum(m_i),
                                          acf = sigma_sq*(rho.func)^(0:(L-1)))

      
      ########################################
      # combine to get data Y_ij(s) = mu(s,T_ij) + X_i(s, T_ij) + epsilon_ij(s)
      ########################################
      Y = mean + ar.error #response function used for evaluation (matrix form)
      Y.star = Y # no extra noise added on
      
    }
    else if(family == "binomial"){
      
      ###################################
      # ar1 error longitudinal direction
      ###################################
      message("binary")
      beta.mat <- t(cbind(beta0, beta1,beta2))
      XX <- cbind(1, Cov)
      m_i2 <- c(0, m_i)
      
      ar.error <- lapply(1:Nsubj, function(ii)
        lapply(1:L, function(s){
          X.mat <- XX[(sum(m_i2[1:ii])+1):sum(m_i[1:ii]),];
          rbin(clsize = m_i[ii], 
               intercepts = 0,
               betas = beta.mat[,s], 
               xformula = ~X.mat, 
               xdata = X.mat,
               cor.matrix = stats::toeplitz((rho.vec[s])^(0:(m_i[ii]-1))),
               link = "logit")$simdata[,c("y")] # "time",
        } 
        )
      )
      
      ar.error <- lapply(ar.error, function(xx) do.call(cbind, xx))
      ar.error <- lapply(ar.error, function(xx) cbind(1:nrow(xx), xx)) # add in time
      ar.error <- do.call(rbind, ar.error)
      time.vec <- ar.error[,1]
      ar.error <- ar.error[,-1] # remove time

      Y <- Y.star <- ar.error 
      X <- NULL
      
      if(true.cov){
        #---------------------------------
        # estimate true covariance matrix
        #---------------------------------
        rownames(Y) <- visits
        colnames(Y) <- times
        
        num.extract <- function(string){ as.numeric(str_extract(string, "(?<=S_)\\d+(?=v)"))}
        ids <- sapply(rownames(Y), function(xx) num.extract(xx))
        rows <- lapply(unique(ids), function(ii) which(ids == ii))
        subj.X <- lapply(rows, function(ii) as.vector((Y[ii,]) )) # X[ii, s.idx] vectorize subject-specific matrices of X
        XX <- do.call(rbind, subj.X)  
        true.Cov <- cor(XX)
      }else{
        true.Cov <- NULL
      }
    }
  }
  
  return(list(data=list(subjID=subjID, Tij=Tij, visitID=visitID, funcArg=s, Y=Y, 
                        Cov=Cov), X=X, Y.star=Y.star, mean=mean, beta=cbind(beta0,beta1,beta2), 
              phi=phi, xi=xi, psi=psi, zeta=zeta, time = time.vec, true.Cov = true.Cov))
  
}#end of function GenerateData