# Install required packages
.libPaths()
user_lib <- Sys.getenv("R_LIBS_USER")

if (!dir.exists(user_lib)) {
  dir.create(user_lib, recursive = TRUE)
}
install.packages("MCMCprecision", lib = user_lib, dependencies = TRUE)
install.packages("MASS", lib = user_lib, dependencies = TRUE)
library(MCMCprecision, lib.loc = user_lib)
library(MASS, lib.loc = user_lib)

# Prepare Data
df <- read.table(file='wine.data', sep=",", header=FALSE)
df$V1 = NULL
X <- as.matrix(df) # convert data into numeric argument

# Set K, set d <- # of dimensions and N <- # of points
K <- 3
d <- dim(X)[2]
N <- dim(X)[1]

# Set Hyperparameters
kappa_0 <- 1000
alpha_0 <- 1
Beta_0 <- 1
mu_0 <- 0

# Set the number of iteration the Gibbs sampler
# is allowed to attempt to find an improvement and
# EM to converge.
R.Gibbs = 100
R.EM = 1000


######################
#####Functions########
######################


log_dMVN <- function(point, mu, Sigma){
  # Calculates log likelihood of a single point under a single
  # Component of out model. Used in "point_likelihood" function
  # and "log_z_cond" function.
  log_density <- sum(dnorm(point, mean = mu, sd = sqrt(Sigma), log = TRUE))
  log_density
} 

log_sum_exp <- function(x){
  # Function for log-sum-exp trick to avoid underflow
  # when normalizing probabilities.
  # https://gregorygundersen.com/blog/2020/02/09/log-sum-exp/
  c <- max(x)
  c + log(sum(exp(x - c)))
}

########Objective Function#######
Obj <- function(data, mu, Sigma, pi.mix){
  # Calculates overall log likelihood, using all data.
  sum(apply(data, 1, point_likelihood, mu<-mu, Sigma<-Sigma, pi.mix<-pi.mix))
}

calc.BIC <- function(obj, mu, Sigma, pi.mix){
  log_likelihood = obj
  C <- length(mu) + length(Sigma) + length(pi.mix) - 1
  # minus 1 because pi.mix must sum to 1 -> only K - 1 are free
  BIC <- C*log(N) - 2*log_likelihood
  BIC
}

point_likelihood <- function(point, mu, Sigma, pi.mix){
  # Log likelihood of a single point under GMM, used in function "Obj".
  likelihood <- 0
  for (k in 1:K){
    likelihood <- likelihood + pi.mix[k]*exp(log_dMVN(point, mu[k,], Sigma[k,]))
  }  
  log(likelihood)
}

log_z_cond <- function(k, point, mu, Sigma, pi.mix){
  # P(z_i = k | x_i, theta). 
  # Calculates omega_ik given x_i and model parameters.
  log_numerator <- log(pi.mix[k]) + log_dMVN(point, mu[k,], Sigma[k,])
  mvn_densities <- sapply(1:K, function(i) log_dMVN(point, mu[i,], Sigma[i,]))
  log_denom_terms <- log(pi.mix) + mvn_densities
  # Use log_sum_exp trick
  log_denom <- log_sum_exp(log_denom_terms)
  log_out <- log_numerator - log_denom
  out <- exp(log_out)
}

sample_pi <- function(z){
  # For sampling pi in Gibbs sampler.
  # Samples from Dirichlet(alpha | others)
  alpha_k.vec <- numeric(K)
  for (i in 1:K){
    alpha_k.vec[i] <- alpha_0 + sum(z == i)
  }
  pi_samples <- rdirichlet(1, alpha_k.vec)
  pi_samples
}


sample_z <- function(data, mu, Sigma, pi.mix){
  # For sampling Z in Gibbs sampler.
  # p(z | others)
  samples <- numeric(length(data[,1]))
  mu <- mu[order(mu[,1], decreasing=F),] # Attempt to avoid the label-switching problem.
  for (i in 1:length(data[,1])){
    probs <- numeric(K)
    for (k in 1:K){
      probs[k] <- log_z_cond(k, data[i,], mu, Sigma, pi.mix)
    }
    samples[i] <- sample(1:K, 1, replace=FALSE, probs)
  }
  samples
}

sample_means <- function(X, Sigma, z){
  # For sampling mu in Gibbs sampler.
  # p(mu | others).
  samples <- matrix(0, K, d) 
  for (k in 1:K){
    indexes <- which(z == k, arr.ind=TRUE)
    n_k <- sum(z == k)
    if (n_k==1){
      data_k <- X[indexes,]
      sample_mean <- data_k
      kappa.hat <- 1/(kappa_0^(-1) + n_k)
      mu.hat <- kappa.hat*((kappa_0^(-1))*mu_0 + n_k*sample_mean)
      Sigma <- diag(kappa.hat*Sigma[k,], d, d)
      samples[k,] <- mvrnorm(n=1, mu=mu.hat, Sigma)
    }
    else if (n_k==0){
      sample_mean <- rep(mu_0, d)
      kappa.hat <- 1/(kappa_0^(-1) + n_k)
      mu.hat <- kappa.hat*((kappa_0^(-1))*mu_0 + n_k*sample_mean)
      Sigma <- diag(kappa.hat*Sigma[k,], d, d)
      samples[k,] <- mvrnorm(n=1, mu=mu.hat, Sigma)
    }
    else{
      data_k <- X[indexes,]
      sample_mean <- colSums(data_k)/(dim(data_k)[1])
      kappa.hat <- 1/(kappa_0^(-1) + n_k)
      mu.hat <- kappa.hat*((kappa_0^(-1))*mu_0 + n_k*sample_mean)
      Sigma <- diag(kappa.hat*Sigma[k,], d, d) # Note to Dr. Goh: I am not sure what we are doing with
      samples[k,] <- mvrnorm(n=1, mu=mu.hat, Sigma) # Sigma here. Do you remember?
    }
  }
  samples
}

sample_variance <- function(X, z){
  # For sampling Sigma in Gibbs sampler.
  # p(Sigma | others)
  sample_var <- matrix(0, k, d)
  for (k in 1:K){
    indexes <- which(z == k, arr.ind<-TRUE)
    n_k <- sum(z == k)
    if (n_k>1){
      data_k <- X[indexes,]
      sample_mean <- colSums(data_k)/n_k
      alpha.hat <- rep(alpha_0 + n_k/2, d)
      Beta.term1 <- colSums((sweep(data_k,2,sample_mean))^2)
      Beta.term2 <- (n_k/(1 + kappa_0*n_k))*((sample_mean - mu_0)^2)
      Beta.hat <- Beta_0 + .5*(Beta.term1 + Beta.term2)
      for (dim in 1:d){
        sample_var[k,dim] <- 1/rgamma(1,shape=alpha.hat[dim], rate=Beta.hat[dim])
      }
    }else{sample_var[k,] <- 1/rgamma(d,shape=alpha_0, rate=Beta_0)}
  }
  sample_var
}


##############################
############EM################
##############################

# Performs the EM algorithm.

# Takes in all model parameters and outputs
# a results object containing optimized parameters,
# as well as the objective function evaluated
# at convergence.
EM <- function(mu, Sigma, alpha, z){

  mu0 <- mu # for first convergence check
  alpha0 <- alpha # for first convergence check
  Sigma0 <- Sigma # for first convergence check

  for(r in 1:R.EM){
    ###########################
    #E-Step ###################
    ###########################
    num.w <- matrix(0,N,K)
    log_num.w <- matrix(0,N,K)
    for(i in 1:N){
      for(k in 1:K){
        log_num.w[i,k] <- sum(dnorm(X[i,], mu[k,], 
                                    sqrt(Sigma[k,]),
                                    log=TRUE)) + log(alpha[k])
      }}
    log_denom <- apply(log_num.w, 1, log_sum_exp)
    w <- exp(log_num.w - matrix(rep(log_denom, each=K), N, K, byrow=TRUE))
    
    ############################
    #M-Step ####################
    ############################
    
    #1. alpha update
    N.k <- colSums(w, na.rm=TRUE)
    alpha <- N.k/N
    
    #2. mu & sigma update
    for(k in 1:K){
      w.ik <- matrix(w[,k],N,d)
      mu[k,] <- colSums(w.ik*X)/N.k[k]
      mu.k <- matrix(mu[k,], N, d, byrow=TRUE)
      Sigma[k,] <- colSums(w.ik*(X-mu.k)^2)/N.k[k]
      Sigma[k,which(Sigma[k,]==0)] <- 0.1^5
    }
    if (max(abs(mu0-mu))<0.1^5 & 
        max(abs(alpha0-alpha))<0.1^5 &
        max(abs(Sigma0-Sigma))<0.1^5){
      break
    }
    mu0 <- mu
    Sigma0 <- Sigma
    alpha0 <- alpha
    if(r==R.EM){print("Warning: Algorithm does not converge")}
  }
  Obj.after <- Obj(data=X, mu, Sigma, alpha)
  results <- list("means"=mu, "variances"=Sigma,
                 "pi.mix"=alpha, "z"=apply(w,1,which.max),
                 "Obj"=Obj.after)
}


########################################
###############Gibbs####################
########################################

Gibbs <- function(start_obj, mu.init, sig.init, pi.init, z.init, Gibbs_improvement){
  # Gibbs sampler. 
  stop.all <- TRUE # Bool which indicates whether

  # Initialize sample chains
  z_samples <- matrix(0, N, R.Gibbs)
  pi_samples <- matrix(0, K, R.Gibbs)
  mean_samples <- array(rep(0,K*d*R.Gibbs), c(K,d,R.Gibbs))
  variance_samples <- array(rep(0,K*d*R.Gibbs), c(K,d,R.Gibbs))

  # Set first sample to output of EM algorithm
  mean_samples[,,1] <- mu.init
  variance_samples[,,1] <- sig.init
  z_samples[,1] <- z.init
  pi_samples[,1] <- pi.init

  for (g in 2:R.Gibbs){

    if (Gibbs_improvement==TRUE){print("Gibbs has improved our solution")}
    print(sprintf("Gibbs iteration: %d", g))

    z_samples[,g] <- sample_z(X, mean_samples[,,g-1],
                              variance_samples[,,g-1], pi_samples[,g-1])
    pi_samples[,g] <- sample_pi(z_samples[,g])
    mean_samples[,,g] <- sample_means(X, variance_samples[,,g-1],
                                     z_samples[,g])
    variance_samples[,,g] <- sample_variance(X, z_samples[,g])

    Obj.gibbs <- Obj(data=X, mu=mean_samples[,,g],
                    Sigma=variance_samples[,,g],
                    pi.mix=pi_samples[,g])

    last_index <- g

    if (Obj.gibbs > start_obj){
      Gibbs_improvement <- TRUE # Indicates that Gibbs has found a better soln.
      stop.all <- FALSE # If FALSE, outer while loop with run EM again.
      break
    }
  }
  results <- list("z_samples"=z_samples[,last_index],
                 "pi_samples"=pi_samples[,last_index],
                 "mean_samples"=mean_samples[,,last_index],
                 "variance_samples"=variance_samples[,,last_index],
                 "Obj"=Obj.gibbs,
                 "stop"=stop.all,
                 "Gibbs_improvement"=Gibbs_improvement)
}




#########################
###### MAIN LOGIC #######
#########################


#############################
#   Initialize Parameters  #
#############################

z <- sample(1:K,N,replace<-TRUE)
mu <- matrix(0,K,d)
Sigma <- matrix(1,K,d) #set initial values of all sigma to be one
for(k in 1:K){
  mu[k,] <- apply(X[z==k,],2,mean)}
pi.mix <- rep(1/K,K) # equal weights at the beginning

EM_results <- EM(mu, Sigma, pi.mix, z) # First run of EM
Gibbs_improvement <- FALSE # Bool to indicate if Gibbs has improved our solution.

while(TRUE){

  EM_means <- EM_results$means
  EM_variances <- EM_results$variances
  EM_pi.mix <- EM_results$pi.mix
  EM_z <- EM_results$z

  EM_obj <- Obj(X, mu=EM_means, Sigma=EM_variances,
                pi.mix=EM_pi.mix)

  print("EM obj is")
  print(EM_obj)

  # Allow Gibbs to search for improvement
  new.attempt <- Gibbs(start_obj<-EM_obj, mu.init<-EM_means,
                      sig.init<-EM_variances, pi.init<-EM_pi.mix,
                      z.init <- EM_z, Gibbs_improvement<-Gibbs_improvement)
  
  # Bool to track if Gibbs has found an improvement over EM
  Gibbs_improvement <- new.attempt$Gibbs_improvement 

  if (new.attempt$stop==TRUE){
    ######ALGORITHM END#######

    # Met if the Gibbs sampler ran through R iterations
    # without finding a better solution.
    # Prints out final parameters.
    final.means <- EM_means
    final.variances <- EM_variances
    final.pi.mix <- EM_pi.mix
    final.z <- EM_z
    final.obj <- EM_obj
    
    BIC = calc.BIC(final.obj, final.means,
                   final.variances, final.pi.mix)
    print(sprintf("BIC is %f", BIC)) # Return BIC
    
    # Print results
    print(final.means)
    print(final.variances)
    print(final.pi.mix)
    print(final.z)
    print(final.obj)
    
    break
  }
  else{EM_results <- EM(mu=new.attempt$mean_samples,
                         Sigma=new.attempt$variance_samples,
                         alpha=new.attempt$pi_samples,
                         z=new.attempt$z_samples) # Continue with new EM
      }
}

