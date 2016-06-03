library(mvtnorm)

####################
# Data preparation #
####################
data <- read.csv("~/Projects/Various/Bishop/breast-cancer-wisconsin.data",
                 header=FALSE,na.strings = "?")
data[,1] <- rep(1, nrow(data))

data <- data[100,]

names(data) <- c("Bias",
                 "ClumpThickness",
                 "UniformityCellSize",
                 "UniformityCellShape",
                 "MarginalAdhesion",
                 "SingleEpithelialCellSize",
                 "BareNuclei",
                 "BlandChromatin",              
                 "NormalNucleoli",              
                 "Mitoses",                     
                 "Class")

data$Class[data$Class == 4] <- 1
data$Class[data$Class == 2] <- 0

for (col in names(data)){
    NAs <- is.na(data[,col])
    data[NAs, col] <- mean(data[,col], na.rm = T)
}

pairs(data[,1:9],col=data$Class+1)

X <- as.matrix(data[,-ncol(data)])
C <- data[,ncol(data)]

############################
# MCMC Logistic Regression #
############################

# Initialization
samples <- list()
s <- list()
s$sigma <- 1
s$w <- rep(0, 10)

for (i in 1:5000){
   t <- list() 
   
   # Samples a sigma
   # Proposal distribution: gamma 1 - 1
   x_sigmas <- rgamma(100, 1, 1)
   p_sigma <- pgamma(x_sigmas, 0.01, 0.01) 
   p_w_g_sigma <- sapply(x_sigmas, function(sig) prod(pnorm(s$w, 0, sig)))
   thresholds <- p_sigma * p_w_g_sigma
   
   p_x_proposal <- pgamma(x_sigmas, 1, 1)
   coeff <- 1.5 * max(thresholds) / max(p_x_proposal)
   u <- runif(100, 0, coeff * p_x_proposal)
   x <- x_sigmas[u < thresholds]
   if (length(x) < 1) {
       cat("Sampling failed for Sigma!\n")
       next()
   } else {
       t$sigma <- x[[1]]
   }
   
   # Samples a w 
   # Propsal distribution: N(0, t$sigma)
   x_w <- matrix(rnorm(10 * 1000, 0, t$sigma), ncol = 10)
   
   p_w <- apply(log(pnorm(x_w, 0, t$sigma)), 1, sum)
   p_w_sig <- apply(x_w, 1, function(w) sum(log(pnorm(w, 0, t$sigma))))
   p_t_w <- apply(x_w, 1, function(w){
       proj <- as.numeric(X %*% w)
       ps <- 1 / (1 + exp(-proj))
       likelihoods <- log(ps^C * (1 - ps)^(1 - C))
       sum(likelihoods)
   })
   thresholds <- exp(p_w + p_w_sig + p_t_w)
   
   p_w_proposal <- apply(pnorm(x_w, 0, t$sigma), 1, prod)
   log_coeff <-  max(thresholds) / max(p_w_proposal)
   u <- runif(1000, 0, log_coeff * p_w_proposal)
   x <- x_w[u < thresholds,]
   if (length(x) < 1) {
       stop("Sampling failed for W!\n")
   } else {
       t$w <- x[1,]
   }
   
   samples[[length(samples) + 1]] <- t
   s <- t
   
}

cat(length(samples), " samples generated\n")