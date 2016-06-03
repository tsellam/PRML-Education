library(dplyr)
library(ggplot2)

####################
# Data preparation #
####################
data <- read.csv("~/Projects/Various/Bishop/breast-cancer-wisconsin.data",
                 header=FALSE,na.strings = "?")
data[,1] <- NULL

names(data) <- c("ClumpThickness",
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

pairs(data[,1:9],col=data$Class)

###########################
# Multi-layer perceptron! #
###########################
sigma <- function(x) 1 / (1 + exp(-x))
d_sigma <- function(x) sigma(x) * (1 - sigma(x))

loss <- function(w, D, O){
    a <- D %*% w
    y <- sigma(a)
    err <- apply(cbind(y,O), 1,
        function(r) r[2] * log(r[1]) + (1 - r[2]) * log(1 - r[1])
    ) 
    -err
}
d_loss <- function(w, D, O){
    a <- D %*% w
    Y <- sigma(a)
    spr <- Y - O
    spr %*% D
}

D <- as.matrix(data[,1:9])
D <- cbind(D, rep(1, nrow(D)))
colnames(D)[10] <- "Intercept"
O <- as.numeric(data[,"Class"])

# Test
w <- c(0.531259214,0.006879861,0.330098192,0.239278440,
       0.067566791,0.406755479,0.409317156,0.146323109,
       0.548835040,-9.672751)

multilayer_perceptron <- function(D, O, network, nsteps=50000, step=0.01){
    
    # Inialization
    rnd <- matrix(runif(nrow(network) * ncol(network)) * 2 -1,
                  nrow = nrow(network),
                  ncol = ncol(network))
    weights <- rnd * network
     
    for (i in 1:nsteps){
        
        # Gets a random tuple
        row_i <- sample(x=nrow(D), 1)
        tuple <- D[row_i,]
        obj <- O[row_i] 
        
        # Forward propagation to get predictions
        y <- rep(0, nrow(network))
        for (i in 1:nrow(network)){
            inputs <- which(network[i,] != 0)
            if (length(inputs) == 0){
                y[i] <- tuple[i]
            } else {
                y[i] <- sigma(weights[i,inputs] %*% y[inputs])
            }
        }
        
        # Backward propagation to get gradient
        grad <- matrix(0, nrow = nrow(network), ncol = ncol(network))
        err <- rep(0, nrow(network))
        for (i in nrow(network):1){
            outputs <- which(network[,i] != 0)    
            inputs  <- which(network[i,] != 0)  
            if (length(outputs) == 0){
                err[i]   <- y[i] - obj
                grad[i,inputs] <- y[inputs] * err[i]
            } else{
                err[i] <- d_sigma(weights[i,inputs] %*% y[inputs]) *
                            weights[outputs,i] %*% err[outputs]
                grad[i,inputs] <- y[inputs] * err[i]
            }
        }
        
        # Updates weight vector
        weights <- weights - step * grad
    }
    
    weights
}

network <- matrix(c(
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,
    1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,
    1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,
    1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,
    1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,0),
    nrow=16, byrow = T)

net_predict <- function(tuple, weights){
    tmp <- rep(0, nrow(weights))
    for (i in 1:nrow(weights)){
        inputs <- which(weights[i,] != 0)
        if (length(inputs) == 0){
            tmp[i] <- tuple[i]
        } else {
            tmp[i] <- sigma(weights[i,inputs] %*% tmp[inputs])
        }
    }
    return(tmp[length(tmp)])
}

WW <- multilayer_perceptron(D, O, network)
data$pred <- apply(D, 1, net_predict, WW)
data$pred <- ifelse(data$pred > 0.5, 1, 0)

# Evaluates
misPerc <-  1 - (sum(data$pred == data$Class) / nrow(data))
cat("NN:", misPerc, "\n")

data$pred <- NULL

########################
# Logistic regression! #
########################
log_reg <- function(data, j_class = 10, nsteps = 25){
    
    D <- as.matrix(data[,-j_class])
    D <- cbind(D, rep(1, nrow(D)))
    t <- data[,j_class]
    w <- rep(0, ncol(D)) 
    
    for (i in 1:nsteps){
        # Gets the R and S matrix
        y <- sigma(D %*% w)
        R <- diag(as.numeric(y*(1 - y)))
        S <- t(D) %*% R
        
        # Gets the objective vector
        z = D %*% w - solve(R, y - t)
        
        #Upates
        w <- solve(S %*% D, t(D) %*% R %*% z)
    }
    
    return(w)
}

out <- log_reg(data)
coeff_logreg <- out[1:nrow(out)-1,]
offset_logreg <- as.numeric(t(out[nrow(out),]))

# Predicts
scores <- as.matrix(data[,-10]) %*% coeff_logreg + rep(offset_logreg, nrow(data))
data$pred <- sapply(scores, function(s) ifelse(s > 0, 1, 0)) 

# Evaluates
misPerc <-  1 - (sum(data$pred == data$Class) / nrow(data))
cat("Logistic Regression:", misPerc, "\n")