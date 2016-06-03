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
data$Class[data$Class == 2] <- -1

for (col in names(data)){
    NAs <- is.na(data[,col])
    data[NAs, col] <- mean(data[,col], na.rm = T)
}
D <- as.matrix(data[,1:9])
D <- cbind(D, rep(1, nrow(D)))
colnames(D)[10] <- "Intercept"
O <- as.numeric(data[,"Class"])

################
# Kernel Stuff #
################
ker_identity <- function(u, v) as.numeric(u %*% v)
ker_polyn <- function(u, v) (as.numeric(u %*% v) + 1)^2
ker_gauss <- function(u, v, sd = .3) exp(-as.numeric((u - v) %*% (u - v)) / (2 * sd^2))

gram <- function(data, kernel, ...){
    dots <- apply(data, 1, function(u){
        apply(data, 1, function(v){
            kernel(u, v, ...)
        })
    })
    as.matrix(dots, nrow = nrow(data), ncol = ncol(data))
}

#########################
# SVM training with SGD #
#########################
SGD_SVM <- function(data, target, C = 1, nsteps=50000, step=0.1){
    
    # Intialization
    w <- runif(ncol(data))
    
    for (i in 1:nsteps){
        # Gets a random tuple
        row_i <- sample(x=nrow(data), 1)
        tuple <- data[row_i,]
        obj <- target[row_i] 
        
        # Computes gradient
        grad_fit <- if (obj * as.numeric(w %*% tuple) < 0) -1 * obj * tuple else 0
        grad_reg <- 1 / (2 * C * nsteps) * w
        
        # Updates w
        w <- w - step * (grad_fit + grad_reg)
    }
    
    return(w)
}

# Trains
out <- SGD_SVM(D, O)
coeff_logreg <- out[names(out) != "Intercept"]
offset_logreg <- out[names(out) == "Intercept"]

# Predicts
scores <- as.matrix(data[,-10]) %*% coeff_logreg + rep(offset_logreg, nrow(data))
data$pred <- sapply(scores, function(s) ifelse(s > 0, 1, -1)) 

# Evaluates
misPerc <-  1 - (sum(data$pred == data$Class) / nrow(data))
cat("SVM with Hinge Loss:", misPerc, "\n")

##########################################
# SVM training with Platt's SMO algorithm #
##########################################
SGD_SMO <- function(data, target, C = 1, nsteps=50000, kernel = ker_identity){
    
    a <- runif(nrow(data))
    n_iter <- 0
    
    cat("Prepares the Gram matrix...")
    K <- gram(data, kernel)    
    cat("Done...")
    
    for (n_iter in 1:nsteps){
        
        # Gets a couple of constraints which violate the KKT conditions
        check_KKTs <- sapply(1:nrow(data), function(i){
            # Gets prediction for tuple
            tuple <- data[i,]
            y <- sum(a * target * apply(data, 1, kernel, tuple))
            # Checks inequalities
            if (a[i] == 0){
                a[i] * y >= 1
            } else if (a[i] > 0 & a[i] < C){
                a[i] * y >= .99 & a[i] * y <= 1.01
            } else if (a[i] == C){
                a[i] * y <= 1
            }})
        
        if (all(check_KKTs)){
            cat("Algorithm convereged!\n")
            return(a)
        }
        
        violated <- which(!check_KKTs)
        m <- violated[1]
        n <- violated[2] 
        
        # Reoptimizes for these guys
        # Gets old values
        t_1_old <- data[m,]
        o_1 <- target[m]
        y_1 <- sum(a * target * apply(data, 1, kernel, t_1_old))
        e_1 <- y_1 - o_1 
        
        t_2_old <- data[n,] 
        o_2 <- target[n]
        y_2 <- sum(a * target * apply(data, 1, kernel, t_2_old))
        e_2 <- y_2 - o_2
        
        # Optimizes a_2
        der <- 2 * kernel(t_1_old, t_2_old) -
                   kernel(t_1_old, t_1_old) -
                   kernel(t_2_old, t_2_old)
        a_2_new <- a[n] -  y_2 * (e_2 - e_1) / der
        
        # Clips a_2
        if (o_1 != o_2){
            L <- max(0, a[n] - a[m])
            H <- min(C, C + a[n] - a[m])
        } else {
            L <- max(0, a[n] + a[m] - C)
            H <- min(C, a[n] + a[m])
        }
        a_2_new <- if (a_2_new < L){
            L
        } else if (a_2_new > H){
            H
        } else {
            a_2_new
        }
        
        # Gets a_1 
        a_1_new <- a[m] + o_1 * o_2 * (a[n] - a_2_new)
        
        # Updates
        a[m] <- a_1_new
        a[n] <- a_2_new
    }
        
}
out <- SGD_SMO(D, O)