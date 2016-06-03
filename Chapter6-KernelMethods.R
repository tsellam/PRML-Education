library(dplyr)
library(mvtnorm)

data <- read.csv("~/Projects/Various/Bishop/abalone.data", header=FALSE)
names(data) <- c("sex",
                 "length",
                 "diameter",
                 "height",
                 "whole_weight",
                 "shucked_weight",
                 "viscera_weight",
                 "shell_weight",
                 "target")
data <- sample_n(data, nrow(data), replace = FALSE)
data <- mutate(data, sex = as.numeric(sex))
data <- as.data.frame(scale(data))
train <- data[1:3000,]
test  <- data[3001:nrow(data),]

M_train <- as.matrix(select(train, -target))
M_test  <- as.matrix(select(test, -target))

###########################
# Plain linear regression #
###########################
model_LR <- lm(target ~., data=train)

pred_LR  <- predict(model_LR, select(test, -target))
err_LR <- sum(abs(pred_LR - test$target)) / nrow(test)
cat("Plain linear regression, test set:", err_LR, "\n")




###################
# Kenel functions #
###################
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

##########################
# Naradaya-Watson Kernel #
##########################
dens_smooth <- function(x, data, sd=1){
    diffs <- apply(data, 1, dmvnorm, mean = x, sigma = sd * diag(length(x)))
    (1 / nrow(data)) * sum(diffs)
}

predict_NWK <- function(x, data, target, sd=1){
   local_kernels <- apply(data, 1, ker_gauss, x, sd = 0.3) 
   local_kernels %*% target / sum(local_kernels)
}

pred <- apply(M_test, 1, predict_NWK, M_train, train$target)
err <- sum(abs(pred - test$target)) / nrow(test)
cat("Kernel regression:", err, "\n")

####################
# Gaussian Process #
####################
GP_predict <- function(data, target, kernel=ker_identity){
    cat("Computing Gram Matrix...")
    K <- gram(data, kernel)    
    cat("Done...")
    
    cat("Solving dual equations...")
    a <- solve(K + lambda * diag(nrow(K)), target)
    cat("Done...\n")
    
}

##########################
# Dual Linear Regression #
##########################
dual_lm <- function(data, target, lambda = 0.001, kernel=ker_identity){
    cat("Computing Gram Matrix...")
    K <- gram(data, kernel)    
    cat("Done...")
    
    cat("Solving dual equations...")
    a <- solve(K + lambda * diag(nrow(K)), target)
    cat("Done...\n")
    
    return(a)
}

pred_dual_lm <- function(x, a, data, kernel){
    k_x <- apply(data, 1, kernel, x)
    as.numeric(a %*% k_x)
}

a <- dual_lm(M_train, train$target, , kernel = ker_identity)
pred_id <- apply(M_test, 1, pred_dual_lm, a, M_train, ker_identity)
err_LR <- sum(abs(pred_id - test$target)) / nrow(test)
cat("Dual linear regression, test set:", err_LR, "\n")

a <- dual_lm(M_train, train$target, kernel = ker_polyn)
pred_id <- apply(M_test, 1, pred_dual_lm, a, M_train, ker_polyn)
err_LR <- sum(abs(pred_id - test$target)) / nrow(test)
cat("Dual linear regression, polynomial kernel:", err_LR, "\n")

a <- dual_lm(M_train, train$target, kernel = ker_gauss)
pred_id <- apply(M_test, 1, pred_dual_lm, a, M_train, ker_gauss)
err_LR <- sum(abs(pred_id - test$target)) / nrow(test)
cat("Dual linear regression, gaussian kernel:", err_LR, "\n")

