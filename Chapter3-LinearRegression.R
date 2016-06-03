library(tidyr)
library(dplyr)
library(ggplot2)
library(MASS)
select <- dplyr::select

data <- read.csv("~/Projects/Various/Bishop/abalone.data", header=FALSE)
names(data) <- c("sex",
                 "length",
                 "diameter",
                 "height",
                 "whole_weight",
                 "shucked_weight",
                 "viscera_weight",
                 "shell_weight",
                 "rings")
M <- as.matrix(data.frame(lapply(data, as.numeric)))
M <- scale(M[,-1])

M_train <- M[1:floor(nrow(M) * 0.8),]
T_train <- M_train[, ncol(M), drop=FALSE]
colnames(T_train) <- "truth"
M_train <- M_train[, -ncol(M)]

M_test <- M[ceiling(nrow(M)*0.8):nrow(M),]
T_test <- M_test[, ncol(M), drop=FALSE]
colnames(T_test) <- "truth"
M_test <- M_test[, -ncol(M)]


# Algorithm 1: least squares/MLE
get_weights_mle <- function(M, T, lambda=0){
    
    inv_left_side <- lambda * diag(ncol(M)) + t(M) %*% M
    right_side <- t(M) %*% T
    
    # TO DO: speed up with QR decomposition
    out <- solve(inv_left_side, right_side) 
    
    colnames(out) <- c(paste0("mle", lambda))
    rownames(out) <- colnames(M)    
    out
}
w_mle <- get_weights_mle(M_train, T_train)
w_mle_reg <- get_weights_mle(M_train, T_train, 10)


# Algorithm 2: gradient descent
get_weights_SGD <- function(M, Ts, lambda=0, step=0.005, MAX=50000){
    
    w <- rep(0.1, ncol(M))
    w_values <- matrix(nrow=MAX, ncol=ncol(M))
    
    # Scanning and updating the weight
    for (i in 1:MAX){
       
        # Gets random tuple
        t <- sample(1:nrow(M), 1)
        Mn <- M[t,]
        Tn <- Ts[t]
         
        # Computes gradient
        delta_regul <- (lambda / MAX) * w
        delta_error <- (Tn - w %*% Mn) * Mn
        
        # Updates weight vector
        w <- w + step * (delta_regul + delta_error)
        w_values[i,] <- w
    }
    
    # Pretty output
    out <- as.matrix(w, nrow=1) 
    colnames(out) <- c(paste0("SGD", lambda))
    rownames(out) <- colnames(M)
    
    out
}
w_SGD <- get_weights_SGD(M_train, T_train)


# Algorithm 3: Bayesian Linear Regression
bayes_linreg <- function(M, Ts, MAX= 500){
   
    # Initialization
    alpha <- 0.1
    beta  <- 0.1
    history <- matrix(nrow=MAX, ncol=2)
    
    # Precomputations
    M_gram <- t(M) %*% M
    lambdas <- eigen(M_gram)$values
    
    # Iterations
    for (i in 1:MAX){
        # Gets posterior precision matrix
        p_prec <- alpha * diag(ncol(M)) + beta * M_gram
        # Gets posterior mean
        r_side <- t(M) %*% Ts
        p_mean <- beta * solve(p_prec, r_side)
        
        # Gets degrees of freedom
        L <- sum(sapply(lambdas, function(l) l / (alpha + l)))
        
        # Gets errors
        errors <- (Ts - M %*% p_mean)
        total_errors <- sum(errors^2)
        
        # Updates hyperparameters
        alpha <- drop(L / (t(p_mean) %*% p_mean))
        beta <- drop(( nrow(M) - L ) / total_errors)
        history[i, ] <- c(alpha, beta)
    }
    
    cat("Alpha:", alpha, "\n")
    cat("Beta", beta, "\n")
    
    # Pretty output
    out <- as.matrix(p_mean, nrow=1) 
    colnames(out) <- c(paste0("Bayes"))
    rownames(out) <- colnames(M)
    return(out)
}
w_bayes <- bayes_linreg(M_train, T_train)

#############
# Benchmark #
#############
estimators <- cbind(w_mle, w_mle_reg, w_SGD, w_bayes)

# Plotting regression lines
coef <- as.data.frame(estimators) %>%
        mutate(predictor = rownames(estimators)) %>%
        gather_("method", "coefficient", colnames(estimators))

one_predictions <- data.frame(cbind(M_train, T_train)) %>%
            gather_("predictor", "x", colnames(M_train)) %>%
            inner_join(coef) %>%
            mutate(prediction = coefficient * x) %>%
            select(-coefficient)

graph <- ggplot(one_predictions, aes(x=x, y=prediction,
                             color=method, fill=method)) +
            geom_line() +
            geom_point(aes(y=truth), color="blue", alpha=0.3) +
            facet_wrap(~predictor, scales = "free")
    
print(graph)

# Test error
predictions <- M_test %*% estimators

errors <- as.data.frame(cbind(predictions, T_test)) %>%
            gather_("method", "prediction", colnames(predictions)) %>%
            group_by(method) %>%
            summarize(errors = sum((truth - prediction)^2)/nrow(T_train))

print(errors)