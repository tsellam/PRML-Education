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

###########################################
# Discriminant model 1: linear regression #
###########################################
# Recoding of the data
data$class1 <- data$Class
data$class0 <- 1 - data$Class

# Classification
prob1 <- lm(class1~., data[,-c(10,12)])$fitted.values
prob0 <- lm(class0~., data[,-c(10,11)])$fitted.values
data$Predict <- apply(cbind(prob0, prob1), 1, which.max) - 1

# Scoring
misLineReg <-  1 - (sum(data$Predict == data$Class) / nrow(data))
cat("Linear Regression:", misLineReg, "\n")

# Plotting errors:
data$class1 - prob1 -> errors1
hist(as.numeric(errors1), breaks = 100)

# Cleanup
data <- select(data,-class0, -class1, -Predict)

#############################
# Discriminant model 2: LDA #
#############################
lda <- function(data, j_class){
    # Computes the two centroids and the variance matrix
    m0 <- colMeans(data[data$Class == 0, -j_class], na.rm = T)
    n0 <- nrow(data[data$Class == 0, -j_class])
    m1 <- colMeans(data[data$Class == 1, -j_class], na.rm = T)
    n1 <- nrow(data[data$Class == 1, -j_class])
    S <- var(data[data$Class == 0, -j_class]) * n0 / (n0 + n1) +
         var(data[data$Class == 1, -j_class]) * n1 / (n0 + n1)
    
    # Returns the coefficients 
    solve(S, m1 - m0)
}

# Prediction
coeff <- lda(data, 10)
threshold <- coeff %*% colMeans(data[,-10], na.rm = T)
data$projection <- apply(data[,-10], 1, function(row) row %*% coeff)
data$prediction <- ifelse(data$projection > rep(threshold, length(data$projection)), 1, 0)

# Scoring
misLDA <-  1 - (sum(data$prediction == data$Class) / nrow(data))
cat("Linear Discriminant Analysis:", misLDA, "\n")

# Plotting
plot <- ggplot(data, aes(x = projection, fill=factor(Class))) + 
        geom_bar() + 
        geom_vline(xintercept=threshold)
print(plot)

# Cleaning
data <- select(data, -projection, -prediction)

#####################################
# Discriminant model 3 : Perceptron #
#####################################
perceptron <- function(data, steps=1000){
    w <- rep(1, ncol(data) - 1)
    for (n in 1:steps){
        tuple <- sample_n(data, 1)
        x <- as.numeric(select(tuple, -PerClass))
        t <- as.numeric(select(tuple, PerClass))
        predict <- ifelse(x %*% w > 0, 1, -1)
        if (predict * t < 0)
            w <- w + x * t
    }
    return(w)
}
data$Intercept <- rep(1, nrow(data))
data$PerClass <- ifelse(data$Class == 1, 1, -1)
coeff_per <- perceptron(select(data, -Class))

data$Predict <- as.matrix(select(data, -Class, -PerClass)) %*% coeff_per
data$Predict <- ifelse(data$Predict > 0, 1, 0)

# Scoring
misPerc <-  1 - (sum(data$Predict == data$Class) / nrow(data))
cat("Perceptron:", misPerc, "\n")

# Cleaning
data <- select(data, -Intercept, -PerClass, -Predict)

####################
# Generative Model #
####################
sigma <- function(a) 1 / (1 + exp(-a))

generativeFit <- function(data, j_class = 10){
    m0 <- colMeans(data[data$Class == 0, -j_class], na.rm = T)
    n0 <- nrow(data[data$Class == 0, -j_class])
    m1 <- colMeans(data[data$Class == 1, -j_class], na.rm = T)
    n1 <- nrow(data[data$Class == 1, -j_class])
    S <- var(data[data$Class == 0, -j_class]) * n0 / (n0 + n1) +
         var(data[data$Class == 1, -j_class]) * n1 / (n0 + n1)
    list(m0 = m0, m1 = m1, n0 = n0, n1 = n1, S = S)
}

model <- generativeFit(data)

# Checks similarity with LDA
P <- solve(model$S)
coeff_gen <- P %*% (model$m1 - model$m0)
offset_gen <- -(1/2) * t(model$m1) %*% (P %*% model$m1) +
               (1/2) * t(model$m0) %*% (P %*% model$m0) +
               log(model$n1 / model$n0)

# Predicts
scores <- as.matrix(data[,-10]) %*% coeff_gen + rep(offset_gen, nrow(data))
data$pred <- sapply(scores, function(s) ifelse(s > 0, 1, 0)) 
# Alternative predictions
rlda <- MASS::lda(Class~., data)
data$altpred<- predict(rlda, data[,-10])$class

# Evaluates
misPerc <-  1 - (sum(data$pred == data$Class) / nrow(data))
cat("Generative LDA:", misPerc, "\n")
# Alternative
misPerc <-  1 - (sum(data$altpred == data$Class) / nrow(data))
cat("Generative LDA - alt:", misPerc, "\n")

# Cleaning up
data <- select(data, -pred, -altpred)

########################
# Logistic regression! #
########################
log_reg <- function(data, j_class = 10, nsteps = 100){
    
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