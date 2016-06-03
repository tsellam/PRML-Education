library(mvtnorm)

#########################
# Loads and prints data #
#########################
data<- read.csv("~/Projects/Various/Bishop/forestfires.csv", header=TRUE)
data <- data[,-c(1:4, 12)]
data <- na.omit(data)
data <- scale(data)
pairs(data, pch=".")

#############
# Algorithm #
#############
EMM <- function(data, n, imax = 25){
    
    # INITIALIZATION
    km <- kmeans(data, n)
    centers   <- lapply(1:n, function(i) km$centers[i,])
    variances <- lapply(1:n, function(i){
        items  <- data[km$cluster == i,]
        center <- centers[[i]]
        # Warning! Pivots rows and columns
        diffs <- apply(items, 1, function(x) (x - center))
        diffs %*% t(diffs) / nrow(items)
    })
    coeffs <- as.list(km$size / nrow(data))
    
   # ITERATIONS
    for (it in 1:imax){
        # E-step
        resp <- apply(data, 1, function(row){
            densities <- sapply(1:n, function(i){
                coeffs[[i]] * dmvnorm(row, centers[[i]], variances[[i]])
            })
            densities / sum(densities)
        })
        resp <- t(resp)
         
        # M-step
        sum_contrib <- apply(resp, 2, sum)
        centers <- lapply(1:n, function(i){
            resp[,i] %*% as.matrix(data) / sum_contrib
        })
        variances <- lapply(1:n, function(i){
            center <- centers[[i]]
            diffs <- t(apply(data, 1, function(x) (x - center)))
            coeffs <- resp[,i] / sum_contrib[[i]]
            t(diffs) %*% diag(coeffs) %*% diffs
        })
        coeffs <- sum_contrib / nrow(data)
    }
    
    cat("Done!\n")
    list(centers=centers, variances=variances, coeffs=coeffs)
}

out <- EMM(data, n=2)

#####################################
# Simulates data from distributions #
#####################################
simul_sets <- lapply(1:length(out$centers), function(n){
    simul <- rmvnorm(ceiling(out$coeffs[[n]] * 1000), out$centers[[n]], out$variances[[n]])
    class <- rep(n+1, ceiling(out$coeffs[[n]] * 1000))
    cbind(simul, class)
})
simul_data <-  do.call("rbind", simul_sets)

all_data <- cbind(data, rep(1, nrow(data)))
all_data <- rbind(all_data, simul_data)
all_data <- all_data[nrow(all_data):1,]
              
pairs(all_data[,1:(ncol(all_data)-1)], col=all_data[,ncol(all_data)], pch = ".")
