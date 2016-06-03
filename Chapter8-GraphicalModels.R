library(dplyr)

# Utils
h <- function(l) paste(sort(unique(l)), collapse = "")

######################
# Prepares the graph #
######################
FG <- matrix(
    c(0, 0, 0, 0, 0, 1, 0, 0,
      0, 0, 0, 0, 0, 1, 1, 1,
      0, 0, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 0, 1, 0,
      1, 1, 0, 0, 0, 0, 0, 0,
      0, 1, 0, 0, 1, 0, 0, 0,
      0, 1, 1, 1, 0, 0, 0, 0),
    byrow = T, nrow = 8
)
rownames(FG) <- 1:nrow(FG)
colnames(FG) <- 1:ncol(FG)

V_nodes <- 1:5
F_nodes <- 6:8

potentials <- list() 
potentials[[h(c(1,2))]] <- data.frame(
    x1 = c(0,0,1,1),
    x2 = c(0,1,0,1),
    p  = c(0.5, 1, 0.5, 0)
)

potentials[[h(c(2,3,4))]] <- data.frame(
    x2 = c(1,1,1,1,0,0,0,0),
    x3 = c(0,0,1,1,0,0,1,1),
    x4 = c(0,1,0,1,0,1,0,1),
    p = c(0.05, 0.20, 0.05, 0.08,
          0.05, 0.20, 0.05, 0.32)
)

potentials[[h(c(2,5))]] <- data.frame(
    x2 = c(0, 1, 0, 1),
    x5 = c(0, 0, 1, 1),
    p  = c(0, 0.5, 1, 0.5)
)

##############
# Algorithms #
##############

max_sum <- function(graph, V_nodes, F_nodes, potentials){
    
    # First, backtracks to get order of nodes
    target <- V_nodes[length(V_nodes)]
    to_explore <- c(target)
    
    while(length(to_explore) < nrow(graph)){
        old <- to_explore
        to_explore <- c()
        for (p in old){
            parents <- as.numeric(rownames(graph)[graph[, p] == 1])
            parents <- setdiff(parents, to_explore)
            parents <- setdiff(parents, old)
            to_explore <- c(to_explore, parents, p)
        }
    }
    
    # Reorders matrix
    graph_fwd <- graph
    graph_bwd <- graph
    for (i in 1:length(to_explore)){
        node <- to_explore[i]
        graph_fwd[node, to_explore[1:i]] <- 0
        graph_bwd[node, to_explore[i:length(to_explore)]] <- 0
    }
    
    
    messages <- data.frame(from=c(), to=c(), v=c(), p=c())
    
    # First phase of the message passing algorithm
    for (it in 1:2){
        if (it==1){
            cat("*********GOING FORWARD")
            graph <- graph_fwd 
        } else {
            cat("\n\n********* GOING BACKWARDS")
            graph <- graph_bwd 
            to_explore <- rev(to_explore)
        }
        
        for (n in to_explore){
            
            cat("\n***Exploring node", n, "\n")
            
            in_messages <- filter(messages, to == n)
            dest <- which(graph[n,] == 1)
            origin <- which(graph[,n] == 1)
            if (length(dest) < 1) next
            
            cat("Incoming messages:\n")
            print(in_messages)
            
            # First case: variable -> factor node
            if (n %in% V_nodes){
                
                # If leaf node:
                if (length(origin) < 1){
                    cat("Leaf node detected")
                    messages <- rbind(messages,
                         data.frame(from = c(n, n), to = c(dest, dest),
                                    v=c(0, 1), p=c(1, 1)))
                } else {
                    # Otherwise
                    for (one_dest in dest){
                        new_message <- in_messages %>%
                                        filter(!(from %in% one_dest)) %>%
                                        group_by(v) %>%
                                        summarise(p = prod(p)) %>%
                                        ungroup
                        
                        new_message <- new_message %>%
                                        mutate(from = rep(n, nrow(new_message)),
                                               to   = rep(one_dest, nrow(new_message)))
                                               
                        print(new_message)
                        messages <- rbind(messages,new_message)
                    }
                }   
                
            # Second case: factor -> variable node
            } else if (n %in% F_nodes){
                # Gets potentials
                k <- h(c(in_messages$from, dest))
                tab <- potentials[[k]]
                
                for (one_dest in dest){
                
                    # Does the marginalization
                    for (x in origin){
                        join_col <- paste0("x", x) 
                        
                        tab <- in_messages %>%
                                filter(from == x) %>%
                                select(v, p) %>%
                                rename_(.dots = setNames(list("p", "v"), c("marg_p", join_col))) %>%
                                inner_join(tab, by = c(join_col))
                        
                        tab <- tab %>%
                                mutate(p = p * marg_p) %>%
                                select(-marg_p)
                    }
                    
                    new_message <- tab %>%
                                    group_by_(.dots = paste0("x", one_dest)) %>%
                                    summarize(p = sum(p))
                        
                    nrows <- nrow(new_message)
                    new_message <- new_message %>%
                                    mutate(from = rep(n, nrows),
                                           to   = rep(one_dest, nrows)) %>%
                                    rename_(.dots = setNames(paste0("x", one_dest), "v"))
                    
                    print(new_message)
                    
                    messages <- rbind(messages,new_message)
                }
            }
            
        }
        
    }
    cat("\n\n\nDone!!!!\n\n")
   
    probas <- messages %>%
                group_by(to, v) %>%
                filter(to %in% V_nodes) %>%
                summarize(p = prod(p))
    
    print(probas)
    
}

max_sum(FG, V_nodes, F_nodes, potentials)