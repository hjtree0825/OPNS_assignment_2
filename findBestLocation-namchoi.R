find_best_location <- function(num.sites, assemblers, competetors, union.rate.fn, beta, num.tries){
  
  # Obtain the location parameter
  C_bar_tri <- function(assemblers, suppliers){
    euc_dist <- sum(sqrt((assemblers[1] - suppliers[1])^2 + (assemblers[2] - suppliers[2])^2))
    union_rate <- union.rate.fn(suppliers[1], suppliers[2])
    c_bar_tri <- beta[[1]] * euc_dist + beta[[2]] * union_rate
    
    return(c_bar_tri)
  }
  
  # Find LSE based on the closed form we derived from Part I
  Exp_N <- function(assemblers, suppliers){
    lse <- t(data.frame("lse" = 1:nrow(suppliers)))
    for (i in 1:nrow(suppliers)){
      lse[i] <- -1*C_bar_tri(assemblers, suppliers[i,])
    }
    LSE <- -1*log(sum(exp(unlist(lse))))
    
    return(LSE)
  }
  
  Minimum <- function(suppliers){
    df_sup_input <- matrix(suppliers, nc = 2)
    df_sup_total <- rbind(df_sup_input, as.matrix(competetors))
    min <- t(data.frame("min" = 1:nrow(assemblers)))
    for (a in 1:nrow(assemblers)){
      min[a] <- Exp_N(assemblers[a,], df_sup_total)
    }
    min_sum <- sum(min)
    
    return(min_sum)
  }
  
  opt_sol <- seq(num.tries) %>%
    llply(function(l){
      optim(
        runif(2 * num.sites),
        Minimum,
        control=list(maxit=10000)
      )
    },
    .parallel = TRUE
    )
  
  '%!in%' <- function(x,y) !( '%in%'(x,y) )
  coordinates_cost <- unlist(opt_sol)[names(unlist(opt_sol)) %!in%
                                        c("counts.function", "counts.gradient", "convergence")] %>%
    matrix(ncol= 2*num.sites+1, byrow=T) %>% data.frame()
  colnames(coordinates_cost) <- c(1:(2*num.sites), "cost")
  opt_location <- coordinates_cost %>% filter(cost==min(cost)) %>% select(-cost) %>% matrix(ncol=2)
  
  return(opt_location)
  
}