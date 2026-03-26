# Simulate Hakwes Processes ####

# in the case of Hawkes process, we cannot fix a constant value for rho_bar, since:
# - if is too large the probability can be always low, and so we simulate more and
#   more points to find the actual points
# - if is too small, the probability can be >1 and must be calculated again

# the solution could be to constantly update the value of rho_bar so that is incremented
# each time we accept a new point.

simulate_Hawkes_t <- function(mu_func, mu_bar, T, g_func){
  t=0
  all_t <- c(NULL)
  delta_t_1 <- rexp(n = 1, rate = mu_bar)
  t_1 <- t + delta_t_1
  prob_1 <- mu_func(t_1) / mu_bar
  if (runif(1)<=prob_1){
    all_t <- c(all_t, t_1)
  }
  t <- t_1
  
  lambda_bar <- mu_bar + g_func(0) # since g is decreasing, the max is g(0)
  while (t<=T){
    
    delta_t <- rexp(n = 1, rate = lambda_bar)
    candidate_t <- t + delta_t
    if (candidate_t>T){
      break
    }
    triggering_factor <- sum(g_func(candidate_t - all_t))
    lambda_real <- mu_func(candidate_t) + triggering_factor
    prob <- lambda_real / lambda_bar # add if to check lambda real<lambda_bar

    if (runif(1) <= prob){
      all_t <- c(all_t,candidate_t) # we first add the new t in the vector
      lambda_bar <- mu_bar + triggering_factor + g_func(0) 
      # then the new maximum value will be the max value for mu and g + the 
      # contribute of the triggering factor, that since g is decreasing 
      # we have:
      
      # delta_t high, the influence of the previous t is lower, since we have
      # a bigger value in the exp because t*>>ti and the new element in the sum
      # does not equiparate the decay for the previous ti and probably we 
      # will reduce the upperbound, became closer to just mu_bar, 
      
      # delta_t low, where we still have
      # influence of the previous point, i.e the difference between t* and ti is not
      # that big, so the exp(t*-ti) will not be so low as before and the 
      # new point t* will compensate the decay of t*-ti since we add the new element
      # in the sum t*-t^(*-1) that allows the upperbound to increase 
    } else{
      lambda_bar <- mu_bar + triggering_factor + g_func(0) 
      # we update lambda_bar even if fails because means that we are distant from
      # the upperbound and if we continue we just continue to refuse point 
      # (computationally inefficient) but if we update the lambda_bar we would 
      # decrease the upperbound but still considering the max jump that the next
      # time could have, i.e. g(0)
    }
    t <- candidate_t
  }
  return(sort(all_t))
}




# Simulation ####

set.seed(123)
T_end<- 2000
mu_bar <- 0.5
alpha <- 0.8 # must be lower than 1
beta <- 1.2

# function
mu_f <- function(t) { return(mu_bar) } # background constant
g_f <- function(t) { return(alpha * beta *exp(-beta * t)) } # exponenetial kernel

# simulation
events <- simulate_Hawkes_t(
  mu_func = mu_f, 
  mu_bar = mu_bar, 
  T = T_end, 
  g_func = g_f)

# results
cat("Number events generated:", length(events), "\n")
head(events)
# Plot counting process N(t)
plot(c(0, events, T_end), 0:(length(events)+1), type="s", 
     main="Simulation Hawkes process", 
     xlab="Time t", ylab="N(t)")



# Inference ####


# We define the negative log likelihood

nll_hawkes <- function(par, events, T) {
  
  # par contains mu, alpha and beta
  mu <- par[1]
  alpha <- par[2]
  beta <- par[3]
  
  # Mu and Beta must be > 0. Alpha must be >= 0 and < 1.
  if (mu <= 0 || alpha < 0 || alpha >= 1 || beta <= 0) {
    return(1e+10) # since these are conditions for the process, we return an high value
  }
  
  N <- length(events)
  
  # The first parte of the integral is simply integral(mu) between 0 and T in dt
  # that since mu is constant gives simply mu*T
  integral_term <- mu * T
  
  # Second part of the integral luckily is a closed form solution given by
  # (sum(alpha * (1 - exp(-beta * (T - t_i)))))
  integral_term <- integral_term + sum(alpha * (1 - exp(-beta * (T - events))))

  
  # We evaluate now the first part of the likelihood that is sum(log(lambda(t)))
  sum_log_intensity <- 0
  
  for (i in 1:N) {
    t_i <- events[i]
      
    # We evaluate the triggering factor
    triggering_factor <- 0
    # All the previous time for that i are simply: events[1:(i-1)]
    t_j_vec <- events[1:(i-1)] 
      
    # Kernel: alpha * beta * exp(-beta * (t_i - t_j))
    kernel_values <- alpha * beta * exp(-beta * (t_i - t_j_vec))
    triggering_factor <- sum(kernel_values)
      
    # we add the mu background intensity in t_i
    lambda_t_i <- mu + triggering_factor
    # then we add to the other values in log
    sum_log_intensity <- sum_log_intensity + log(lambda_t_i)
  }
  
  
  # Since here we want to find the MLE but we want to minimize the elements we use
  # the negative log likelihood
  Neg_LL <- integral_term - sum_log_intensity
  
  return(Neg_LL)
}


initial_params <- c(mu = 2.5, alpha = 0.1, beta = 1) 
# We use optim as optimization function
mle_results <- optim(
  par = initial_params, 
  fn = nll_hawkes, 
  events = events, 
  T = T_end, 
  method = "BFGS",
  control = list(fnscale = 1)
)


cat("\nInitial parameters:\n")
print(initial_params)
cat("Real parameters: mu =", mu_bar, ", alpha =", alpha, ", beta =", beta, "\n")
cat("\nStimed parameters:\n")
print(mle_results$par)
cat("\nNeg_LL Minima, i.e MLE:", mle_results$value, "\n")
