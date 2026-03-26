# Simulate REM ####
require(mgcv)

# Set the number of individuals in the system
p <- 50
set.seed(1234)

# Set coefficients according to the covariate type
b1 <- 0.5
b2 <- 0.5
b3 <- 0 # Basically I remove reciprocity. So the covariate are estimated but do not influence
        # the lp value

# Lambda0 has an impact when dealing with the time-based definition
lambda0 <- 0.005 #use 0.0005 to see some reliable result in dat matrix

##############################################################################
# STEP 1: Simulate the data
##############################################################################

# --------------------------------------------
# Build initial risk set (all possible directed pairs)
# --------------------------------------------

# Create all possible pairs
riskset <- cbind(rep(1:p, each = p), rep(1:p, p))

# Remove self-pairs
riskset <- riskset[riskset[, 1] != riskset[, 2], ]

# Lenght of riskset
l <- nrow(riskset)

# Parameter for the triggering kernel
alpha <- 0.5
beta <- 0.1 # beta is chosen low to actually have the self-excitation effect that last longer
mu <- 1

# These values are needed for the final plot
target_s <- 14
target_r <- 2
idx_target <- which(riskset[,1] == target_s & riskset[,2] == target_r)
track_data <- matrix(NA, nrow = 5000, ncol = 3)
colnames(track_data) <- c("tm", "lambda_14_2", "sum_hazard")


# --------------------------------------------
# Initialize covariates for the hazard function (Logic: Hawkes)
# --------------------------------------------

# Sender activity time: history of all event times with sender involved
s_act_hist <- rep(list(-Inf), p)

# Receiver popularity time: history of all event times with receiver involved
r_pop_hist <- rep(list(-Inf), p)

# Reciprocity time: history of all times in which a reciprocal event occurred
rec_hist <- rep(list(-Inf), l)

# Initialize data for simulated events
dat <- NULL

# Initialize clock to 0
tm <- 0

# Decide the last event
event <- 0
total <- 5000

# All the covariates at the start will be 0
# I analyze each list that contain p list (one for each node) of all the previous activity as a 
# sender, receiver and reciprocity
s_cov_init <- numeric(l)
r_cov_init <- numeric(l)
rec_cov_init <- numeric(l)
for (k in 1:l){
  
  # Extract sender e receiver k
  s_id <- riskset[k, 1]
  r_id <- riskset[k, 2]
  
  # The history for each node is evaluated below
  s_id_hist <- s_act_hist[[s_id]]
  r_id_hist <- r_pop_hist[[r_id]]
  
  # We then create our vector of covariate for all the possible dyads (obviously
  # if we have a value for sender 1 will be the same for all the possible receiver
  # so that we will have for p=5 the first 4 value of the vector equal)
  s_cov_init[k] <- mu + alpha*beta*sum(exp(-beta*(tm - s_id_hist)))
  r_cov_init[k] <- mu + alpha*beta*sum(exp(-beta*(tm - r_id_hist)))
  
  # For the reciprocity we directly select the dyad k since are obviously all different
  rec_id_hist <- rec_hist[[k]]
  rec_cov_init[k] <- mu + alpha*beta*sum(exp(-beta*(tm - rec_id_hist)))
}


# Compute initial linear predictor with covariates
lp <- (b1* s_cov_init + b2*r_cov_init + b3*rec_cov_init)
hazard <- exp(lp)*lambda0

# --------------------------------------------
# Keep going until we reach number of events = total
# --------------------------------------------

while (event < total) {
  
  # Simulate time until next event 
  # According to exponential (sum of hazard since is the total hazard for all dyads)
  tm <- tm + rexp(1, sum(hazard))
  
  # Small digression to update the values for the final plot
  track_data[event + 1, 1] <- tm                  # Time
  track_data[event + 1, 2] <- hazard[idx_target]  # Specific lambda
  track_data[event + 1, 3] <- sum(hazard)         # Sum(lambda)
  
  # Select which dyad (sender-receiver pair) experiences the event
  # according to Multinomial (hazard / sum(hazard)) since want to identify which
  # actual dyads generated that event
  event.id <- sample(length(hazard), 1, prob = hazard / sum(hazard))
  
  # Identify sender (s) and receiver (r)
  s <- riskset[event.id, 1]
  r <- riskset[event.id, 2]
  
  # Identify reciprocal pair index (receiver sending to sender)
  event.id.rec <- which(riskset[, 1] == r & riskset[, 2] == s)
  
  # Evaluate the covariates for the specific sender s and reciver r and event.id for that s and r
  s_cov_curr <- mu + alpha*beta*sum(exp(-beta*(tm - s_act_hist[[s]])))
  r_cov_curr <- mu + alpha*beta*sum(exp(-beta*(tm - r_pop_hist[[r]])))
  rec_cov_curr <- mu + alpha*beta*sum(exp(-beta*(tm - rec_hist[[event.id]]))) 
  
  # Store relevant information about the event
  dat <- rbind(dat, c(
    tm,
    s, r,
    s_cov_curr,
    r_cov_curr,
    rec_cov_curr
  ))
  
  # Update covariates based on the event
  
  s_act_hist[[s]] <- c(s_act_hist[[s]], tm)
  
  r_pop_hist[[r]] <- c(r_pop_hist[[r]], tm)
  
  # I evaluate the value of the cov for the event for example 1->3 but I insert 
  # the time of occurrence inside 3->1. In this way I each time that happens 
  # an event I actaully check when happened the reciprocal
  rec_hist[[event.id.rec]] <- c(rec_hist[[event.id.rec]], tm)
  
  # Recalculate the hazard function with updated covariates
  s_cov_next <- numeric(l)
  r_cov_next <- numeric(l)
  rec_cov_next <- numeric(l)
  for (k in 1:l){
    
    # Extract sender e receiver k
    s_id <- riskset[k, 1]
    r_id <- riskset[k, 2]
    
    # The history for each node is evaluated below
    s_id_hist <- s_act_hist[[s_id]]
    r_id_hist <- r_pop_hist[[r_id]]
    
    
    # We then create our vector of covariate for all the possible dyads (obviously
    # if we have a value for sender 1 will be the same for all the possible receiver
    # so that we will have for p=5 the first 4 value of the vector equal)
    s_cov_next[k] <- mu + alpha*beta*sum(exp(-beta*(tm - s_id_hist)))
    r_cov_next[k] <- mu + alpha*beta*sum(exp(-beta*(tm - r_id_hist)))
    
    # For the reciprocity we directly select the dyads k since are obviously all different
    rec_id_hist <- rec_hist[[k]]
    rec_cov_next[k] <- mu + alpha*beta*sum(exp(-beta*(tm - rec_id_hist)))
  }
  
  
  lp <- (b1*s_cov_next + b2*r_cov_next + b3*rec_cov_next)
  
  
  hazard <- exp(lp)*lambda0
  event <- event + 1
}

# --------------------------------------------
# Finalize: name the columns of the event log
# --------------------------------------------

colnames(dat) <- c("tm", 
                   "solicitor", 
                   "receiver", 
                   "s_act_cov", 
                   "r_pop_cov",
                   "rec_cov"
)
n.event <- nrow(dat)
dat
tail(dat)
# Compare how the events are spread (for sender)
division<- table(dat[,"solicitor"])
division


##############################################################################
# STEP 2: APPLY HAWKES ON SINGLE SENDER
##############################################################################

# Identify a sender
sender <- 14
# In the actor 14 you clearly see the effect: self excitation and beta small
# In the actor 4 you see that has been estimated a value for beta high: short memory
# the reason of the difference is the stochasticity of the process
# the sender 30 has acted more time as sender-> increasing intensity

# We only extract the time where we have the sender selected
events <- sort(dat[dat[,"solicitor"] == sender, "tm"])
sender_counts <- length(events)
T_end <- max(dat[,"tm"])

cat("We only look for sender=", sender, "\n")
cat("Total events in which acted as a sender:", sender_counts, "\n")

# We define the (negative) log-likelihood for a Hawkes process
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
  
  # The first part of the integral is simply integral(mu) between 0 and T in dt
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


initial_params <- c(mu = 2, alpha = 0.1, beta = 1) 
# We use optim as optimization function
opt_res <- optim(
  par = initial_params, 
  fn = nll_hawkes, 
  events = events, 
  T = T_end, 
  method = "BFGS",
  control = list(fnscale = 1)
)

cat("\n--- Results estimate Hawkes ---\n")
cat("Mu (background rate):", opt_res$par[1], "\n")
cat("Alpha:", opt_res$par[2], "\n")
cat("Beta (decay):", opt_res$par[3], "\n")
cat("----------------------------------\n")


# Plot the intensity
# I create a vector of time to create a smoother plot of the density
t_grid <- seq(0, T_end, length.out = 20000)
lambda_grid <- numeric(length(t_grid))

est_mu <- opt_res$par[1]
est_alpha <- opt_res$par[2]
est_beta <- opt_res$par[3]

for (j in 1:length(t_grid)) {
  t_val <- t_grid[j]
  # I take all the events happened before
  past_events <- events[events < t_val]
  
  if (length(past_events) > 0) {
    decay <- sum(est_alpha *est_beta* exp(-est_beta * (t_val - past_events)))
    lambda_grid[j] <- est_mu + decay
  } else {
    lambda_grid[j] <- est_mu
  }
}

plot(t_grid, lambda_grid, type = "l", col = "blue", lwd = 2,
     main = paste("Hawkes intensity for actor:", sender),
     ylab = "Intensity", xlab = "Time")

points(events, rep(min(lambda_grid), length(events)), col = "red", pch = 16)
legend("topright", legend=c("Estimated intensity", "Real events"),
       col=c("blue", "red"), lty=c(1, NA), pch=c(NA, 16))

# Even if the process seems to have the almost everywhere a small jump, relatively
# to the background, the jump is big

# Obviously we cannot have the same parameters of the Hawkes process covariate, due
# to the fact that the events are not generated relatively to a single actor, but
# are endogenous and depends on all the network. Indeed, we can still see the 
# behaviour of the events that more or less resemble the idea behind the 
# Hawkes process

# Remind that is just a way to give the idea of self-excitation, does not make sense to say
# "this sender have small memory" but "the MLE estimate given the events of this
# sender estimated value that remind an Hawkes process with small memory"




##############################################################################
# STEP 3: SOME PLOTS
##############################################################################


# --------------------------------------------
# Plot for bursts
# --------------------------------------------


library(ggplot2)
library(dplyr)

# We select first of all the first three an last three in the total number of events
activity_ranking <- as.data.frame(dat) %>%
  count(solicitor) %>%
  arrange(desc(n))
top_3_ids <- head(activity_ranking$solicitor, 3)
bottom_3_ids <- tail(activity_ranking$solicitor, 3)
target_ids <- c(top_3_ids, bottom_3_ids)

# We build the dataframe with time gap and a way to show the difference in time gap
dat_raster <- as.data.frame(dat) %>%
  filter(solicitor %in% target_ids) %>% 
  group_by(solicitor) %>%
  arrange(tm) %>%
  mutate(
    # The difference wrt the last event
    time_gap = tm - lag(tm, default = 0),
    # Basically we concentrate the difference in the interval [0,5]
    gap_visual = ifelse(time_gap > 5, 5, time_gap) 
  ) %>%
  ungroup()

# We write this to make sure that in the plot are ordered with the top first
dat_raster$solicitor <- factor(dat_raster$solicitor, levels = rev(target_ids))

# Build the plot
p1 <- ggplot(dat_raster, aes(x = tm, y = solicitor, color = gap_visual)) +
  # shape = 124 is '|'
  geom_point(shape = 124, size = 4) + 
  
  # Color scale: red (gap=0, self-exciting) -> Blue (high gap, decay)
  scale_color_gradientn(
    colors = c("#e74c3c", "#f1c40f", "#3498db"),
    values = c(0, 0.2, 1), # so that we have red for gap<0.2
    name = "Waiting time\n(red = burst)"
  ) +
  
  # We zoom only on the first 25 seconds
  coord_cartesian(xlim = c(0, 90)) + 
  
  theme_minimal(base_size = 14) +
  labs(
    title = "Seeing the Hawkes bursts",
    subtitle = "Red zone indicate self-excitation",
    x = "Time",
    y = "Sender"
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "gray90"),
    legend.position = "bottom"
)
print(p1)



# --------------------------------------------
# Plot for covariate change
# --------------------------------------------


# Select the node that you want to analize
focus_node <- 14 
subset_dat <- as.data.frame(dat)
subset_dat <- subset_dat[subset_dat$solicitor == focus_node, ]

# To visualize the continuous decay, we create an artificial grid
time_grid <- seq(0, max(subset_dat$tm), length.out = 1000)

# Built the s_act_cov to the node to be continuous
s_cov_trace <- sapply(time_grid, function(t_curr) {
  past_events <- subset_dat$tm[subset_dat$tm < t_curr]
  if(length(past_events) == 0) return(mu)
  mu + alpha * beta * sum(exp(-beta * (t_curr - past_events)))
})

df_trace <- data.frame(Time = time_grid, Intensity_Component = s_cov_trace)

p2 <- ggplot() +
  # Continuous line for decay 
  geom_line(data = df_trace, aes(x = Time, y = Intensity_Component), color = "black") +
  # Points are the events
  geom_point(data = subset_dat, aes(x = tm, y = s_act_cov), shape = 1, size = 3) +
  geom_segment(data = subset_dat, aes(x = tm, xend = tm, y = 0, yend = s_act_cov), 
               linetype = "dashed", color = "gray") +
  theme_classic() +
  labs(title = " Covariate dynamics",
       subtitle = paste("Self-excitation dynamics for node", focus_node),
       x = "Time", y = "Sender activity covariate")

print(p2)



# --------------------------------------------
# Plot to show the cumulative number of events
# --------------------------------------------

# Cumulative count for each node
dat_step <- as.data.frame(dat) %>%
  group_by(solicitor) %>%
  arrange(tm) %>%
  mutate(cumulative_events = row_number()) %>%
  ungroup()

# Show only the most and least active
selected_senders <- c(30, 47) # 30 is the one with most events, 47 is the one with less events
dat_step_filtered <- dat_step %>% filter(solicitor %in% selected_senders)

p3 <- ggplot(dat_step_filtered, aes(x = tm, y = cumulative_events, color = factor(solicitor))) +
  geom_step(size = 1) + # geom step remain the same if there are no events and then jump of 1
  theme_bw() +
  labs(
    title = "Number of events as sender",
    subtitle = "How slope changes due to Hawkes effect",
    x = "Time",
    y = "Cumulative number of events",
    color = "Sender"
  ) +
  theme(legend.position = "bottom")
print(p3)



# --------------------------------------------
# Plot for difference with only count (zoomed)
# --------------------------------------------


# Select top sender
sender_stats <- as.data.frame(dat) %>% count(solicitor) %>% arrange(desc(n))
focus_node <- sender_stats$solicitor[1] # select top sender (actor 30)
subset_dat <- as.data.frame(dat) %>% filter(solicitor == focus_node)
event_times <- sort(subset_dat$tm)

# We must use alpha * beta as jump size since in our simulation, without 
# decay effect we do not have +1 but +a*b
jump_size <- alpha * beta 

# Create the grid
T_end <- max(dat[,"tm"])
time_grid <- seq(0, T_end, length.out = 3000)

hawkes_curve <- numeric(length(time_grid))
no_decay_curve <- numeric(length(time_grid))

for(j in 1:length(time_grid)) {
  t_curr <- time_grid[j]
  past_events <- event_times[event_times < t_curr]
  
  if(length(past_events) > 0) {
    # Hawkes formula
    hawkes_curve[j] <- mu + sum(jump_size * exp(-beta * (t_curr - past_events)))
    # No decay formula
    no_decay_curve[j] <- mu + (length(past_events) * jump_size)
  } else {
    hawkes_curve[j] <- mu
    no_decay_curve[j] <- mu
  }
}

df_global_check <- data.frame(
  Time = time_grid,
  Hawkes_Value = hawkes_curve,
  No_Decay_Value = no_decay_curve
)

# To make the plot more understandable we must zoom it, since the global plot
# have as a problem the fact that at the end the two lines are two distant
max_time_zoom <- 20

p4 <- ggplot() +
  # Red line: no decay
  geom_line(data = df_global_check %>% filter(Time < max_time_zoom), 
            aes(x = Time, y = No_Decay_Value, color = "No decay"), 
            size = 0.8, linetype = "dashed") +
  # Blue line: decay
  geom_line(data = df_global_check %>% filter(Time < max_time_zoom), 
            aes(x = Time, y = Hawkes_Value, color = "Hawkes covariate value"), 
            size = 1) +
  # Events
  geom_point(data = subset_dat %>% filter(tm < max_time_zoom), 
             aes(x = tm, y = s_act_cov), 
             shape = 1, size = 2, color = "black") +
  
  scale_color_manual(values = c("No decay" = "#e74c3c", 
                                "Hawkes covariate value" = "#2980b9")) +
  theme_minimal(base_size = 14) +
  labs(
    title = paste("Effect decay: sender", focus_node),
    subtitle = paste("Zoom on first", max_time_zoom, "seconds"),
    x = "Time",
    y = "Value covariate"
  ) +
  theme(legend.position = "bottom")
print(p4)


# --------------------------------------------
# Plot for the "dilution" effect
# --------------------------------------------

# We convert the table saved before in dataframe. The objective of the table saved
# was to show the difference between the value of the covariate and the hazard function
# for a specific dyad: while the covariate remains the same each time the actor selected
# (in this case 14) act as a sender, the value of the hazard function obviously is different
# for each dyad
df_track <- as.data.frame(track_data)

# Evaluate the probability in each time in which any events occurred
df_track$true_prob <- df_track$lambda_14_2 / df_track$sum_hazard

# 2. Evaluate the covariate for each instant
events_sender_14 <- dat[dat[,"solicitor"] == 14, "tm"]
df_track$covariate_value <- sapply(df_track$tm, function(t_curr) {
  past_events <- events_sender_14[events_sender_14 < t_curr]
  if(length(past_events) == 0) return(mu)
  mu + alpha * beta * sum(exp(-beta * (t_curr - past_events)))
})

# Use a scale factor to align the two curves visually
scale_factor <- max(df_track$covariate_value) / max(df_track$true_prob)

p5 <- ggplot(df_track, aes(x = tm)) +
  
  # Blue line: Covariate values for any dyads in which actor 14 is sender
  geom_line(aes(y = covariate_value, color = "Hawkes covariate (sender 14)"), size = 1) +
  
  # Red line: True probability (scaled)
  geom_line(aes(y = true_prob * scale_factor, color = "True probability (14->2)"), 
            size = 0.8, linetype = "dashed") +
  scale_y_continuous(
    name = "Covariate value",
    sec.axis = sec_axis(~ . / scale_factor, name = "True event probability")
  ) +
  
  scale_color_manual(values = c("Hawkes covariate (sender 14)" = "#2980b9", 
                                "True probability (14->2)" = "#e74c3c")) +
  
  theme_minimal(base_size = 14) +
  labs(
    title = "The dilution effect",
    subtitle = "Comparing covariate value vs true probability for dyad 14->2",
    x = "Time",
    color = "Metric"
  ) +
  theme(legend.position = "bottom")

print(p5)

# This approach is a way to inject the Hawkes process in the REM, but since the
# intensity of the triggering kernel in the Hawkes process is directly added to the
# intensity of the process the effect is defintely stronger, while here the effect
# of the self-excitation is not strong enough to increase the intensity as in the Hawkes.
# Probably the solution (more complicated) is to do the opposite: try to build an 
# Hawkes process, that account for sender and receiver activity (multivariate Hawkes process)
