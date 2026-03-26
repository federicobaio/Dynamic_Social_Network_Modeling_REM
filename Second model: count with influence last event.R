# Simulate REM with Count and Decay Covariate ####
require(mgcv)

# Set the number of individuals in the system
p <- 50
set.seed(1234)

# Set coefficients according to the covariate type
b1 <- 0.01
b2 <- 0.01 # reduced values to avoid explosion
b3 <- 0 

# Lambda0 has an impact when dealing with the time-based definition
lambda0 <- 0.005 

##############################################################################
# STEP 1: Simulate the data
##############################################################################

# --------------------------------------------
# Build initial risk set (all possible directed pairs)
# --------------------------------------------

# Create all possible pairs
riskset <- cbind(rep(1:p, each = p), rep(1:p, p))

# Remove self-pairs (no one sends a letter to themselves)
riskset <- riskset[riskset[, 1] != riskset[, 2], ]

# Lenght of riskset
l <- nrow(riskset)

# --------------------------------------------
# Initialize covariates (Logic: Count + last event influence)
# --------------------------------------------

s_counts <- numeric(p)       # Count as sender
s_last_tm <- rep(-Inf, p)    # Last event

r_counts <- numeric(p)       # Count as receiver
r_last_tm <- rep(-Inf, p)    

rec_counts <- numeric(l)     # Count of dyads
rec_last_tm <- rep(-Inf, l)  

# Initialize data for simulated events
dat <- NULL

# Initialize clock to 0
tm <- 0

# Decide the last event
event <- 0
total <- 5000

s_cov_init <- numeric(l)
r_cov_init <- numeric(l)
rec_cov_init <- numeric(l)
for (k in 1:l){
  
  # extract sender e receiver k
  s_id <- riskset[k, 1]
  r_id <- riskset[k, 2]
  
  #Sender Activity
  s_cov_init[k] <- s_counts[s_id] *exp(-(tm-s_last_tm[s_id]))
  
  # Receiver Popularity
  r_cov_init[k] <- r_counts[r_id] *exp(-(tm-r_last_tm[r_id]))
  
  # For the reciprocity we directly use k
  rec_cov_init[k] <- rec_counts[k] * exp(-(tm - rec_last_tm[k]))
}

lp <- (b1 * s_cov_init + b2 * r_cov_init + b3 * rec_cov_init)
hazard <- exp(lp) * lambda0

# --------------------------------------------
# Keep going until we reach number of events = total
# --------------------------------------------

while (event < total) {
  
  # Step 1: Simulate time until next event 
  tm <- tm + rexp(1, sum(hazard))
  
  # Step 2: Select which dyad experiences the event
  event.id <- sample(length(hazard), 1, prob = hazard / sum(hazard))
  
  # Identify sender (s) and receiver (r)
  s <- riskset[event.id, 1]
  r <- riskset[event.id, 2]
  
  # Identify reciprocal pair index (receiver sending to sender)
  event.id.rec <- which(riskset[, 1] == r & riskset[, 2] == s)
  
  # Save current values before update (these are the values that triggered the event)
  val_s <- s_counts[s] * exp(-(tm-s_last_tm[s]))
  val_r <- r_counts[r] * exp(-(tm-r_last_tm[r]))
  val_rec <- rec_counts[event.id] * exp(-(tm-rec_last_tm[event.id]))
  
  dat <- rbind(dat, c(
    tm,
    s, r,
    val_s,
    val_r,
    val_rec
  ))
  
  # Update Sender
  s_counts[s] <- s_counts[s] + 1
  s_last_tm[s] <- tm
  
  # Update Receiver
  r_counts[r] <- r_counts[r] + 1
  r_last_tm[r] <- tm
  
  # Update Reciprocity (Note: Reciprocity tracks the previous r->s dyad)
  rec_counts[event.id.rec] <- rec_counts[event.id.rec] + 1
  rec_last_tm[event.id.rec] <- tm
  
  # Recalculate the hazard function with updated covariates
  s_cov_next <- numeric(l)
  r_cov_next <- numeric(l)
  rec_cov_next <- numeric(l)
  for (k in 1:l){
    
    # Extract sender e receiver k
    s_id <- riskset[k, 1]
    r_id <- riskset[k, 2]
    # Sender Activity
    s_cov_next[k] <- s_counts[s_id] *exp(-(tm-s_last_tm[s_id]))
    # Receiver Popularity
    r_cov_next[k] <- r_counts[r_id] *exp(-(tm-r_last_tm[r_id]))
    
    # For the reciprocity we directly use k
    rec_cov_next[k] <- rec_counts[k] * exp(-(tm - rec_last_tm[k]))
  }
  
  lp <- (b1 * s_cov_next + b2 * r_cov_next + b3 * rec_cov_next)
  hazard <- exp(lp) * lambda0
  
  event <- event + 1
}

colnames(dat) <- c("tm", 
                   "solicitor", 
                   "receiver", 
                   "s_act_cov", 
                   "r_pop_cov", 
                   "rec_cov")
n.event <- nrow(dat)
dat
tail(dat)
# Compare how the events are spread
division<- table(dat[,"solicitor"])
division




##############################################################################
# STEP 2: PLOT 
##############################################################################

require(dplyr)
require(ggplot2)

df_dat <- as.data.frame(dat)

# We select our sender
sender <- 14
# Extract events in which the node selected act as a sender
events_node <- df_dat %>% filter(solicitor == sender)
event_times <- sort(events_node$tm)

# Create the grid needed to get a smooth plot
time_grid <- seq(0, max(df_dat$tm), length.out = 5000)
decay_curve <- numeric(length(time_grid))
no_decay_curve <- numeric(length(time_grid))

curr_count <- 0
last_t_val <- -Inf
ev_idx <- 1

# Loop to get the continous line
for(j in 1:length(time_grid)) {
  t_val <- time_grid[j]
  
  while(ev_idx <= length(event_times) && event_times[ev_idx] <= t_val) {
    curr_count <- curr_count + 1
    last_t_val <- event_times[ev_idx]
    ev_idx <- ev_idx + 1
  }
  
  # Step curve
  no_decay_curve[j] <- curr_count
  
  # Decay effect
  if(last_t_val == -Inf) {
    decay_curve[j] <- 0
  } else {
    dt <- t_val - last_t_val
    decay_curve[j] <- curr_count * exp(-dt)
  }
}

df_compare <- data.frame(Time = time_grid, With_Decay = decay_curve, No_Decay = no_decay_curve)
max_time_zoom <- min(max(df_dat$tm), 50) 

ggplot() +
  # Red: Only count
  geom_line(data = df_compare %>% filter(Time < max_time_zoom), 
            aes(x = Time, y = No_Decay, color = "Count"), 
            size = 0.8, linetype="dashed") +
  # Blu: Real value covariate
  geom_line(data = df_compare %>% filter(Time < max_time_zoom), 
            aes(x = Time, y = With_Decay, color = "Count * exp(-t)"), 
            size = 1) +
  scale_color_manual(values = c("Count" = "red", "Count * exp(-t)" = "blue")) +
  theme_minimal() +
  labs(
    title = paste("Zoomed covariate value for sender", sender),
    subtitle = "Fast decay (implicit Beta=1)",
    y = "Covariate value", x = "Time"
  ) +
  theme(legend.position = "bottom")