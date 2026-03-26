# Simulate REM with Pure Count (No Decay) ####
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
# Build initial risk set 
# --------------------------------------------

# Create all possible pairs
riskset <- cbind(rep(1:p, each = p), rep(1:p, p))

# Remove self-pairs (no one sends a letter to themselves)
riskset <- riskset[riskset[, 1] != riskset[, 2], ]

# Lenght of riskset
l <- nrow(riskset)


# --------------------------------------------
# Initialize covariates (Logic: Pure Count)
# --------------------------------------------

s_counts <- numeric(p)       # Count as sender

r_counts <- numeric(p)       # Count as receiver

rec_counts <- numeric(l)     # Count of dyads

# Initialize data 
dat <- NULL

# Initialize clock to 0
tm <- 0

# Decide the last event
event <- 0
total <- 5000

# Initialize first covariate
s_cov_init <- numeric(l)
r_cov_init <- numeric(l)
rec_cov_init <- numeric(l)

s_id <- riskset[, 1]
r_id <- riskset[, 2]

s_cov_init<- s_counts[s_id] 
r_cov_init<- r_counts[r_id] 
rec_cov_init <- rec_counts 

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
  
  # Save current values BEFORE update (these are the values that triggered the event)
  val_s <- s_counts[s] 
  val_r <- r_counts[r] 
  val_rec <- rec_counts[event.id] 
  
  dat <- rbind(dat, c(
    tm,
    s, r,
    val_s,
    val_r,
    val_rec
  ))
  
  # Update Counts
  s_counts[s] <- s_counts[s] + 1
  r_counts[r] <- r_counts[r] + 1
  rec_counts[event.id.rec] <- rec_counts[event.id.rec] + 1
  
  # Recalculate Hazard
  s_cov_next <- s_counts[riskset[, 1]]
  r_cov_next <- r_counts[riskset[, 2]]
  rec_cov_next <- rec_counts
  
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


##############################################################################
# STEP 2: PLOT
##############################################################################

# We can see that in this way the covariate does not decrease (not influence of 
# when the last occurence happened)
require(ggplot2)
require(dplyr)
df_dat <- as.data.frame(dat)
sender <- 14
events_node <- df_dat %>% filter(solicitor == sender)

ggplot(events_node, aes(x = tm, y = s_act_cov)) +
  geom_step(color = "red", size = 1) +
  geom_point(shape=1) +
  theme_minimal() +
  labs(
    title = paste("Only count", sender),
    subtitle = paste("count of events in which sender ", sender, " act as a sender"),
    x = "Time", 
    y = "Value covariate"
  )