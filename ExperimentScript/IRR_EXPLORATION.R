# To calculate inter-rater reliability for the observed behavior onsets and offsets 
# between two raters in R, you can use the `irr` package. Specifically, you can 
# use the `irr::kappa2` function to compute Cohen's Kappa, which is suitable for 
# comparing two raters' categorical agreements.

# Script to calcualte IRR between two coded files. The user needs to provide 
# file paths to coded ELAN files exported to excel. The rest of the script runs 
# automatically, the output will be displayed in the console window. Note down 
# the average observed kappa value and report it in the table.

# Created by Dr Ola Dopierala, July 2023
# Updated August 2023 to allow user to drop variables (levels)
# Updated October 2023 to allow user to only assess first few minutes (E-CMS)


# Load the required package (if not installed, install first - install(packagename), e.g., install(dplyr))
library(dplyr)
library(irr)
library(readxl)

################################################################################
                        # NEEDS USER INPUT: 
################################################################################

# Exclude data from this timepoint forward (in minutes). Set to NULL if not excluding any timepoint.
ExclusionTimepoint_minutes <- 10  # e.g., to exclude data from 10 minutes onwards

# Convert timepoint in seconds to timepoint in minutes
# Do not change or enter any values here
ExclusionTimepoint_seconds <- ExclusionTimepoint_minutes * 60

if (!is.null(ExclusionTimepoint_minutes)) {
  ExclusionTimepoint_seconds <- ExclusionTimepoint_minutes * 60
} else {
  ExclusionTimepoint_seconds <- NULL
}


# Add the file paths to the two files you want to compare: 
# Rater 1 
rater1_data <- read_excel("EXPL005_Cam2_Sync_Motor_ DC_updated (Nov 2).xlsx", 
                          col_names = FALSE, skip=2) #first two rows don't contain data

# Rater 2 
rater2_data <- read_excel("EXPL005_Cam2_LL_Motor_New.xlsx", 
                          col_names = FALSE, skip=2) #first two rows don't contain data


# If you want to remove some variables (not run analyses on them), list them here 
# e.g., c('CIDS','TO')
#LevelsToExclude <- c('IG','IF' )# comment this line out by putting # before if you want to keep all levels! 


################################################################################
                        # Runs automatically
################################################################################

# checking number of columns and assign column names
# rater 1
if (ncol(rater1_data) == 6) {
colnames(rater1_data) <- c("beh", "blank", "onset", "offset", "duration", "notes")
} else {
colnames(rater1_data) <- c("beh", "blank", "onset", "offset", "duration")
}

if (!is.null(ExclusionTimepoint_seconds)) {
  rater1_data <- rater1_data %>% filter(onset < ExclusionTimepoint_seconds)
}

# rater 2
if (ncol(rater2_data) == 6) {
  colnames(rater2_data) <- c("beh", "blank", "onset", "offset", "duration", "notes")
} else {
  colnames(rater2_data) <- c("beh", "blank", "onset", "offset", "duration")
}

if (!is.null(ExclusionTimepoint_seconds)) {
  rater2_data <- rater2_data %>% filter(onset < ExclusionTimepoint_seconds)
}

# If exclusion timepoint is set, use it as total duration. Else, determine max duration from the coding files.
total_duration <- ifelse(!is.null(ExclusionTimepoint_seconds), ExclusionTimepoint_seconds, 
                         max(c(rater2_data$offset, rater1_data$offset)))


# getting data to right format
# rater 1
rater1_data$beh <- as.factor(rater1_data$beh)
rater1_data$onset <- as.numeric(rater1_data$onset)
rater1_data$offset <- as.numeric(rater1_data$offset)
# rater 2
rater2_data$beh <- as.factor(rater2_data$beh)
rater2_data$onset <- as.numeric(rater2_data$onset)
rater2_data$offset <- as.numeric(rater2_data$offset)

# Organising data into subsets depending on behaviour observed
# Rater 1
# Get unique levels of the 'beh' factor
unique_beh_levels1 <- unique(rater1_data$beh)

# Creating an empty list to store the subsets
subsets1 <- list()

# Looping through each level of the 'beh' factor
for (level in unique_beh_levels1) {
  # Create a subset based on the current level
  subset_data1 <- subset(rater1_data, beh == level)
  
  # Storing the subset in the 'subsets' list
  subsets1[[level]] <- subset_data1
}

# Rater 2
unique_beh_levels2 <- unique(rater1_data$beh)

# Creating an empty list to store the subsets
subsets2 <- list()

# Looping through each level of the 'beh' factor
for (level in unique_beh_levels2) {
  # Create a subset based on the current level
  subset_data2 <- subset(rater2_data, beh == level)
  
  # Storing the subset in the 'subsets' list
  subsets2[[level]] <- subset_data2
}


# Check if LevelsToExclude variable exists
if (exists("LevelsToExclude")) {
  unique_beh_levels1 <- droplevels(unique_beh_levels1, exclude = LevelsToExclude)
  unique_beh_levels2 <- droplevels(unique_beh_levels2, exclude = LevelsToExclude)
}


# CONVERT ONSETS AND OFFSETS INTO TIME SERIES

# Determining the max duration of either coding file
total_duration <- max(c(rater2_data$offset, rater1_data$offset))

# Creating a time vector with 100ms resolution
time_points <- seq(0, total_duration, by = 0.1)

# Initializing a list to store results for each level of 'beh'
kappa_results <- list()

# Initialize a vector to store the kappa values for each subset
kappa_values <- c()

# Looping through each level of the 'beh' factor
for (level in unique_beh_levels1) {
  # Subset data for the current level of 'beh'
  subset_rater1 <- subset(rater1_data, beh == level)
  subset_rater2 <- subset(rater2_data, beh == level)
  
  # RATER 1
  # Initializing a new data frame to store the time points and behavior observations
  rater1_observed_data <- data.frame(Time = time_points, BehaviorObserved = 0)
  
  # For each instance of the behavior, marking the corresponding time points as 1
  for (i in 1:nrow(subset_rater1)) {
    onset_time <- subset_rater1$onset[i]
    offset_time <- subset_rater1$offset[i]
    rater1_observed_data$BehaviorObserved[time_points >= onset_time & time_points <= offset_time] <- 1
  }
  
  # RATER 2
  # Initializing a new data frame to store the time points and behavior observations
  rater2_observed_data <- data.frame(Time = time_points, BehaviorObserved = 0)
  
  # For each instance of the behavior, marking the corresponding time points as 1
  for (i in 1:nrow(subset_rater2)) {
    onset_time <- subset_rater2$onset[i]
    offset_time <- subset_rater2$offset[i]
    rater2_observed_data$BehaviorObserved[time_points >= onset_time & time_points <= offset_time] <- 1
  }
  
  # Merging the onsets and offsets for both raters
  merged_data <- data.frame(
    rater1 = rater1_observed_data$BehaviorObserved,
    rater2 = rater2_observed_data$BehaviorObserved
  )
  
  # Calculating Cohen's Kappa
  kappa_result <- kappa2(merged_data)
  
  # Storing the result for the current level of 'beh' in the list
  kappa_results[[level]] <- kappa_result
}

# Displaying the results for each level of 'beh'
for (level in unique_beh_levels1) {
  print(paste("Kappa result for level", level, ":"))
  print(kappa_results[[level]])
  
  # Store the result for the current level of 'beh' in the vector
  kappa_values <- c(kappa_values, kappa_result$value)
  
}

# Calculate the average observed kappa value
average_kappa <- mean(kappa_values)

# Display the average kappa value
print(paste("Average observed kappa value:", average_kappa))

# The `kappa2` function will provide you with a Kappa value representing the 
# inter-rater reliability. The closer the Kappa value is to 1, the higher the 
# agreement between the raters.
