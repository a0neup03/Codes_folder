if (!requireNamespace("survival", quietly = TRUE)) {
  install.packages("survival")
}
library(survival)

# Function to calculate required events using Schoenfeld's method
calculate_events <- function(power, alpha, HR, p) {
  Z_alpha <- qnorm(1 - alpha / 2)
  Z_beta <- qnorm(power)
  P <- p * (1-p)
  events <- ((Z_alpha + Z_beta)^2) / (P * (log(HR))^2)
  return(ceiling(events))
}

# Function to calculate total sample size based on event rate
calculate_sample_size <- function(events, event_rate) {
  sample_size <- events / event_rate
  return(ceiling(sample_size))
}

# Create results dataframe function
create_results_df <- function(type = "events") {
  # Parameters
  alpha <- 0.05
  powers <- c(0.70, 0.80, 0.90)
  HR_values <- c(1.5, 2.0, 3.0)
  prevalence_values <- c(0.2, 0.3, 0.4)
  event_rates <- c(0.4,0.5, 0.7)
  
  # Create empty dataframe
  results <- data.frame(
    Power = numeric(),
    Prevalence = numeric(),
    HR_1.5 = numeric(),
    HR_2.0 = numeric(),
    HR_3.0 = numeric()
  )
  
  # Calculate results
  for (power in powers) {
    for (prev in prevalence_values) {
      row <- c(power * 100, prev * 100)  # Convert to percentages
      for (hr in HR_values) {
        events <- calculate_events(power, alpha, hr, prev)
        if (type == "events") {
          row <- c(row, events)
        } else if (type == "sample_size_40") {
          row <- c(row, calculate_sample_size(events, 0.4))
        } else if (type == "sample_size_70") {
          row <- c(row, calculate_sample_size(events, 0.7))
        }
        else if (type=="sample_size_50"){
          row <- c(row, calculate_sample_size(events, 0.5))
        }
      }
      results <- rbind(results, row)
    }
  }
  
  colnames(results) <- c("Power", "Prevalence", "HR = 1.5", "HR = 2.0", "HR = 3.0")
  return(results)
}

events_table <- create_results_df("events")
sample_size_40 <- create_results_df("sample_size_40")
sample_size_50 <- create_results_df("sample_size_50")

sample_size_70 <- create_results_df("sample_size_70")

print_formatted_table <- function(df, title) {
  cat("\n", title, "\n")
  print(df, row.names = FALSE)
}

print_formatted_table(events_table, 
                      "Table 1: Required Number of Events for Different Power Levels by Hazard Ratio and Biomarker Prevalence")
print_formatted_table(sample_size_40, 
                      "\nTable 2A: Required Sample Size for Different Power Levels (40% Event Rate)")
print_formatted_table(sample_size_70, 
                      "\nTable 2B: Required Sample Size for Different Power Levels (70% Event Rate)")

write.csv(events_table, "required_events.csv", row.names = FALSE)
write.csv(sample_size_40, "sample_size_40_events.csv", row.names = FALSE)
write.csv(sample_size_70, "sample_size_70_events.csv", row.names = FALSE)



#second question

# Parameters
alpha <- 0.05
powers <- c(0.70, 0.80, 0.90)     # Three power levels
HR_values <- c(1.5, 2.0, 3.0)
prevalence_values <- c(0.2, 0.3, 0.4)

# Create list to store results for each power level
all_results <- list()

for (power in powers) {
  # Create matrix for events
  results_events <- matrix(NA, nrow = length(prevalence_values), ncol = length(HR_values))
  
  # Calculate events for all combinations
  for (i in seq_along(prevalence_values)) {
    for (j in seq_along(HR_values)) {
      # Calculate required events
      results_events[i,j] <- calculate_events(power, alpha, HR_values[j], prevalence_values[i])
    }
  }
  
  # Format results as data frame
  results_events <- as.data.frame(results_events)
  
  # Add row and column names
  colnames(results_events) <- paste("HR =", HR_values)
  rownames(results_events) <- paste0(prevalence_values * 100, "% Prevalence")
  
  # Store results
  all_results[[paste0("power_", power)]] <- results_events
}

# Print results for all power levels
for (power in powers) {
  cat("\n\nRequired Number of Events for Power =", power * 100, "%\n")
  print(all_results[[paste0("power_", power)]])
}

# Create a summary table comparing events across power levels
cat("\n\nComparison of Required Events Across Power Levels:\n")
cat("\nFor 20% Prevalence:\n")
events_comparison <- data.frame(
  "Power" = powers,
  "HR_1.5" = sapply(powers, function(p) all_results[[paste0("power_", p)]][1,1]),
  "HR_2.0" = sapply(powers, function(p) all_results[[paste0("power_", p)]][1,2]),
  "HR_3.0" = sapply(powers, function(p) all_results[[paste0("power_", p)]][1,3])
)
print(events_comparison)

# Calculate percentage increase in events needed
cat("\nPercentage Increase in Events from 70% to 90% Power:\n")
percent_increase <- function(low_power, high_power, hr_col) {
  increase <- (events_comparison[events_comparison$Power == high_power, hr_col] - 
                 events_comparison[events_comparison$Power == low_power, hr_col]) / 
    events_comparison[events_comparison$Power == low_power, hr_col] * 100
  return(round(increase, 1))
}

for (hr in c("HR_1.5", "HR_2.0", "HR_3.0")) {
  cat(sprintf("\nFor %s: %.1f%% increase", 
              gsub("_", " = ", hr), 
              percent_increase(0.7, 0.9, hr)))
}

