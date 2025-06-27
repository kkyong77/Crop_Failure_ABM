# --- Load Required Libraries ---
library(dplyr)
library(purrr)
library(lubridate)
library(readr)
library(ranger)

# --- Load Supporting Scripts and Data ---
source("R/AgentClass.R")  # updated R6 class

agent_profile <- read_csv("data/agent_profile.csv")
agent_crop_history <- read_csv("data/agent_crop_history.csv")
agent_precip_list <- readRDS("data/precip/agent_precip_list.rds")

# --- Initialize Agents ---
agents <- vector("list", length = nrow(agent_profile))
names(agents) <- agent_profile$FBndID

for (i in seq_len(nrow(agent_profile))) {
  row <- agent_profile[i, ]
  
  agent <- Agent$new(
    id = row$FBndID,
    insurance_coverage = row$PER_INSUR,
    rotation_index = row$mean_RCI,
    lat = row$lat,
    lon = row$lon,
    acres = row$Acres
  )
  
  # Set precipitation time series
  agent$set_precip_ts(agent_precip_list[[row$FBndID]])
  
  # Add crop history
  for (yr in 10:22) {
    colname <- paste0("maj", yr)
    crop <- agent_crop_history[agent_crop_history$FBndID == row$FBndID, colname]
    if (length(crop) > 0) {
      agent$add_history(year = 2000 + yr, crop = crop[[1]])
    }
  }
  
  agents[[row$FBndID]] <- agent
}

# --- Simulation Loop ---
sim_years <- 2018:2020

crop_hist_all <- bind_rows(lapply(agents, function(agent) {
  do.call(rbind, lapply(agent$history, function(h) {
    data.frame(
      year = h$year,
      crop = h$crop,
      plantable = h$plantable %||% NA,
      plant_week = h$plant_week %||% NA,
      FBndID = agent$id,
      prob_plantable=NA,
      prob_max_rf = NA 
    )
  }))
}))

## clean up rownames
rownames(crop_hist_all)<-NULL

for (yr in sim_years) {
  message("Simulating year: ", yr)
  
  # --- Train model using all previous years (2011 to current - 1) ---
  crop_train <- crop_hist_all %>%
    filter(year < yr) %>%
    group_by(FBndID) %>%
    arrange(year) %>%
    mutate(
      Crop_Lag1 = lag(crop, 1),
      Crop_Lag2 = lag(crop, 2),
      Crop_Lag3 = lag(crop, 3)
    ) %>%
    ungroup() %>%
    left_join(agent_profile %>% select(FBndID, CropSumry), by = "FBndID") %>%
    filter(!is.na(Crop_Lag1), !is.na(Crop_Lag2), !is.na(Crop_Lag3)) %>%
    mutate(crop=as.factor(crop),
           Crop_Lag1=as.factor(Crop_Lag1),
           Crop_Lag2=as.factor(Crop_Lag2),
           Crop_Lag3=as.factor(Crop_Lag3),
           CropSumry=as.factor(CropSumry)
           )
  
  library(ranger)
  set.seed(123)
  
  rf_model <- ranger(
    formula = crop ~ Crop_Lag1 + Crop_Lag2 + Crop_Lag3 + CropSumry,
    data = crop_train %>%
      select(crop, Crop_Lag1, Crop_Lag2, Crop_Lag3, CropSumry),
    num.trees = 500,
    probability = TRUE,
    importance = "impurity"
  )
  
  total_agents <- length(agents)
  
  i <- 1
  for (fbnd_id in names(agents)) {
    agent <- agents[[fbnd_id]]
    
    message(" - Processing agent: ", fbnd_id)
    message(sprintf(" - Progress: %.1f%%", 100 * i / total_agents))
    i=i+1
    prev_crop1 <- agent$history[[as.character(yr - 1)]]$crop
    prev_crop2 <- agent$history[[as.character(yr - 2)]]$crop
    prev_crop3 <- agent$history[[as.character(yr - 3)]]$crop
    
    crop_summary <- agent_profile %>%
      filter(FBndID == fbnd_id) %>%
      pull(CropSumry)
    
    if (!any(is.na(c(prev_crop1, prev_crop2, prev_crop3))) && !is.null(agent$precip_ts)) {
      agent$simulate_year(
        year,
        rf_model,
        prev_y1_crop = prev_crop1,
        prev_y2_crop = prev_crop2,
        prev_y3_crop = prev_crop3,
        crop_summary = crop_summary
      )
    }
  }
  
  
  # --- Extract New Simulated Year Data ---
  new_hist <- bind_rows(lapply(agents, function(agent) {
    h <- agent$history[[as.character(yr)]]
    if (!is.null(h)) {
      df <- data.frame(
        year = h$year,
        crop = h$crop,
        plantable = h$plantable %||% NA,
        plant_week = h$plant_week %||% NA,
        FBndID = agent$id,
        prob_plantable=h,
        prob_plantable = h$prob_plantable %||% NA,
        prob_max_rf = h$prob_max_rf %||% NA
        
      )
    } else {
      NULL
    }
  }))
  
  # Append new year to cumulative history
  crop_hist_all <- bind_rows(crop_hist_all, new_hist)
 
}

# --- Save Output ---
write_csv(crop_history_all, "output/simulated_agent_crop_history.csv")
saveRDS(agents, "output/simulated_agents_list.rds")