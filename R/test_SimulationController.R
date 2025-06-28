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


set.seed(42)  # for reproducibility
sampled_agents <- agent_profile %>%
  sample_n(100)

# Filter only data for sampled agents
agent_crop_history <- agent_crop_history %>%
  filter(FBndID %in% sampled_agents$FBndID)

agent_precip_list <- agent_precip_list[sampled_agents$FBndID]

agents <- vector("list", length = nrow(sampled_agents))
names(agents) <- sampled_agents$FBndID

for (i in seq_len(nrow(sampled_agents))) {
  row <- sampled_agents[i, ]
  
  agent <- Agent$new(
    id = row$FBndID,
    insurance_coverage = row$PER_INSUR,
    rotation_index = row$mean_RCI,
    lat = row$lat,
    lon = row$lon,
    acres = row$Acres
  )
  
  agent$set_precip_ts(agent_precip_list[[row$FBndID]])
  
  for (yr in 10:22) {
    colname <- paste0("maj", yr)
    crop <- agent_crop_history[agent_crop_history$FBndID == row$FBndID, colname]
    if (length(crop) > 0) {
      agent$add_history(year = 2000 + yr, crop = crop[[1]])
    }
  }
  
  agents[[row$FBndID]] <- agent
}
## testing 
beta_0=1
beta_1= -0.12
b_p=500

# --- Simulation Loop ---
sim_years <- 2019:2019

crop_hist_all <- bind_rows(lapply(agents, function(agent) {
  do.call(rbind, lapply(agent$history, function(h) {
    data.frame(
      year = as.integer(h$year),
      crop = h$crop,
      plantable = h$plantable %||% NA,
      plant_week = as.numeric(h$plant_week %||% NA),
      FBndID = agent$id,
      prob_plantable=as.numeric(NA),
      prob_max_rf = as.numeric(NA) 
    )
  }))
}))
# 
# crop_hist_all <- crop_hist_all %>%
#   mutate(
#     year = as.integer(year),
#     plant_week=as.numeric(plant_week),
#     prob_plantable = as.numeric(prob_plantable),  # Coerce to double
#     prob_max_rf = as.numeric(prob_max_rf)  # Coerce to double
#     
#     
#   )

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
        yr,
        rf_model,
        prev_y1_crop = prev_crop1,
        prev_y2_crop = prev_crop2,
        prev_y3_crop = prev_crop3,
        crop_summary = crop_summary,
        beta_0 = beta_0,
        beta_1 =beta_1,
        b_p = b_p
      )
    }
  }
  
  
  # --- Extract New Simulated Year Data ---
  new_hist <- bind_rows(lapply(agents, function(agent) {
    h <- agent$history[[as.character(yr)]]
    if (!is.null(h)) {
      data.frame(
        year = h$year,
        crop = h$crop,
        plantable = h$plantable %||% NA,
        plant_week = h$plant_week %||% NA,
        FBndID = agent$id,
        prob_plantable = h$prob_plantable %||% NA,
        prob_max_rf = h$prob_max_rf %||% NA
      )
    } else {
      NULL
    }
  }))
  
  crop_hist_all <- crop_hist_all %>%
    mutate(year = as.integer(year))
  
   # Suppose each row is uniquely identified by agent_id and year
  crop_hist_all <- crop_hist_all %>%
    rows_update(new_hist, by = c("FBndID", "year"))
 
}

# --- Save Output ---
write_csv(crop_hist_all, "model/output/simulated_agent_crop_history.csv")
saveRDS(agents, "model/output/simulated_agents_list.rds")

library(ggplot2)
library(ggplot2)

crop_hist_all %>%
  filter(year == 2019) %>%
  ggplot(aes(x = crop, fill = crop)) +   # fill mapped to crop to get legend
  geom_bar() +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.3, size = 3) +
  theme_minimal() +
  labs(
    title = "Simulated Crop Distribution in 2019",
    x = "Crop Type",
    y = "Count",
    fill = "Crop Type"  # legend title
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"  # you can change to "bottom" or "top"
  )


