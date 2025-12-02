library(dplyr)
library(lubridate)
library(readr)
library(ranger)
library(rlang)
library(tidyr)

source("R/AgentClass_model2.R")  # Your Agent class with simulate_year method

# Load your data
agent_profile <- read_csv("data/agent_profile.csv")
agent_crop_history <- read_csv("data/agent_crop_history.csv")
agent_precip_list <- readRDS("data/precip/agent_precip_list.rds")
agent_moist_list <- readRDS("data/climate_moist/agent_moist_list.rds")

# Sample agents if needed (or use all)
set.seed(42)
sampled_agents <- agent_profile 
#%>% sample_n(1000)

agent_crop_history <- agent_crop_history %>%
  filter(FBndID %in% sampled_agents$FBndID)
agent_moist_list <- agent_moist_list[sampled_agents$FBndID]

# Initialize Agents
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
  agent$set_moist_ts(agent_moist_list[[row$FBndID]])
  
  for (yr in 10:22) {
    colname <- paste0("maj", yr)
    crop <- agent_crop_history[agent_crop_history$FBndID == row$FBndID, colname]
    rf_crop<-"NA"
    if (length(crop) > 0) agent$add_history(year = 2000 + yr, crop = crop[[1]])
  }
  agents[[row$FBndID]] <- agent
}

# Parameters
beta_0 <- 1
beta_1 <- -0.1
b_d <- -5

sim_years <- 2018:2019

# Initialize crop history dataframe
crop_hist_all <- bind_rows(lapply(agents, function(agent) {
  do.call(rbind, lapply(agent$history, function(h) {
    data.frame(
      year = as.integer(h$year),
      crop = h$crop,
      rf_crop=as.character(NA),
      acres=as.numeric(h$acres),
      plantable = h$plantable %||% NA,
      plant_week = as.numeric(h$plant_week %||% NA),
      FBndID = agent$id,
      prob_plantable=as.numeric(NA),
      prob_max_rf = as.numeric(NA) 
    )
  }))
}))
# 

rownames(crop_hist_all) <- NULL

for (yr in sim_years) {
  message("Simulating year: ", yr)
  
  # Train RF model on all previous years
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
    mutate(
      crop = as.factor(crop),
      Crop_Lag1 = as.factor(Crop_Lag1),
      Crop_Lag2 = as.factor(Crop_Lag2),
      Crop_Lag3 = as.factor(Crop_Lag3),
      CropSumry = as.factor(CropSumry)
    )
  
  set.seed(123)
  
  rf_model <- ranger(
    formula = crop ~ Crop_Lag1 + Crop_Lag2 + Crop_Lag3 + CropSumry,
    data = crop_train %>% select(crop, Crop_Lag1, Crop_Lag2, Crop_Lag3, CropSumry),
    num.trees = 500,
    probability = TRUE,
    importance = "impurity"
  )
  
  # Step 1: Predict crop and planting probabilities for all agents for the year
  pred_list <- list()
  
  for (fbnd_id in names(agents)) {
    agent <- agents[[fbnd_id]]
    
    prev_crop1 <- agent$history[[as.character(yr - 1)]]$crop
    prev_crop2 <- agent$history[[as.character(yr - 2)]]$crop
    prev_crop3 <- agent$history[[as.character(yr - 3)]]$crop
    
    crop_summary <- agent_profile %>%
      filter(FBndID == fbnd_id) %>%
      pull(CropSumry)
    
    if (!any(is.na(c(prev_crop1, prev_crop2, prev_crop3))) && !is.null(agent$moist_ts)) {
      # Prepare input for RF prediction
      input <- data.frame(
        Crop_Lag1 = factor(prev_crop1, levels = levels(crop_train$Crop_Lag1)),
        Crop_Lag2 = factor(prev_crop2, levels = levels(crop_train$Crop_Lag2)),
        Crop_Lag3 = factor(prev_crop3, levels = levels(crop_train$Crop_Lag3)),
        CropSumry = factor(crop_summary, levels = levels(crop_train$CropSumry))
      )
      
      rf_probs <- predict(rf_model, data = input, type = "response")$predictions[1, ]
      crop_pred <- names(which.max(rf_probs))
      agent
      planting_windows <- list(
        "C" = 14:26,  # corn weeks
        "B" = 16:28   # soybean weeks
      )
      
      # If predicted crop is corn or soybean (C or B), compute plantability each week
      if (crop_pred %in% c("C", "B")) {
        weeks_to_check <- planting_windows[[crop_pred]]
        
        moist_df <- agent$moist_ts %>%
          filter(lubridate::year(date) == yr) %>%
          mutate(week = lubridate::week(date)) %>%
          group_by(week) %>%
          summarise(weekly_moist = sum(moist), .groups = "drop")
        
        for (w in weeks_to_check) {
          p1 <- moist_df %>% filter(week == w - 1) %>% pull(weekly_moist)
          p2 <- moist_df %>% filter(week == w) %>% pull(weekly_moist)
          p1 <- ifelse(length(p1) == 0, 0, p1)
          p2 <- ifelse(length(p2) == 0, 0, p2)
          twomoist <- p1 + p2
          prob_plantable <- 1 / (1 + exp(-(beta_0 + beta_1 * (twomoist - b_d))))
          
          pred_list[[paste(fbnd_id, w, sep = "_")]] <- data.frame(
            FBndID = fbnd_id,
            crop_pred = crop_pred,
            prob_plantable = prob_plantable,
            plant_week = w,
            max_rf_prob = max(rf_probs),
            acres = agent$acres,
            stringsAsFactors = FALSE
          )
        }
      } else {
        # For other crops, just assign planting at the first week of a default window (e.g., 14)
        pred_list[[fbnd_id]] <- data.frame(
          FBndID = fbnd_id,
          crop_pred = crop_pred,
          prob_plantable = NA,   # no probability calculated
          plant_week = NA,       # assume  don't know 
          max_rf_prob = max(rf_probs),
          acres = agent$acres,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  
  # Bind all rows (for all weeks for corn/soybean, and single row for others)
  crop_week_df <- do.call(rbind, pred_list)
  
  # Step 2: Apply planting constraints per crop and week using actual acres
  # Define max percent area per crop per week allowed to plant
  #max_weekly_percent <- 0.1  # example: max 15% of the area of each crop acres can plant per week
  
  ## for C/B
  selected_agents <- list()
  
  ## except for C/B
  unselected_agents <- list()
  
  crops <- unique(crop_week_df$crop_pred)
  weeks <- unique(crop_week_df$plant_week)
  weeks <- weeks[!is.na(weeks)]
  
  
  for (crop_i in crops) {
    
    crop_df <- crop_week_df %>% 
      filter(crop_pred == crop_i)
    
    if (!(crop_i %in% names(planting_windows))) {
      unselected_agents[[paste(crop_i)]] <- crop_df
      
    }
    
    else{
      
      ## initializing m and s
      # mu=22
      # s=1.3
      
      mu <- list(
        "C" = 20,  # corn weeks
        "B" = 22   # soybean weeks
      )
      
      s <- list(
        "C" = 0.88,  # corn weeks
        "B" =1.04    # soybean weeks
      )
 
      # Step 1: Calculate logistic CDF over the planting weeks
      weeks <- planting_windows[[crop_i]]
      cdf_values <- 1 / (1 + exp(-(weeks - mu[[crop_i]]) / s[[crop_i]]))
      
      # Step 2: Convert to PDF (difference between CDF values)
      pdf_values <- c(cdf_values[1], diff(cdf_values))  # First week is just CDF value
      
      # Step 3: Normalize PDF so it sums to 1 (if needed)
      pdf_values <- pdf_values / sum(pdf_values)
      
      # Step 4: Multiply by total available acres for that crop
      total_acres_crop <- crop_df %>%
        distinct(FBndID, .keep_all = TRUE) %>%
        summarise(total_acres = sum(acres, na.rm = TRUE)) %>%
        pull(total_acres)
      
      weekly_acres_allowed <- pdf_values * total_acres_crop
      
      names(weekly_acres_allowed) <- weeks
      
      
      already_selected_ids <- c()
      
      message("Crop Type: ", crop_i)
      
      for (w in weeks) {
        week_df <- crop_df %>%
          filter(plant_week == w, !(FBndID %in% already_selected_ids)) %>%
          filter(prob_plantable > 0.5) %>%
          arrange(desc(prob_plantable)) %>%
          mutate(cumsum_acres = cumsum(acres)) %>%
          filter(cumsum_acres <= weekly_acres_allowed[as.character(w)])
        
         message("Week ", w, " — Selected: ", nrow(week_df), " — Limit: ", round(weekly_acres_allowed[as.character(w)], 2))
        
        already_selected_ids <- union(already_selected_ids, week_df$FBndID)
        
        week_df <- week_df %>%
          mutate(plant_week = w)
        
        selected_agents[[paste(crop_i, w, sep = "_")]] <- week_df
      }
      
    }
  }
  
  ## saving the planting decision for C/B
  selected_df <- bind_rows(selected_agents) %>%
    select("FBndID","crop_pred","prob_plantable","plant_week","max_rf_prob","acres")
  
  any(duplicated(selected_df$FBndID))  # Returns TRUE if duplicates exist
  
  rownames(selected_df) <- selected_df$FBndID
  
  unselected_df <- bind_rows(unselected_agents)
  
  # Now rbind safely
  combined_df <- rbind(selected_df, unselected_df)
  
  # Update agent histories
  for (fbnd_id in names(agents)) {
    agent <- agents[[fbnd_id]]
    pred_row <- combined_df %>% filter(FBndID == fbnd_id)
    
    if (nrow(pred_row) == 1) {
      agent$history[[as.character(yr)]] <- list(
        year = yr,
        crop = pred_row$crop_pred,
        rf_crop=pred_row$crop_pred,
        plant_week = pred_row$plant_week,
        plantable = TRUE,
        prob_plantable = pred_row$prob_plantable,
        prob_max_rf = pred_row$max_rf_prob,
        acres = agent$acres
      )
    } else{
      rf_pred <- bind_rows(pred_list) %>%
        filter(FBndID == fbnd_id)
      # Not selected to plant, mark as idle
      agent$history[[as.character(yr)]] <- list(
        year = yr,
        crop = "I",
        rf_crop=rf_pred$crop_pred[1],
        plant_week = NA,
        plantable = FALSE,
        prob_plantable = NA,
        prob_max_rf = rf_pred$max_rf_prob[1],
        acres = agent$acres
      )
    }
  }
  
  # Construct new_hist after agent updates
  new_hist <- bind_rows(lapply(agents, function(agent) {
    h <- agent$history[[as.character(yr)]]
    if (!is.null(h)) {
      as.data.frame(h, stringsAsFactors = FALSE) %>% mutate(FBndID = agent$id)
    } else {
      NULL
    }
  }))
  
  # Ensure year is integer (to match existing crop_hist_all format)
  new_hist <- new_hist %>%
    mutate(year = as.integer(year))
  
  # Update existing records
  crop_hist_all <- crop_hist_all %>%
    rows_update(new_hist, by = c("FBndID", "year"))
  
  
}

# --- Save Outputs ---
write_csv(crop_hist_all, "model/output/model2_1_simulated_agent_crop_history.csv")
saveRDS(agents, "model/output/model2_1_simulated_agents_list.rds")



## model evaluation ####
## planted area per week for corn and soybean by year #################

library(dplyr)


crop_hist_sim<-crop_hist_all %>%
  filter(year %in%sim_years)


crop_c_b_cum <- crop_hist_sim %>%
  filter(crop %in% c("C", "B"), plantable == TRUE) %>%
  group_by(year, crop, plant_week) %>%
  summarise(weekly_acres = sum(acres, na.rm = TRUE), .groups = "drop") %>%
  arrange(year, crop, plant_week) %>%
  group_by(year, crop) %>%
  mutate(
    cum_acres = cumsum(weekly_acres),
    total_acres = sum(weekly_acres),
    cum_pct = cum_acres / total_acres * 100
  ) %>%
  ungroup()

crop_c_b_cum <- crop_c_b_cum %>%
  mutate(source = "Simulated") %>%
  select(year, crop, plant_week, cum_pct, source)

library(ggplot2)

ggplot(crop_c_b_cum, aes(x = plant_week, y = cum_pct, color = as.factor(year))) +
  geom_line(size = 1) +
  geom_point(size = 2, alpha = 0.7) +
  facet_wrap(~ crop) +
  labs(
    title = "Cumulative Percent Planted by Week",
    x = "Plant Week",
    y = "Cumulative Percent (%)",
    color = "Crop"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "bottom"
  )


# Prepare observed corn data ##############
data_corn_cum <- data_corn %>%
  filter(year %in% c(2018, 2019), !is.na(weekNum), !is.na(percentPlanted)) %>%
  mutate(
    crop = "C",  # Assuming "C" is the label for corn in your model
    source = "Observed",
    cum_pct = percentPlanted,
    plant_week = weekNum
  ) %>%
  select(year, crop, plant_week, cum_pct, source)

# Prepare observed soybean data ##############
data_soy_cum <- data_soy %>%
  filter(year %in% c(2018, 2019), !is.na(weekNum), !is.na(percentPlanted)) %>%
  mutate(
    crop = "B",  # Assuming "C" is the label for corn in your model
    source = "Observed",
    cum_pct = percentPlanted,
    plant_week = weekNum
  ) %>%
  select(year, crop, plant_week, cum_pct, source)


## combining the model/observed planted area 

combined_c_b_cdf <- bind_rows(data_corn_cum,data_soy_cum, crop_c_b_cum)

p1<-ggplot(combined_c_b_cdf %>% filter(crop == "C"), 
           aes(x = plant_week, y = cum_pct, color = factor(year), linetype = source)) +
  geom_line(size = 1.2) +
  geom_point(size = 2, alpha = 0.7, aes(shape = source)) +
  labs(
    title = "Cumulative Corn Planting Progress",
    x = "Planting Week",
    y = "Cumulative Percent Planted",
    color = "Year",
    linetype = "Data Source",
    shape = "Data Source"
  ) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "bottom"
  )

print(p1)

p2<-ggplot(combined_c_b_cdf %>% filter(crop == "B"), 
           aes(x = plant_week, y = cum_pct, color = factor(year), linetype = source)) +
  geom_line(size = 1.2) +
  geom_point(size = 2, alpha = 0.7, aes(shape = source)) +
  labs(
    title = "Cumulative Soybean Planting Progress",
    x = "Planting Week",
    y = "Cumulative Percent Planted",
    color = "Year",
    linetype = "Data Source",
    shape = "Data Source"
  ) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "bottom"
  )

## print the model comparison 
library(cowplot)


# Combine plots side by side with labels
combined_plots <- plot_grid(p1, p2, labels = c("A", "B"), ncol = 2, align = "hv")

print(combined_plots)

ggsave("model/figures/model2_1_crop_planting_comparison_sim_obs_2018_2019.jpg",
       plot = combined_plots,
       width = 14, height = 6, dpi = 300)


compute_rmse <- function(sim, obs) {
  sqrt(mean((sim - obs)^2, na.rm = TRUE))
}


cdf_wide <- combined_c_b_cdf %>%
  filter(crop == "C", year == 2018) %>%
  pivot_wider(names_from = source, values_from = cum_pct)

rmse <- compute_rmse(cdf_wide$Simulated, cdf_wide$Observed)

cat("Corn,RMSE:", rmse, "\n")


## compute the observed crop type data ######
crop_obs_long <- agent_crop_history %>%
  pivot_longer(cols = maj10:maj22, names_to = "Year", values_to = "Crop_Type") %>%
  mutate(
    Year = as.integer(sub("maj", "", Year)) + 2000,
    Crop_Type = as.factor(Crop_Type)
  ) %>%
  select(FBndID, Year, Crop_Type, Acres)

# Compute percent acreage by crop type and year
crop_obs_percent <- crop_obs_long %>%
  group_by(Year, Crop_Type) %>%
  summarise(Total_Acres = sum(Acres, na.rm = TRUE), .groups = "drop") %>%
  group_by(Year) %>%
  mutate(Percent = 100 * Total_Acres / sum(Total_Acres)) %>%
  ungroup()

obs_df <- crop_obs_percent %>%
  rename(year = Year, crop = Crop_Type, percent = Percent) %>%
  mutate(source = "Observed")


library(dplyr)

crop_area_percent <- crop_hist_sim %>%
  group_by(year, crop) %>%
  summarise(total_acres = sum(acres, na.rm = TRUE), .groups = "drop_last") %>%
  mutate(year_total = sum(total_acres)) %>%
  ungroup() %>%
  mutate(percent = 100 * total_acres / year_total)

sim_df <- crop_area_percent %>%
  rename(year = year, crop =crop, percent = percent) %>%
  mutate(source = "Simulated")

combined_df <- bind_rows(obs_df, sim_df)

p3<-ggplot(combined_df %>% filter(year >= 2010), 
           aes(x = year, y = percent, color = crop, linetype = source)) +
  geom_line(linewidth = 1) +
  geom_point(aes(shape = source), size = 2) +
  scale_color_manual(
    values = c(
      "B" = "green", 
      "C" = "#FFD700",  # bright yellow-gold
      "W" = "orange", 
      "P" = "skyblue",
      "R" = "black",
      "O" = "grey",
      "I" = "brown"
    )
  ) +
  scale_x_continuous(breaks = seq(2010, max(combined_df$year), by = 2)) +  # here!
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  labs(
    title = "Observed vs Simulated Percent Area by Crop",
    x = "Year",
    y = "Percent of Area",
    color = "Crop Type",
    linetype = "Source",
    shape = "Source"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(p3)

ggsave("model/figures/model2_1_crop_percent_comparison_sim_obs_2018_2019.jpg",
       plot = p3,
       width = 14, height = 6, dpi = 300)


library(ggplot2)
p4<-crop_hist_all %>%
  filter(year %in%sim_years) %>%
  ggplot(aes(x=crop, fill = crop)) +   # fill mapped to crop to get legend
  geom_bar() +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.3, size = 3) +
  scale_fill_manual(
    values = c(
      "B" = "green", 
      "C" = "yellow", 
      "W" = "orange", 
      "P" = "skyblue",
      "R" = "black",
      "O" = "grey",
      "I" = "brown"
    )
  ) +  # <-- added plus sign here
  facet_wrap(~ year) +
  theme_minimal() +
  labs(
    title = "Simulated Crop Distribution",
    x = "Crop Type",
    y = "Count",
    fill = "Crop Type"  # legend title
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

print(p4)

ggsave("model/figures/model2_1_crop_pixel_comparison_sim_obs_2018_2019.jpg",
       plot = p4,
       width = 14, height = 6, dpi = 300)


## reading model output ###

crop_hist_sim<-read_csv("model/crop_hist_sim_model2.1.csv")

## spatial mapping of the crop choice 
library(sf)
library(dplyr)
library(ggplot2)

# Load parcel shapefile
acpf_prb_sf <- st_read("gis/prb_merged_acpf.shp")

# Crop fill colors
crop_colors <- c(
  "B" = "green", 
  "C" = "#FFD700",  # bright yellow-gold
  "W" = "orange", 
  "P" = "skyblue",
  "R" = "black",
  "O" = "grey",
  "I" = "brown"
)

# Years to map
predict_years <- 2018:2019

# Loop over years
for (yr in predict_years) {
  data_yr <- crop_hist_sim %>% filter(year == yr)
  
  # Merge shapefile with crop prediction
  map_data <- left_join(acpf_prb_sf, data_yr, by = "FBndID")
  
  # Plot
  p <- ggplot(map_data) +
    geom_sf(aes(fill = crop), color = NA) +
    scale_fill_manual(values = crop_colors, na.value = "white") +
    theme_minimal() +
    labs(
      title = paste("Simulated Crop Choice -", yr),
      x = "Longitude", y = "Latitude", fill = "Crop Type"
    )
  
  print(p)
  
  ggsave(
    filename = paste0("model/figures/crop_choice_model_", yr, ".jpg"),
    plot = p, width = 8, height = 6, dpi = 300
  )
}


## spatial mapping of the planting week 
library(viridis)  # for color scales

# Loop over years
for (yr in predict_years) {
  data_yr <- crop_hist_sim %>% filter(year == yr)
  
  # Merge shapefile with crop prediction
  map_data <- left_join(acpf_prb_sf, data_yr, by = "FBndID")
  
  # Plot planting week with distinct color
  p <- ggplot(map_data) +
    geom_sf(aes(fill = factor(plant_week)), color = NA) +
    scale_fill_viridis_d(
      option = "plasma",  # or "turbo", "magma", "viridis"
      na.value = "grey90",
      name = "Planting Week"
    ) +
    theme_minimal() +
    labs(
      title = paste("Planting Week of Crop Choice -", yr),
      x = "Longitude", y = "Latitude"
    )
  
  print(p)
  
  ggsave(
    filename = paste0("model/figures/crop_planting_weeks_model_", yr, ".jpg"),
    plot = p, width = 8, height = 6, dpi = 300
  )
}




### observed crop data  ###
library(sf)
library(dplyr)
library(ggplot2)

# Load parcel shapefile
acpf_prb_sf <- st_read("gis/prb_merged_acpf.shp")

# Crop fill colors
crop_colors <- c(
  "B" = "green", 
  "C" = "#FFD700",  # bright yellow-gold
  "W" = "orange", 
  "P" = "skyblue",
  "R" = "black",
  "O" = "grey",
  "I" = "brown"
)

# Years to map
predict_years <- 2018:2019

# Loop over years
for (yr in predict_years) {
  data_yr <- crop_hist_sim %>% filter(year == yr)
  
  # Merge shapefile with crop prediction
  map_data <- left_join(acpf_prb_sf, data_yr, by = "FBndID")
  
  # Plot
  p <- ggplot(map_data) +
    geom_sf(aes(fill = crop), color = NA) +
    scale_fill_manual(values = crop_colors, na.value = "white") +
    theme_minimal() +
    labs(
      title = paste("Simulated Crop Choice -", yr),
      x = "Longitude", y = "Latitude", fill = "Crop Type"
    )
  
  print(p)
  
  ggsave(
    filename = paste0("model/figures/crop_choice_model_", yr, ".jpg"),
    plot = p, width = 8, height = 6, dpi = 300
  )
}


library(sf)
library(dplyr)
library(ggplot2)

# Load parcel shapefile
acpf_prb_sf <- st_read("gis/prb_merged_acpf.shp")

# Crop fill colors
crop_colors <- c(
  "B" = "green", 
  "C" = "#FFD700",  # bright yellow-gold
  "W" = "orange", 
  "P" = "skyblue",
  "R" = "black",
  "O" = "grey",
  "I" = "brown"
)


crop_hist<-read_csv("data/acpf_prb_reclassfied_0523.csv")

crop_data <- crop_hist %>%
  dplyr::select(FBndID,lat,lon, maj10:maj22,CropSumry)

# Years to map
# Loop over years
  # Merge shapefile with crop prediction
  hist_crop_map <-merge(acpf_prb_sf[,"FBndID"], crop_data, by = "FBndID") %>%
                  select("FBndID","maj18","maj19")
  
  # Plot
  p <- ggplot(hist_crop_map) +
    geom_sf(aes(fill = maj19), color = NA) +
    scale_fill_manual(values = crop_colors, na.value = "white") +
    theme_minimal() +
    labs(
      title = paste("Observed Crop Choice -", 2019),
      x = "Longitude", y = "Latitude", fill = "Crop Type"
    )
  
  print(p)
  
  ggsave(filename = paste0("model/figures/hist_crop_choice_", 2019, ".jpg"),
    plot = p, width = 8, height = 6, dpi = 300)


