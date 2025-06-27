Agent <- R6::R6Class("Agent",
                     public = list(
                       id = NULL,
                       insurance_coverage = NULL,
                       rotation_index = NULL,
                       lat = NULL,
                       lon = NULL,
                       acres = NULL,
                       history = list(),
                       precip_ts = NULL,
                       
                       initialize = function(id, insurance_coverage, rotation_index, lat = NA, lon = NA, acres = NA) {
                         self$id <- id
                         self$insurance_coverage <- insurance_coverage
                         self$rotation_index <- rotation_index
                         self$lat <- lat
                         self$lon <- lon
                         self$acres <- acres
                       },
                       
                       add_history = function(year, crop) {
                         self$history[[as.character(year)]] <- list(year = year, crop = crop)
                       },
                       
                       get_history_df = function() {
                         do.call(rbind, lapply(self$history, function(h) as.data.frame(h)))
                       },
                       
                       set_precip_ts = function(precip_df) {
                         self$precip_ts <- precip_df
                       },
                       
                       logistic_prob = function(precip, beta_0 = 1, beta_1 = -0.1, b_p = 250) {
                         logit <- beta_0 + beta_1 * (precip - b_p)
                         prob <- 1 / (1 + exp(-logit))
                         return(prob)
                       },
                       
                       simulate_year = function(year, rf_model, prev_y1_crop, prev_y2_crop, prev_y3_crop, crop_summary,
                                                beta_0 = 1, beta_1 = -0.1, b_p = 250,
                                                planting_start_week = 14,
                                                planting_end_week = 26) {
                         
                         precip_window <- self$precip_ts %>%
                           filter(lubridate::year(date) == year) %>%
                           mutate(week = lubridate::week(date)) %>%
                           filter(week >= (planting_start_week - 1) & week <= planting_end_week) %>%
                           group_by(week) %>%
                           summarise(weekly_precip = sum(precip), .groups = "drop") %>%
                           arrange(week)
                         
                         plantable <- FALSE
                         plant_week <- NA
                         candidate_weeks <- planting_start_week:planting_end_week
                         
                         for (w in candidate_weeks) {
                           prev_week_precip <- precip_window %>% filter(week == w - 1) %>% pull(weekly_precip)
                           curr_week_precip <- precip_window %>% filter(week == w) %>% pull(weekly_precip)
                           
                           if (length(prev_week_precip) == 0) prev_week_precip <- 0
                           if (length(curr_week_precip) == 0) curr_week_precip <- 0
                           
                           two_week_cum <- prev_week_precip + curr_week_precip
                           prob_plant <- self$logistic_prob(two_week_cum, beta_0, beta_1, b_p)
                           
                           if (prob_plant > 0.5) {
                             plantable <- TRUE
                             plant_week <- w
                             break
                           }
                         }
                         
                         # Prepare input for ranger model
                         input <- data.frame(
                           Crop_Lag1 = factor(prev_y1_crop, levels = levels(crop_train$Crop_Lag1)),
                           Crop_Lag2 = factor(prev_y2_crop, levels = levels(crop_train$Crop_Lag2)),
                           Crop_Lag3 = factor(prev_y3_crop, levels = levels(crop_train$Crop_Lag3)),
                           CropSumry = factor(crop_summary, levels = levels(crop_train$CropSumry))
                         )
                         
                         prob <- predict(rf_model, data = input, type = "response")$predictions[1, ]
                         max_crop_prob <- max(prob, na.rm = TRUE)
                         
                         crop <- if (!plantable) "I" else names(which.max(prob))
                         
                         self$history[[as.character(year)]] <- list(
                           year = year,
                           crop = crop,
                           plantable = plantable,
                           plant_week = plant_week,
                           prob_plantable=prob_plant,
                           prob_max_rf = max_crop_prob
                         )
                       }
                     )
)
