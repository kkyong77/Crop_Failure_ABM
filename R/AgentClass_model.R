## including field scale soil moisture data 
Agent <- R6::R6Class("Agent",
                     public = list(
                       id = NULL,
                       insurance_coverage = NULL,
                       rotation_index = NULL,
                       lat = NULL,
                       lon = NULL,
                       acres = NULL,
                       geoid = NULL,
                       history = list(),
                       precip_ts = NULL,
                       moist_ts=NULL,
                       sm_ts=NULL,
                       
                       initialize = function(id, insurance_coverage, rotation_index, lat = NA, lon = NA, acres = NA,geoid=NA) {
                         self$id <- id
                         self$insurance_coverage <- insurance_coverage
                         self$rotation_index <- rotation_index
                         self$lat <- lat
                         self$lon <- lon
                         self$acres <- acres
                         self$geoid<-geoid
                       },
                       
                       add_history = function(year, crop, acres = NA) {
                         self$history[[as.character(year)]] <- list(
                           year = year,
                           crop = crop,
                           rf_crop=rf_crop,
                           acres=acres
                         )
                       },
                       
                       get_history_df = function() {
                         do.call(rbind, lapply(self$history, function(h) as.data.frame(h)))
                       },
                       
                       set_precip_ts = function(precip_df) {
                         self$precip_ts <- precip_df
                       },
                       
                       set_moist_ts = function(moist_df) {
                         self$moist_ts <-moist_df
                       },
                       
                       
                       set_sm_ts = function(sm_df) {
                         self$sm_ts <-sm_df
                       },
                       
                       
                       simulate_planting = function(year, crop, plant_week, plantable, prob_plantable, prob_max_rf = NA) {
                         self$history[[as.character(year)]] <- list(
                           year = year,
                           crop = crop,
                           rf_crop=rf_crop,
                           plant_week = plant_week,
                           plantable = plantable,
                           prob_plantable = prob_plantable,
                           prob_max_rf = prob_max_rf,
                           acres = self$acres,
                           county=self$geoid,
                           rci=self$rotation_index,
                           insurance=self$insurance_coverage
                         )
                       }
                     )
)
