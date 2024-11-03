#########################################################################################
# Replication Code for GOV 2020 class excersises
# Author: Kentaro Nakamura
# Last Update: 2024/10/27
#########################################################################################
rm(list = ls())

#set wordking directory
getwd()
dir_path <- "/n/holyscratch01/imai_lab/Users/knakamura/GeoTwitter" #change to your directory
setwd(dir_path)
getwd()

## Install & load packages 
detachAllPackages <- function() {
  basic.packages <- c("package:stats","package:graphics","package:grDevices","package:utils","package:datasets","package:methods","package:base")
  package.list <- search()[ifelse(unlist(gregexpr("package:",search()))==1,TRUE,FALSE)]
  package.list <- setdiff(package.list,basic.packages)
  if (length(package.list)>0)  for (package in package.list) detach(package, character.only=TRUE)
}
detachAllPackages()
list.of.packages <- c("tidyverse", "haven", "lme4", "CausalImpact", "patchwork", "lubridate", "bsts", "sf", "doParallel", "foreach")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]; if(length(new.packages)){install.packages(new.packages,dependencies=TRUE)}
loaded.packages <- lapply(list.of.packages, require, character.only = TRUE)
rm(list.of.packages,new.packages,loaded.packages,detachAllPackages)

## loading data
data <- read_dta("simil_dgm_final.dta") #required for STEP1 and STEP2(a)
constituencies <- st_read("btw17_geometrie_wahlkreise_vg250_shp/Geometrie_Wahlkreise_19DBT_VG250.shp") #constituency-level map (needed for visualization)

##################################################################################
# STEP 1: Replication of Original Article
##################################################################################
## <---------------- Data Preprocessing ----------------> ##
head(data)
## A tibble: 6 × 11
#simil const party period   day month  year partyid edate      state_number state
#<dbl> <dbl> <chr>  <dbl> <dbl> <dbl> <dbl>   <dbl> <date>            <dbl> <dbl>
#1  0.226      3 afd     3405     1     5  2017       1 2017-05-01            1     1
#2  0.356      3 afd     3315    31     1  2017       1 2017-01-31            1     1
#3  0.428      3 afd     3185    23     9  2016       1 2016-09-23            1     1
#4  0.308      3 afd     2837    11    10  2015       1 2015-10-11            1     1
#5 -0.0274     3 afd     3363    20     3  2017       1 2017-03-20            1     1
#6 -0.0259     3 afd     2886    29    11  2015       1 2015-11-29            1     1


# Create 'event' column
data <- data %>%
  dplyr::mutate(
    # Events defined in the paper (that will be changed later)
    event = case_when(
      day == 13 & month == 11 & year == 2015 ~ 1,
      day == 1 & month == 1 & year == 2016 ~ 1,
      day == 22 & month == 3 & year == 2016 ~ 1,
      day == 14 & month == 7 & year == 2016 ~ 1,
      day == 19 & month == 12 & year == 2016 ~ 1,
      day == 22 & month == 3 & year == 2017 ~ 1, 
      day == 20 & month == 4 & year == 2017 ~ 1, 
      day == 22 & month == 5 & year == 2017 ~ 1,
      day == 3 & month == 6 & year == 2017 ~ 1,
      day == 16 & month == 8 & year == 2017 ~ 1,
      day == 15 & month == 9 & year == 2017 ~ 1,
      TRUE ~ 0),
    # Dummy for 2016
    year_2016 = ifelse(year == 2017, 1, 0)) %>%
  # Create 'period' variable
  dplyr::arrange(year, month, day) %>%
  dplyr::mutate(period = group_indices(., year, month, day)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    time = period - 1,
    time2 = time^2
  )

# For each events...
for (i in 1:11) {
  # Create transition variables
  data <- data %>%
    mutate(!!paste0("trans", i) := ifelse(period >= c(71, 120, 201, 315, 473, 566, 595, 627, 639, 713, 743)[i] &
                                            period < c(120, 201, 315, 473, 566, 595, 627, 639, 713, 743, Inf)[i], 1, 0))
  # Create recovery variables
  data <- data %>%
    mutate(rec = ifelse(period > c(71, 120, 201, 315, 473, 566, 595, 627, 639, 713, 743)[i] &
                          period < c(120, 201, 315, 473, 566, 595, 627, 639, 713, 743, Inf)[i], 
                        as.integer(factor(period)) - c(71, 120, 201, 315, 473, 566, 595, 627, 639, 713, 743)[i], 0)) %>%
    group_by(period) %>%
    mutate(!!paste0("recov", i) := max(rec)) %>%
    ungroup() %>%
    select(-rec)
}
# Create recov^2 for all events
data <- data %>% dplyr::mutate(mutate(across(starts_with("recov"), ~ .x^2, .names = "{.col}_2"),))

## <---------------- Data Preprocessing (END) ----------------> ##

## <---------------- Discontinuous Growth Model ----------------> ##
# For DGM, read Bliese and Lang (2016). Here, the quantity of interest is not the immediate impact on the discourse (what they call "absolute decline");
# the authors are interested in the relative decline. From Bliese and Lang (2016),
#
#     "That is, it is fundamentally different to expect an intervention designed to curb losses in customer satisfaction to 
#     (a) immediately increase satisfaction and produce a positive trajectory versus (b) hold the decline steady relative to 
#     the pre-intervention and produce a subsequent loss trajectory that is less pronounced than the pre-intervention slope (but still a decline).
#
# In the following, I will replicate the original paper's result and DGM w/ single time-trend.

## Fit discontinuous growth model (Caution: it is computationally expensive)
model <- lme4::lmer(simil ~ time + time2 + trans1 + trans2 + trans3 + trans4 + trans5 + 
                trans6 + trans7 + trans8 + trans9 + trans10 + trans11 +
                recov1 + recov2 + recov3 + recov4 + recov5 + recov6 + recov7 + recov8 + 
                recov9 + recov10 + recov11 + 
                recov1_2 + recov2_2 + recov3_2 + recov4_2 + 
                recov5_2 + recov6_2 + recov7_2 + recov8_2 + recov9_2 + recov10_2 + recov11_2 + 
                year_2016 +
                (1 + time + trans1 + trans2 + trans3 + trans4 + trans5 + trans6 + 
                   trans7 + trans8 + trans9 + trans10 + trans11 + recov1 + recov2 + recov3 + 
                   recov4 + recov5 + recov6 + recov7 + recov8 + recov9 + recov10 + recov11 || const), #constituency level random effect
              data = subset(data, party == "afd"),
              control = lme4::lmerControl(optimizer = "bobyqa")
              )
## Check the outputs (the result in the original paper is in Appendix Table F.2 p.22)
summary(model, correlation = FALSE)

coefficients <- lme4::fixef(model)[4:14]
conf_intervals <- confint(model, level = 0.95, method = "Wald")
as.data.frame(conf_intervals) %>%
  tibble::rownames_to_column(var = "Term") %>%
  dplyr::filter(Term %in% names(coefficients)) %>%
  dplyr::mutate(
    Coefficient = coefficients[Term],
    Term = factor(Term, 
                  levels = c("trans1", "trans2", "trans3", "trans4", "trans5", "trans6", "trans7", "trans8", "trans9", "trans10", "trans11"))
    ) %>%
  dplyr::rename("Lower" = "2.5 %", "Upper" = "97.5 %") %>%
  ggplot2::ggplot(aes(x = Term, y = Coefficient)) +
  ggplot2::geom_point() +
  ggplot2::geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2) +
  ggplot2::theme_minimal() +
  ggplot2::labs(title = "Coefficient Plot with 95% Confidence Intervals", x = "Terms", y = "Coefficients")

# However, you will see that once we remove the quadratic trend, then the discovered effects are reversed.
# This is because the inference is based on the model-based inference, and the model is terrible.
## model w/o quadratic trends
model2 <- lme4::lmer(simil ~ time + time2 + trans1 + trans2 + trans3 + trans4 + trans5 + 
                      trans6 + trans7 + trans8 + trans9 + trans10 + trans11 +
                      recov1 + recov2 + recov3 + recov4 + recov5 + recov6 + recov7 + recov8 + 
                      recov9 + recov10 + recov11 + 
                      year_2016 +
                      (time + trans1 + trans2 + trans3 + trans4 + trans5 + trans6 + 
                         trans7 + trans8 + trans9 + trans10 + trans11 + recov1 + recov2 + recov3 + 
                         recov4 + recov5 + recov6 + recov7 + recov8 + recov9 + recov10 + recov11 || const), #constituency level random effect
                    data = subset(data, party == "afd"),
                    control = lme4::lmerControl(optimizer = "bobyqa")
)

summary(model2, correlation = FALSE)
coefficients2 <- fixef(model2)[4:14]
conf_intervals2 <- confint(model2, level = 0.95, method = "Wald")

conf_intervals2 %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "Term") %>%
  dplyr::filter(Term %in% names(coefficients)) %>%
  dplyr::mutate(
    Coefficient = coefficients2[Term],
    Term = factor(Term, levels = c("trans1", "trans2", "trans3", "trans4", "trans5", "trans6", "trans7", "trans8", "trans9", "trans10", "trans11"))) %>%
  dplyr::rename("Lower" = "2.5 %", "Upper" = "97.5 %") %>%
  ggplot2::ggplot(aes(x = Term, y = Coefficient)) +
  ggplot2::geom_point() +
  ggplot2::geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2) +
  ggplot2::theme_minimal() +
  ggplot2::labs(title = "Coefficient Plot with 95% Confidence Intervals", x = "Terms", y = "Coefficients")

# Is the model accurate? Does it give us the reasonable prediction? One way to compare is 
# to check how well the model predicts in-sample predictions
pred <- predict(model)
event_dates <- as.Date(c("2015-11-13", "2016-01-01", "2016-03-22", "2016-07-14", 
                         "2016-12-19", "2017-03-22", "2017-04-20", "2017-05-22", 
                         "2017-06-03", "2017-08-16", "2017-09-15"))
data.frame(real = data$simil, pred, date = as.Date(data$edate), const = data$const) %>%
  dplyr::group_by(date) %>%
  dplyr::summarise(
    sum_real = sum(real),
    sum_pred_quadratic = sum(pred)
  ) %>%
  tidyr::pivot_longer(!c(date)) %>%
  dplyr::mutate(name = dplyr::case_when(
    name == "sum_pred_linear" ~ "DGM w/ Linear Term",
    name == "sum_pred_quadratic" ~ "DGM w/ Quadratic Term (Original)",
    name == "sum_real" ~ "Real Similarity")
  ) %>%
  ggplot2::ggplot(aes(x = date, y = value, color = name)) +
  ggplot2::geom_line() +
  ggplot2::geom_vline(xintercept = event_dates, linetype = "dashed", color = "black") + # Add vertical lines
  ggplot2::theme_classic() +
  ggplot2::ylab("Similarity") +
  ggplot2::xlab("Dates") +
  ggplot2::theme(legend.position = "bottom")
# So, even though the model seems fitting better if we include quadratic term, it seems not capturing time trend well...
## <---------------- Discontinuous Growth Model (End) ----------------> ##

##################################################################################
# STEP 2: Extensions
##################################################################################
### <---------------- Processing Data ----------------> ###
## Here are the code to select the events. This data is from Global Terrorism Database
EU_terrorism <- read.csv("EU_terrorism.csv") %>%
  dplyr::mutate(date = ymd(substr(eventid, 1, 8)) %>% as.Date()) %>%
  dplyr::filter( (date > as.Date("2015-09-03")) & (date < as.Date("2017-09-24"))) %>%
  dplyr::filter(country_txt != "Germany") %>% #to prevent endogeneity
  dplyr::filter(gname %in% c("Jihadi-inspired extremists", "Muslim extremists", "Islamic State of Iraq and the Levant (ISIL)", "Kurdistan Workers' Party (PKK)")) %>%
  dplyr::filter(nkill > 0)

table(EU_terrorism$date)
#> 2015-11-13 2016-01-07 2016-03-22 2016-06-13 2016-07-14 2016-07-26 2016-09-01 2016-12-23 2017-03-18 
#> 8          1          2          1          1          1          1          1          1 
#> 2017-04-07 2017-04-20 2017-06-19 2017-06-20 2017-06-30 2017-08-17 2017-08-18 2017-08-25 
#> 1          1          1          1          1          2          2          1 

#Note that 2017-06-19 and 2017-06-20 / 2017-08-17 and 2017-08-18, we only use 2017-06-19 / 2017-08-17 as a cutoff.
new_event_data <- data.frame(
  date = as.Date(c("2015-11-13", "2016-01-07", "2016-03-22", "2016-06-13", "2016-07-14", 
                   "2016-07-26", "2016-09-01", "2016-12-23", "2017-03-18", "2017-04-07", 
                   "2017-04-20", "2017-06-20", "2017-06-30", "2017-08-17", "2017-08-25"))
)

# to check the period of each event
event_period <- data %>%
  dplyr::mutate(date = make_date(year, month, day)) %>%
  dplyr::right_join(new_event_data) %>%
  dplyr::select(date, period) %>%
  unique()

# A tibble: 15 × 2
#date       period
#<date>      <int>
#1 2015-11-13     71
#2 2016-01-07    126
#3 2016-03-22    201
#4 2016-06-13    284
#5 2016-07-14    315
#6 2016-07-26    327
#7 2016-09-01    364
#8 2016-12-23    477
#9 2017-03-18    562
#10 2017-04-07    582
#11 2017-04-20    595
#12 2017-06-20    656
#13 2017-06-30    666
#14 2017-08-17    714
#15 2017-08-25    722

#########################################
# Part (a): Model Validation -> Please take a look at Zhiyu's code (if any)
#########################################

#########################################
# Part (b): Estimation of Causal Effects with Better Models
#########################################
### <---------------- Estimation of Causal Effect with State Space Models ----------------> ###
run_causal_impact <- function(event_name, event_date, data, day_window = 20) {
  pre_event_period <- c(event_date - day_window, event_date - 1)  # (window) days before
  post_event_period <- c(event_date, event_date + day_window)     # (window) days after
  
  time_series <- data %>%
    dplyr::mutate(date = as.Date(edate)) %>%
    dplyr::filter(date >= pre_event_period[1] & date <= post_event_period[2]) %>%
    dplyr::arrange(date) %>%
    dplyr::select(simil)
  
  # Run CausalImpact
  pre_event_int <- length(seq(pre_event_period[1], pre_event_period[2], by = "day"))
  post_event_int <- length(seq(pre_event_period[1], post_event_period[2], by = "day"))
  
  if (nrow(time_series) > pre_event_int){
    impact <- CausalImpact::CausalImpact(data = time_series, c(1, pre_event_int), c(pre_event_int + 1, post_event_int), model.args = list(nseasons = 2))
  } else {
    impact <- NA
  }
  return(impact)
}

event_dates <- list(
  "Paris" = as.Date("2015-11-13"), 
  "Paris2" = as.Date("2016-01-07"),
  "Brussels airport" = as.Date("2016-03-22"), 
  "Magnanville" = as.Date("2016-06-14"),
  "Nice" = as.Date("2016-07-14"),
  "Saint-Etienne-du-Rouvray" = as.Date("2016-07-26"),
  "Copenhagen" = as.Date("2016-09-01"),
  "Milan" = as.Date("2016-12-23"), 
  "Paris3" = as.Date("2017-03-18"),
  "Stockholm" = as.Date("2017-04-07"),
  "Paris4, Champs-Elysees" = as.Date("2017-04-20"), 
  "Paris5, Brussels" = as.Date("2017-06-19"), 
  "Linz" = as.Date("2017-06-30"), 
  "Barcelona / Turku" = as.Date("2017-08-16"),
  "Barcelona2" = as.Date("2017-08-25")
)

#Caution: the following code is computationally expensive.
causal_impact_results <- list()
for (event in names(event_dates)) {
    for (const in names(table(data$const))){
      cat("Estimating Event", event, " at Constituency No.", const, "\n")
      numeric_const <- as.numeric(const)
      causal_impact_results[[event]][[const]] <- run_causal_impact(event, event_dates[[event]], subset(data, party == "afd" & const == numeric_const))
    }
}

causal_impact_results_df <- do.call(rbind, lapply(names(causal_impact_results), function(event) {
  do.call(rbind, lapply(names(causal_impact_results[[event]]), function(constituency) {
    entry <- causal_impact_results[[event]][[constituency]]
    
    # Check if entry is a list and has 'summary' component to avoid atomic vectors
    if (is.list(entry) && !is.null(entry$summary)) {
      att <- entry$summary$AbsEffect[1]
      se <- entry$summary$AbsEffect.sd[1]
      
      # Only include non-NA att and se values
      if (!is.na(att) && !is.na(se)) {
        data.frame(
          event = event,
          constituency = as.integer(constituency),
          att = att,
          se = se
        )
      }
    }
  }))
}))

merged_result <- merge(causal_impact_results_df, constituencies, 
                     by.x = "constituency", by.y = "WKR_NR",
                     all.x = TRUE)

#visualization
merged_result %>%
  ggplot() +
  geom_sf(aes(fill = att, geometry = geometry)) +  # Assuming `geometry` column has sf geometries
  scale_fill_viridis_c(option = "plasma", name = "ATT") +  # Choose a color scale for att
  facet_wrap(~ event, ncol = 5) +  # Create a map for each event
  theme_minimal() +
  theme(axis.text = element_blank(), axis.ticks = element_blank())

#save.image(file = "all_data.RData")
