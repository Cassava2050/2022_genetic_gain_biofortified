
# plant_ architecture traits ----------------------------------------------
interested_trials <- c(sel_data$studyName, genetic_gains$trial_name) %>% unique()


# Download data with QBMS -------------------------------------------------
library(QBMS)


# Conection with Cassava BreedBase server ---------------------------------
set_qbms_config("https://cassavabase.org/brapi/v1/calls/",
                path = "", time_out = 300, no_auth = TRUE,
                page_size = 10000,
                engine = "breedbase")



# Select a crop by name ---------------------------------------------------
set_crop("Cassava")

# Select the desire breeding program by name ------------------------------
set_program("CIAT")


# List all year’s trial in the selected program ---------------------------
trials <- list_trials() %>% pull()

trials <- trials[str_starts(trials, "CIAT")] 
trials <- trials[(35:44)] 
trials <- trials[-4]


# Downloading process -----------------------------------------------------
raw_data = list() # create a list to store the data

i = 1
for (trial in trials) {
  
  set_trial(trial)  
  # get observation variable ontology in the selected study/trial
  ontology <- get_trial_obs_ontology()
  # list all environments/locations information in the selected study/trial
  STUDIES <- list_studies()
  complete_studies <- STUDIES %>% 
    filter(str_detect(studyName, "BC"), 
           !str_detect(studyName, "ciat")) %>% pull(studyName)

  
  #for loop to extract the data of trials selected
  for(i in 1:length(complete_studies)) { 
    set_study(complete_studies[i])
    
    cat("\n_______________")
    cat("\nTRIAL:", complete_studies[i], "\n")
    cat("_______________\n")
    
    raw_data[[complete_studies[i] ]] = get_study_data()
  }
}


# Convert list into a data.frame ------------------------------------------
all_raw = data.table::rbindlist(raw_data, fill = TRUE) %>% 
  as_tibble() %>%
  filter(observationLevel == "plot") 


# save raw data -----------------------------------------------------------
folder_output <- here::here("data//")
meta_file_name <- paste0(folder_output, paste0("2023_", "_pheno_PT_", 
                                               Sys.Date(), ".csv", sep = ""))
write.csv(all_raw, file = meta_file_name, row.names = F, na = "")




