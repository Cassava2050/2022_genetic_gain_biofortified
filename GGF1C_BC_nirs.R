rm(list = ls())

source("https://raw.githubusercontent.com/Cassava2050/PPD/main/utilities_tidy.R")
source("GG_functions.R")

## Load the files to check

# Loading BC content ------------------------------------------------------

  trait = "BC"
  
  local_file <- "yes" #
  
  if (local_file == "yes") {
    folder <- here::here("data//")  
    file <- "phenotype_bc_1.csv"
    skip_col <- 3 # double check the number of col skipped
    trial_interest = "GGBC"
    year_interest <- 2023
  }
  
  # 1) load the data
  sel_data <- read_cassavabase(phenotypeFile = paste0(folder, file))
  
  length(unique(sel_data$studyName)) # 27 trials
  
  sel_data$studyYear %>% unique()
  
  # ---- Change columns into standar names ----
  sel_data_kp <- change_colname(sel_data, NA)
  
  
  ## change the column class
  
  obs_col <- c(
    names(sel_data_kp)[str_detect(names(sel_data_kp), "obs_")],
    "use_rep_number", "blockNumber",
    "use_plot_number", "use_plot_width",
    "use_plot_length"
  )
  sel_data_kp %<>%
    mutate(across(all_of(obs_col), as.numeric))
  
  
  # remove - , replace by _
  names(sel_data_kp) = gsub("-", "_", names(sel_data_kp))
  
  ## Duplications in row and cols
  duplicated_plot <- row_col_dup(sel_data_kp) # non applicable to non row-col designs
  
  ## Plot trial layout
  #trial_layout(sel_data_kp) 
  
  sel_data_kp <- sel_data_kp %>% 
    mutate(use_accession_name = recode_factor(use_accession_name, 
                                              COSTENA = "Costena",
                                              VENEZOLANA = "Venezolana")) 
  
  ## Check the clone name
  cloneName_new_old <- check_clone_name(
    clone_list = sel_data_kp$use_accession_name,
    new_names = NA,
    add_check = NULL
  )
  
  trial_standard <- sel_data_kp %>%
    left_join(cloneName_new_old,
              by = c("use_accession_name" = "accession_name_ori")
    ) %>%
    select(-use_accession_name) %>%
    rename(use_accession_name = use_accession_name.y)
  
# Is numeric all trait data?
  is_numeric(trial_data = trial_standard)
  
# remove missing data  
  trial_standard <- trial_standard %>% select(use_trial_name, use_accession_name, obs_betacarotenoid_nirs) %>% 
    drop_na()
  
# Standardize check clone names
  trial_standard <- trial_standard %>% mutate(use_accession_name = recode_factor(use_accession_name,
                                                               `CM4919-1_is_Veronica` = "CM4919-1",
                                                               `COL2215_is_Venezolana` = "COL2215",
                                                               `SMB2446-2_is_Caiseli` = "SMB2446-2",
                                                               `TAI8_is_TAI` = "TAI8"))
  

  
  ## Get the tidy data
  
  meta_info = names(trial_standard)[str_detect(names(trial_standard), "use_")]
  meta_info = gsub("use_", "", meta_info)
  meta_info
  trial_tidy = trial_standard
  names(trial_tidy)= gsub("use_", "", names(trial_standard))
  # observations
  trait_list = names(trial_tidy)[str_detect(names(trial_tidy), "obs_")]
  trait_list = gsub("obs_", "", trait_list)
  trait_list
  names(trial_tidy)= gsub("obs_", "", names(trial_tidy))
  trial_tidy = trial_tidy[c(meta_info, trait_list)]
  
# let's summarize betacarotenoids   
  
  
  

