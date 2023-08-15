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
                                                               `CG1141-1_is_Costena` = "CG1141-1",
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
  
# let's summarize betacarotenoids by NIRS
  trial_tidy <- trial_tidy %>% group_by(accession_name) %>% 
    summarise(betacarotenoid_nirs = mean(betacarotenoid_nirs, na.rm = T))

# load old data from Xiaofei
  old_BC_data <- read.csv("data/BC_AYT_BC_old years_2021-05-10.csv")
  colnames(old_BC_data)[2] <- c("betacarotenoid_nirs")
  
# bind both data sets
  pvals1 <- trial_tidy %>% bind_rows(old_BC_data) %>% 
    group_by(accession_name) %>% 
    summarise(betacarotenoid_nirs = mean(betacarotenoid_nirs, na.rm = T))

# Merging Predicted value with the yearOrigin
  origin <- read_csv("data/crossing_year_update_2023-08-14.csv") 
  origin <- origin[!origin %>% duplicated(), ]
  
    # cross year missing values per accession_name
pvals1 %>% left_join(origin, by = "accession_name") %>% 
    filter(is.na(year)) %>% distinct(accession_name) %>% 
    write.table("clipboard", col.names = T, row.names = F, sep = "\t", na = "")
  
# cross year missing values per family
  pvals1 %>% left_join(origin, by = "accession_name") %>% 
    filter(str_detect(accession_name, "-"), is.na(cross_year)) %>% 
    separate(accession_name, c("family", "offspring_code"), sep = "-") %>% 
    distinct(family) %>% write.table("clipboard", col.names = T, row.names = F, sep = "\t")  

# merge betacarotenes pvals with crossing year data 
stack <- merge(pvals1, origin, by = "accession_name")

# linear model
library(ggpubr)

t1 <- stack %>% 
  filter(year >= 2007) %>% 
  ggplot(aes(x = year, y = betacarotenoid_nirs )) +
  geom_point()+
  stat_regline_equation()+
  # geom_abline(intercept = -320.462, slope = 0.1639, color="blue",
  #             size=1) +
  labs(x = "crossing year", y = paste("Predicted", trait, sep = "_")) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE) +
  theme_xiaofei()    
t1

ggsave(paste("images/", "GG", trait, sep = "_", ".png"),
       plot = t1, units = "in", dpi = 300, width = 6, height = 6)

model <- lm(formula = betacarotenoid_nirs ~ year, data = stack %>% 
             filter(year >= 2007))

# -------------------------------------------------------------------------
# lineal model
gg_model <- agriutilities::parameters_gg(model = model, trait = trait)
gg_model %>% mutate(across(where(is.numeric), round, 4)) %>% 
  t() %>% as.data.frame() %>% 
  write.table("clipboard", col.names = F, row.names = T, sep = "\t")

# manual way
intercept = model$coefficients[1]
slope = model$coefficients[2]
lastYear = (2014*slope)  + intercept
firstYear = (2007*slope)  + intercept

slope/firstYear * 100  

  
 

