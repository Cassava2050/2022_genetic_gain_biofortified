rm(list = ls())

source("https://raw.githubusercontent.com/Cassava2050/PPD/main/utilities_tidy.R")
source("GG_functions.R")

## Load the files to check

local_file <- "yes" #

if (local_file == "yes") {
  folder <- here::here("data//")  
  file <- "phenotype.csv"
  skip_col <- 3 # double check the number of col skipped
  trial_interest = "GGBC"
  year_interest <- 2023
}

# 1) load the data
sel_data <- read_cassavabase(phenotypeFile = paste0(folder, file))

length(unique(sel_data$studyName)) # 27 trials


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
trial_layout(sel_data_kp) 

sel_data_kp <- sel_data_kp %>% 
  mutate(use_accession_name = recode_factor(use_accession_name, COSTENA = "Costena")) 

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

## Add GIS data
trial_standard <- add_GIS(trial_standard)

## Location Map

#![](images/map.png)


accession_rep_ct <- trial_standard %>%
  count(use_trial_name, use_accession_name, use_rep_number)  %>%
  arrange(use_trial_name) %>%
  filter(n>1)
accession_rep_ct 


## Genotypes per trial


conducted_trials <- 
  trial_standard %>% group_by(use_trial_name, use_plant_date,use_harvest_date, use_location) %>% 
  summarise(n_gen = n_distinct(use_accession_name)) %>% 
  mutate(harvesting_time = 
           interval(ymd(use_plant_date), ymd(use_harvest_date)) %>% as.period,
         harvesting_time = paste0(harvesting_time@month, "month ", harvesting_time@day, "day")) %>% 
  ungroup()

conducted_trials %>% View() # there are some trials with non harvest day.
# double check them when I have kept with the valid trials

# plot plant number
plants_plot <- trial_standard %>%
  group_by(use_trial_name) %>%
  count(obs_planted_number_plot) 

plants_plot

## Frequency harvest plant number

plants_harvested <- trial_standard %>%
  group_by(use_trial_name) %>%
  count(obs_harvest_number) %>% arrange(desc(obs_harvest_number))

plants_harvested %>% left_join(plants_plot %>% select(-n), by = "use_trial_name") %>% 
  left_join(conducted_trials %>% select(use_trial_name, n_gen), by = "use_trial_name") %>% 
  arrange(desc(n))

# Is numeric all trait data?
is_numeric(trial_data = trial_standard)

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

# Select the interest traits
trial_tidy <- trial_tidy %>% select(all_of(meta_info), DM_gravity, yield_ha)

# Now I need to remove the trials with no yield data
trait = "yield_ha"

exp <- trial_tidy %>% 
  as_tibble() %>% 
  group_by(trial_name) %>% 
  summarise(n_miss = sum(is.na(!!sym(trait))),
            n  =  n(), 
            percen = round(n_miss/n, 3),
            reps = n_distinct(rep_number)) %>% 
  arrange(desc(percen)) %>% 
  filter(percen != 1 ) %>% 
  pull(trial_name) 

# remove trials with 100% of missing data

trial_tidy <- trial_tidy %>% 
  filter(trial_name %in% exp)

# Boxplots
my_dat_noNA <- trial_tidy[, colSums(is.na(trial_tidy)) < nrow(trial_tidy)]
trait_wanted <- names(my_dat_noNA)[names(my_dat_noNA) %in% trait_list]
for (i in 1:length(trait_wanted)) {
  y_DATA <- my_dat_noNA[[trait_wanted[i]]]
  x_DATA <- my_dat_noNA$trial_name
  my_DATA <- my_dat_noNA
  y_LABEL <- trait_wanted[i]
  x_LABEL <- NULL
  TITLE <- NULL
  y_MAX <- max(y_DATA, na.rm = TRUE) * 1.2
  y_MIN <- 0
  plot_box <- ggplot(my_DATA, aes(x = x_DATA, y = y_DATA)) +
    geom_violin(trim = FALSE, fill = "gray") +
    geom_boxplot(width = 0.2) +
    coord_cartesian(ylim = c(y_MIN, y_MAX)) +
    theme_xiaofei() +
    labs(
      y = y_LABEL, x = x_LABEL,
      title = TITLE
    )
  plot(plot_box)
}

# Load xiaofei's data
genetic_gains <- read.csv("data/phenotype_old_years.csv") %>% 
  select(-cross_year)

# remove obs_ string from header
names(genetic_gains) = gsub("obs_", "", names(genetic_gains))

# pulling trials with no data
exp_old_ata <- genetic_gains %>% 
  group_by(trial_name) %>% 
  summarise(n_miss = sum(is.na(!!sym(trait))),
            n  =  n(), 
            percen = round(n_miss/n, 3),
            reps = n_distinct(rep_number)) %>% 
  arrange(percen) %>% 
  filter(percen != 1 ) %>% 
  pull(trial_name) 

# Keeping trials with yield data
genetic_gains <- genetic_gains %>% filter(trial_name %in% exp_old_ata)

# bind both data frame
trial_tidy_new <- trial_tidy %>% bind_rows(genetic_gains) %>% as_tibble() 
str(trial_tidy_new)

# becoming factor some meta data cols
trial_tidy_new$year = as.factor(trial_tidy_new$year)
trial_tidy_new$rep_number = as.factor(trial_tidy_new$rep_number)


# plot the boxplot across years
trial_tidy_new %>% 
  droplevels() %>% 
  ggplot(aes(x = trial_name, y = !!sym(trait), fill = year)) + 
  geom_boxplot(width = 0.4) +
  coord_cartesian(ylim = c(y_MIN, y_MAX)) +
  theme_xiaofei() 
ggsave("images/yield.png", units = "in", dpi = 300, width = 20, height = 7)

# Number of shared information
shared <- check_connectivity(
  data = trial_tidy_new,
  genotype = "accession_name",
  trial = "trial_name",
  response = NULL,
  return_matrix = TRUE
) 

# shared information matrix
shared_cor(shared, size = 2) +
  theme_xiaofei() +
  labs(x = NULL, y = NULL) +
  theme(legend.position = "none", panel.background = ggplot2::element_blank()) 
ggsave("images/shared_info.png", units = "in", dpi = 300, width = 15, height = 8)


exp <- unique(c(exp, exp_old_ata)) 
# Models
datos <- trial_tidy_new 

# Equation yield_ha ~ rep_number + (1 | accession_name)
equation <- reformulate(c("rep_number", ran("accession_name") ), response = "yield_ha")

# running mixed model
objt <- mult_lme4(data = genetic_gains, equation = equation, var_sub = "trial_name")

# After this step, 12 trial models did not converge and were removed from the analysis
setdiff(exp, trial_converged) %>% as.data.frame() %>% 
  write.table("clipboard", col.names = T, row.names = F, sep = "\t")

# extracting Vg, Ve, H2 & outliers
mt_summ <- mult_summary(objt , gen = "accession_name", y = "yield_ha")
summ_height0 <- mt_summ

# extracting outliers
dataOut <- lapply(objt, lme4_res, T)
dataOut <- data.frame(plyr::ldply(dataOut[], data.frame, .id = "Experiment"))

# trials with outliers
dataOut %>% filter(Classify == "Outlier") %>% pull(Experiment) -> trial_outliers
dataOut %>% filter(Classify == "Outlier") %>% 
  mutate(across(where(is.numeric), round, 2))

# trials with the most extreme outliers
trial_tidy_new %>% filter(trial_name %in% trial_outliers) %>% 
  droplevels() %>% 
  ggplot(aes(x = trial_name, y = !!sym(trait), fill = year)) + 
  geom_boxplot(width = 0.4) +
  labs(x = NULL) +
  theme_xiaofei() 
ggsave("images/outliers_trial.png", units = "in", dpi = 300, width = 12, height = 6)

# putting outliers as NA
dataOut[dataOut$Classify=="Outlier", "yield_ha"  ] <- NA

# run again without outliers

objt <- mult_lme4(data = dataOut, equation = equation, var_sub = "Experiment")
mt_summ2 <- mult_summary(models = objt ,gen = "accession_name", y = "yield_ha")
summ_height <- mt_summ2
print(mt_summ2)

# ploting mt_summ2
mt_summ2 %>% select() %>% pivot_longer()

# there is a trial with H2 = 0 2019103BCPRC_cbia
trial_tidy_new %>% filter(trial_name == "2019103BCPRC_cbia") %>%
  select(plot_name, rep_number, accession_name, yield_ha) %>%
  ggplot(aes(x = reorder(accession_name, yield_ha), y = yield_ha)) +
  geom_boxplot() +
  geom_point(aes(color = rep_number)) +
  labs(x = NULL, title = "2019103BCPRC_cbia") +
  theme_xiaofei()
ggsave("images/heri_0.png", units = "in", dpi = 300, width = 12, height = 7)

# trials kept
exp_valids <- mt_summ2 %>% filter(h2 >=  0.1 ) %>% pull(Experiment)

# phenotypic yield data across trials kept
trial_tidy_new %>% filter(trial_name %in% exp_valids) %>% 
  droplevels() %>% 
  ggplot(aes(x = trial_name, y = !!sym(trait), fill = year)) + 
  geom_boxplot(width = 0.4) +
  labs(x = NULL) +
  theme_xiaofei()
ggsave("images/kept_trials_yield.png", units = "in", dpi = 300, width = 12, height = 7)
  
# Save the tidy data ------------------------------------------------------


cleanDT <- dataOut %>% filter(Experiment %in% exp_valids) %>% select(-res, -Classify)

variables <- trial_tidy_new %>% select(year , trial_name, plot_name, location, accession_name, rep_number )
variables$rep_number <- as.numeric(variables$rep_number)

cleanDT <- 
  merge(cleanDT,
        variables, 
        by.x = c("Experiment", "rep_number", "accession_name"),
        by.y = c("trial_name", "rep_number", "accession_name"), all.x = TRUE)

# Save the tidy data
folder_output <- here::here("output//")
meta_file_name <- paste0(folder_output, paste0("2023_", trial_interest, "_tidy_data4analysis_", trait, Sys.Date(), ".csv", sep = ""))
write.csv(cleanDT, file = meta_file_name, row.names = F, na = "")


## Load the tidy data
trial_set_number = 1

# all files in the folder
list_file = list.files(here::here("output"))

# tidy data of the trials interested
sel_file = list_file[str_detect(list_file, "_tidy_data4analysis_") &
                       str_detect(list_file, trait)&
                       str_detect(list_file, paste(year_interest, trial_interest, sep = "_"))]
# the data we will use
sel_file_use = sel_file[1]

sel_file_use
trial1_tidy = read.csv(here::here("output", sel_file_use), header=TRUE,
                       stringsAsFactors = FALSE,
                       as.is=T,
                       check.names = FALSE)
if(trial_set_number == 1){
  trial_tidy_all = trial1_tidy
}

cleanDT <- trial_tidy_all
str(cleanDT)

# becoming factor some meta data cols
cleanDT$year = as.factor(cleanDT$year)
cleanDT$rep_number = as.factor(cleanDT$rep_number)
cleanDT$Experiment = as.factor(cleanDT$Experiment)
cleanDT$location = as.factor(cleanDT$location)
cleanDT$accession_name = as.factor(cleanDT$accession_name)


# Genetic gain calculation Yield

equation_fixed <-  reformulate(c("1", "accession_name", "year"), response = trait)

library(asreml)

# model 1
model_1 <- asreml(fixed = equation_fixed,
                random = ~ Experiment + location + 
                  
                  # interactions between the respective factors
                  year:location + accession_name:year + accession_name:location + Experiment:rep_number,
                
                # units represents the individual experimental units or plots within the experiment.
                residual = ~units, data = cleanDT) 
aic_1 <- summary(model_1)$aic
bic_1 <- summary(model_1)$bic

# model 2
model_2 <- asreml(fixed = equation_fixed,
                random = ~ Experiment + location + year:location + accession_name:year + 
                  accession_name:location + Experiment:rep_number,
                
                # residuals (units) are directly modeled separately for each level of the factor Experiment.
                residual = ~ dsum(~ units | Experiment), data = cleanDT)
aic_2 <- summary(model_2)$aic
bic_2 <- summary(model_2)$bic

# model 3
model_3 <- asreml(fixed = equation_fixed,
                random = ~ Experiment + location + year:location + accession_name:year + 
                  accession_name:location + at(Experiment):rep_number,
                residual = ~ dsum(~ units | Experiment), data = cleanDT)
model_3 <- update.asreml(model_3)
model_3 <- update.asreml(model_3)
aic_3 <- summary(model_3)$aic
bic_3 <- summary(model_3)$bic

# model 4
model_4 <- asreml(fixed = equation_fixed,
                random = ~ Experiment + Experiment:accession_name + 
                  at(Experiment):rep_number,
                residual = ~ dsum(~ units | Experiment), data = cleanDT)
model_4 <- update.asreml(model_4)
aic_4 <- summary(model_4)$aic
bic_4 <- summary(model_4)$bic

# model 5
model_5 <- asreml(fixed = equation_fixed,
                random = ~ Experiment + diag(Experiment):accession_name + at(Experiment):rep_number,
                residual = ~ dsum(~ units | Experiment), data = cleanDT)
model_5 <- update.asreml(model_5)
aic_5 <- summary(model_5)$aic
bic_5 <- summary(model_5)$bic

# Model parameters quality
tibble(
  Model = c("Model 1", "Model 2", "Model 3", "Model 4", "Model 5"),
  AIC = c(
    aic_1[1],
    aic_2[1],
    aic_3[1],
    aic_4[1],
    aic_5[1]
  ),
  BIC = c(
    bic_1[1],
    bic_2[1],
    bic_3[1],
    bic_4[1],
    bic_5[1]
  )
) %>% arrange(AIC)

# function for computing the variance-covariance
source("extractG.R")


extractG(model_5, gen = "accession_name", env = "Experiment", vc.model = "diag")$VCOV %>%
  diag() %>%
  data.frame(Exp = names(.), VarG = ., row.names = NULL) %>%
  ggplot(aes(x = Exp, y = VarG))+
  geom_bar(stat = "identity")+
  theme(axis.text.x = element_text(hjust = 1 , angle = 75))

# Merging Predicted value with the yearOrigin

model = ""

pvals1 <- predict(model_5, classify = "accession_name")$pvals
origin <- read_csv("data/crossing_year.csv") 
origin <- origin[!origin %>% duplicated(), ]
stack <- merge(pvals1, origin, by = "accession_name")


# load("environment.RData")

library(ggpubr)
t1 <- stack %>% 
  filter(cross_year >= 2007) %>% 
  ggplot(aes(x = cross_year, y = predicted.value )) +
  geom_point()+
  stat_regline_equation()+
  # geom_abline(intercept = -320.462, slope = 0.1639, color="blue",
  #             size=1) +
  labs(x = "crossing year", y = "Predicted yield_ha") +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE) +
  theme_xiaofei()
  
t1
ggsave("images/GG_model4.png", plot = t1, units = "in", dpi = 300, width = 6, height = 6)

model4 <- lm(formula = predicted.value ~ cross_year, data = stack %>% filter(cross_year > 2006)  )
model5 <- lm(formula = predicted.value ~ cross_year, data = stack %>% filter(cross_year > 2006)  )

# -------------------------------------------------------------------------
# model 5
gg_model_5 <- agriutilities::parameters_gg(model = model5, trait = "yield_ha")
gg_model_5 %>% mutate(across(where(is.numeric), round, 4)) %>% 
  t() %>% as.data.frame() %>% 
  write.table("clipboard", col.names = F, row.names = T, sep = "\t")

# manual way
intercept = model5$coefficients[1]
slope = model5$coefficients[2]
lastYear = (2014*slope)  + intercept
firstYear = (2007*slope)  + intercept

slope/firstYear * 100  

# model 4

gg_model_4 <- agriutilities::parameters_gg(model = model4, trait = "yield_ha")
gg_model_4 %>% mutate(across(where(is.numeric), round, 4)) %>% 
  t() %>% as.data.frame() %>% 
  write.table("clipboard", col.names = F, row.names = T, sep = "\t")




