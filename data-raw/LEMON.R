## code to prepare `LEMON` dataset goes here

library(readr)
library(dplyr)
library(stringr)

data_path <- system.file("extdata", package = "componentts")

# cognitive tasks

cvlt <- read_csv(paste0(data_path, "/cognitive/CVLT.csv"))
tapa <- read_csv(paste0(data_path, "/cognitive/TAP-Alertness.csv"))
tapi <- read_csv(paste0(data_path, "/cognitive/TAP-Incompatibility.csv"))
tapwm <- read_csv(paste0(data_path, "/cognitive/TAP-Working-Memory.csv"))
tmt <- read_csv(paste0(data_path, "/cognitive/TMT.csv"))
wst <- read_csv(paste0(data_path, "/cognitive/WST.csv"))
lps <- read_csv(paste0(data_path, "/cognitive/LPS.csv"))
rwt <- read_csv(paste0(data_path, "/cognitive/RWT.csv"))

# clean cvlt

cvlt <- cvlt %>% select(-CVLT_1, -CVLT_7, -CVLT_16) # missing data
cvlt$CVLT_15 <- cvlt$CVLT_15 * -1 # reverse negative measure
cvlt$CVLT_14 <- cvlt$CVLT_14 * -1 # reverse negative measure
cvlt$CVLT_4 <- cvlt$CVLT_4 * -1 # reverse negative measure
cvlt$CVLT_5 <- cvlt$CVLT_5 * -1 # reverse negative measure
cvlt <- cvlt %>% select(-CVLT_8, -CVLT_9, -CVLT_10, -CVLT_11, -CVLT_12) # highly correlated variables

# clean tmt

tmt <- tmt %>% select(-TMT_4, -TMT_8) # comments
tmt <- tmt %>% select(-TMT_2, -TMT_6) # categorical "brain function"
tmt <- tmt %>% mutate(across(where(is.numeric), `-`)) # reverse negative measures

# clean wst

wst <- wst %>% select(ID, WST_1) # keep only WST raw score, other columns are normalized scores

# clean lps

lps <- lps %>% select(-LPS_2)

# clean tapa

tapa[,2:17] <- sapply(tapa[,2:17], FUN = as.numeric) # chars become NA
tapa <- tapa[,-c(8, 10, 13, 15, 17, 18)] # percentages and comments
tapa <- tapa %>% select(-TAP_A_5, -TAP_A_10, -TAP_A_8, -TAP_A_13) # remove means and SDs, keeping medians only
tapa <- tapa %>% select(-TAP_A_1, -TAP_A_2, -TAP_A_3, -TAP_A_4) #remove individual trials
tapa$TAP_A_6 <- tapa$TAP_A_6 * -1 # reverse negative measures
tapa$TAP_A_11 <- tapa$TAP_A_11 * -1 # reverse negative measures

# clean tapi

tapi[,2:27] <- sapply(tapi[,2:27], FUN = as.numeric) #chars become NA
tapi <- tapi[,-c(4, 6, 8, 11, 13, 15, 18, 20, 22, 24, 26, 28, 29)] #percentages and comments
tapi <- tapi %>% select(-TAP_I_1, -TAP_I_8, -TAP_I_15, -TAP_I_4, -TAP_I_11, -TAP_I_18) # remove means and SDs, keep medians only
tapi <- tapi %>% select(-TAP_I_16, -TAP_I_20) #remove aggregate of compatible and incompatible
tapi <- tapi %>% select(-TAP_I_22, -TAP_I_24, -TAP_I_26) # remove F values
tapi <- tapi %>% mutate(across(where(is.numeric), `-`)) # reverse negative measures

# clean tapwm

tapwm <- tapwm[,-c(4, 6, 9, 11, 13)] #percentages and comments
tapwm <- tapwm %>% select(-TAP_WM_1, -TAP_WM_4) #remove mean and SD, keep median only
tapwm <- tapwm %>% select(-TAP_WM_9, -TAP_WM_11) #remove outliers and missed matches, perfect correlate with correct matches
tapwm$TAP_WM_2 <- tapwm$TAP_WM_2 * -1 # reverse negative measures
tapwm$TAP_WM_7 <- tapwm$TAP_WM_7 * -1

#clean rwt

rwt[,2:ncol(rwt)] <- sapply(rwt[,2:ncol(rwt)], FUN = as.numeric) # chars become NA
rwt <- rwt[,-c(2:8, 14:20)] # remove first and second minute data
rwt <- rwt %>% select(-RWT_12, -RWT_24) # remove comments
rwt <- rwt %>% select(-RWT_9, -RWT_21) # remove percentiles
rwt$RWT_10 <- rwt$RWT_10 * -1 #reverse negative measures
rwt$RWT_11 <- rwt$RWT_11 * -1 #reverse negative measures
rwt$RWT_22 <- rwt$RWT_22 * -1 #reverse negative measures
rwt$RWT_23 <- rwt$RWT_23 * -1 #reverse negative measures

cognitive <- cvlt %>% full_join(tmt, by = "ID") %>%
  full_join(wst, by = "ID") %>%
  full_join(lps, by = "ID") %>%
  full_join(tapa, by = "ID") %>%
  full_join(tapi, by = "ID") %>%
  full_join(tapwm, by = "ID") %>%
  full_join(rwt, by = "ID")

cognitive <- cognitive[complete.cases(cognitive),]

# medical
anthro <- read_csv(paste0(data_path, "/medical/Anthropometry_LEMON.csv"))
bp <- read_csv(paste0(data_path, "/medical/Blood_Pressure_LEMON.csv"))
blood <- read_csv(paste0(data_path, "/medical/Blood_Results_LEMON.csv"))

# anthro is good out of the box

#clean bp

bp <- bp %>% select(-pulse1_right, -BP1_right_systole, -BP1_right_diastole) # remove redundant pulse and bp measures

# clean blood

blood <- blood %>% select(-contains("Reference"), -contains("Rerefence"),
                          -Date_Blood_Drawing_LabAnalysis)
blood$CRP_in_mg_l <- as.numeric(blood$CRP_in_mg_l)
blood <- blood %>% select(-HGBK_in_g_dl, -MCHCK_in_g_dl, -MCHK_in_pg,
                          -`HBA1C_in_%`, -HBA1CI_in_mmol_mol) # remove duplicate measures in different units
na_per_col <- sapply(blood, function(x) sum(is.na(x)))
blood <- blood[,-(which(na_per_col > 12))]
colnames(blood)[2:ncol(blood)] <- str_split(colnames(blood)[2:ncol(blood)],
                                            "_", simplify = TRUE)[,1] # shorten variable names

medical <- anthro %>% full_join(bp, by = "ID") %>%
  full_join(blood, by = "ID")

medical <- medical[complete.cases(medical),]

# Demographics

demo <- read_csv(paste0(data_path,
                        "/META_File_IDs_Age_Gender_Education_Drug_Smoke_SKID_LEMON.csv"))
# create Gender char
demo <- demo %>% mutate(Gender = recode(`Gender_ 1=female_2=male`,
                                        `1` = "Female", `2` = "Male"),
                        .after = `Gender_ 1=female_2=male`)

demo <- demo %>% mutate(Age_bin = recode(Age,
                                 "20-25" = "Young",
                                 "25-30" = "Young",
                                 "30-35" = "Young",
                                 "35-40" = "Young",
                                 "55-60" = "Mature",
                                 "60-65" = "Mature",
                                 "65-70" = "Mature",
                                 "70-75" = "Mature",
                                 "75-80" = "Mature"))

demo <- demo %>% select(ID, Gender, Age_bin)

mutual_ID <- (cognitive %>% inner_join(medical, by = "ID") %>%
  inner_join(demo, by = "ID"))$ID

cognitive <- cognitive %>% filter(ID %in% mutual_ID) %>% arrange(ID)
medical <- medical %>% filter(ID %in% mutual_ID) %>% arrange(ID)
demo <- demo %>% filter(ID %in% mutual_ID) %>% arrange(ID)

write_csv(cognitive, paste0(data_path, "/cleaned/cognitive.csv"))
write_csv(medical, paste0(data_path, "/cleaned/medical.csv"))
write_csv(demo, paste0(data_path, "/cleaned/demo.csv"))

LEMON <- list(cognitive = cognitive,
              medical = medical,
              demo = demo)

usethis::use_data(LEMON, overwrite = TRUE)
