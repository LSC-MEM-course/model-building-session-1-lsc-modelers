# setwd("/Users/christinablomquist/Box/box-group-learning-to-talk/Experiments/CI_SP_Study/Data/ET/R_Analysis/")
#setwd("/Users/cblomquist/Box/box-group-learning-to-talk/Experiments/CI_SP_Study/Data/ET/R_Analysis/")
# wd is project directory

library(arm)
library(tidyverse)
library(eyetrackingR)
library(stringr)
library(lme4)
library(lmerTest)
library(apaTables)

# mean-centering function
mc <- function(x) scale(x, scale = F)[,1]


#### ---- Read in Data ---- ####

# Read in ET samples data
samples<- read_csv("data/cblomquist/samples.csv.gz")
# unique(samples$file)

# Read in participant data
part<-read_csv("data/cblomquist/CI_SP_participants.csv")

###############################

##### ---- Clean Eyetracking data ---- ####


## Filter to only good age and vocab matches ##
# Add ParticipantID column to samples table
samp <- samples %>%
  mutate(ParticipantID = str_extract(file, "^[:alpha:]{2}+[:digit:]{1,2}")) %>%
  mutate(ParticipantID = toupper(ParticipantID))
# Create table of ID and group match variables
gr_var1<- part %>%
  select(ParticipantID, VocabMatch, AGMatch, Condition)
colnames(gr_var1)<-c("ParticipantID", "VocabMatch", "AGMatch", "Group")
# Join samples and group variables by ParticipantID
samp_p <- samp %>%
  left_join(gr_var1)
# Filter to good age and vocab matches
samp_filt <- samp_p %>%
  filter(AGMatch == 1 | VocabMatch == 1)
# Get total number of trials
samp_total <- samp_filt %>%
  select(file, trial) %>%
  unique()
head(samp_p)
## Remove blocks with less than .9 accuracy ##
# Get trial messages
all_trial_msg <- samp_filt %>%
  select(file, trial, cond:respCode) %>%
  distinct()
# Create column for accuracy
trial_acc <- all_trial_msg %>%
  filter(cond != 3) %>%
  group_by(file) %>%
  summarise(
    target_pct = sum(respCode == "target", na.rm=T)/n()
  )
# Add accuracy column to samples
samples_acc <- samp_filt %>%
  left_join(trial_acc, by = "file") %>%
  filter(target_pct >= .9)
# Get table of blocks with acc <.9
bad_blocks <- trial_acc %>%
  filter(target_pct < .9)
  # ~No blocks with less than .9 accuracy~

## Filter to samples of interest ##
# Filter for target analysis...
d_t <- samples_acc %>%
  # ...to only samples where response was target
  filter(respCode == "target",
  # ...to remove filler trials (1188/5940 trials)
         cond!=3,
# ...to only samples within time window of interest for target (896-3000) --> window length is 2104 ms
         rel_time >= 896,
         rel_time <= 3000) %>%
  # *Define some new variables...
         # Define time from verb onset (with 200 ms delay)
  mutate(ons_time = (rel_time-896), 
         # Define trackloss
         trackloss = trackloss_post,
         # Define condition labels
         cond_clean = recode(cond, `1` = "Informative", `2` = "Neutral")) %>%
  # *Make sure aoi columns are logical
  mutate(target_aoi = as.logical(target_aoi),
         cohort_aoi = as.logical(cohort_aoi),
         urt_aoi = as.logical(urt_aoi),
         urc_aoi = as.logical(urc_aoi))

# # How many trials were lost due to non-target response?
# d_t_targfilt <- samples_acc %>%
#   filter(respCode == "target") # %>%
# d_t_targfilt_num <- d_t_targfilt %>%
#   select(file, trial) %>%
#   unique()
#   #5940 - 5824 = 116 trials removed
# 
# # How many additional trials removed due to being fillers?
# d_t_fillers <- d_t_targfilt %>%
#   filter(cond != 3) %>%
#   select(file, trial) %>%
#   unique()
#   #66*18 = 1188 total filler trials
#   #5824 - 4706 = 1118 filler trials removed in this step (70 were both filler trials and non-target responses)
# 
# # Total trials removed at this point?
# trial_table <- d_t %>%
#   select(file, trial) %>%
#   unique()
#   #5940 - 4706 = 1234 non-target and/or filler trials

# Filter for cohort analysis
d_c <- samples_acc %>%
  # ...to only samples where response was target
  filter(respCode == "target",
         # ...to remove filler trials 
         cond!=3,
         # ...to only samples within time window of interest for target (896-3000) --> window length is 2104 ms
         rel_time >= 1914,
         rel_time <= 3000) %>%
  # *Define some new variables... 
  # Define time from verb onset (with 200 ms delay)
  mutate(ons_time = (rel_time-1914), 
         # Define trackloss
         trackloss = trackloss_post,
         # Define condition labels
         cond_clean = recode(cond, `1` = "Informative", `2` = "Neutral")) %>%
  # *Make sure aoi columns are logical
  mutate(target_aoi = as.logical(target_aoi),
         cohort_aoi = as.logical(cohort_aoi),
         urt_aoi = as.logical(urt_aoi),
         urc_aoi = as.logical(urc_aoi))

# Format data for eyetrackingR, and remove trials with too much track loss
d_t_etr <- d_t %>%
  make_eyetrackingr_data(participant_column = "file", trackloss_column = "trackloss",
                         time_column = "ons_time", trial_column = "trial", 
                         aoi_columns = c("target_aoi", "cohort_aoi", "urt_aoi", "urc_aoi"), 
                         treat_non_aoi_looks_as_missing = FALSE, item_columns = "target") %>%
  # Clean by trial based on trackloss within analysis window 
  clean_by_trackloss(trial_prop_thresh = .5, # Allowing for 50% track loss
                     window_start_time = 0, 
                     window_end_time = 2104) # 896-3000 ms (3000 extends to max length of target word)
          # 3.8% of trials lost (178/4706) for 896-3000 ms with .5 trackloss threshold


d_c_etr <- d_c %>%
  make_eyetrackingr_data(participant_column = "file", trackloss_column = "trackloss",
                         time_column = "ons_time", trial_column = "trial", 
                         aoi_columns = c("target_aoi", "cohort_aoi", "urt_aoi", "urc_aoi"), 
                         treat_non_aoi_looks_as_missing = FALSE, item_columns = "target") %>%
  # Clean by trial based on trackloss within analysis window 
  clean_by_trackloss(trial_prop_thresh = .5, # Allowing for 50% track loss
                     window_start_time = 0, 
                     window_end_time = 1086) # 1914-3000 ms (3000 extends to max length of target word)
          # 6.4% of trials lost (302/4706) for 1914-3000 ms with .5 trackloss threshold

# Clean for trial loss by participant 
  # Remove blocks with <9 trials per condition (half of the 18 possible trials)
# Clean target data #
d_t_etr_calc <- d_t_etr %>%
  select(file, trial, cond_clean, ParticipantID, Group) %>%
  distinct()
p <-unique(d_t_etr_calc$file)
d_t_etr_calc$NumP <- NA
d_t_etr_calc$NumN <- NA
for (i in 1:length(p)){
  f <- p[i]
  d_t_etr_calc$NumP <- ifelse(d_t_etr_calc$file == f, nrow(subset(d_t_etr_calc, d_t_etr_calc$file == f & d_t_etr_calc$cond_clean == "Informative")), d_t_etr_calc$NumP)
  d_t_etr_calc$NumN <- ifelse(d_t_etr_calc$file == f, nrow(subset(d_t_etr_calc, d_t_etr_calc$file == f & d_t_etr_calc$cond_clean == "Neutral")), d_t_etr_calc$NumN)
}
# Get list of bad (removed) blocks
d_t_etr_b <- d_t_etr_calc %>%
  filter(NumP<9 | NumN<9) %>%
  select (file, ParticipantID, Group, NumP, NumN) %>%
  distinct()
# When <9 (and trial thresh = .5), no blocks removed

# Select only good blocks
d_t_etr_g<- d_t_etr_calc %>%
  filter(NumP>=9 & NumN>=9) %>%
  select (file, ParticipantID, Group, NumP, NumN) %>%
  distinct()
# Get data only from good blocks
d_t_etr_p<- d_t_etr %>%
  inner_join(d_t_etr_g)

# Clean cohort data #
d_c_etr_calc <- d_c_etr %>%
  select(file, trial, cond_clean, ParticipantID, Group) %>%
  distinct()
p <-unique(d_c_etr_calc$file)
d_c_etr_calc$NumP <- NA
d_c_etr_calc$NumN <- NA
for (i in 1:length(p)){
  f <- p[i]
  d_c_etr_calc$NumP <- ifelse(d_c_etr_calc$file == f, nrow(subset(d_c_etr_calc, d_c_etr_calc$file == f & d_c_etr_calc$cond_clean == "Informative")), d_c_etr_calc$NumP)
  d_c_etr_calc$NumN <- ifelse(d_c_etr_calc$file == f, nrow(subset(d_c_etr_calc, d_c_etr_calc$file == f & d_c_etr_calc$cond_clean == "Neutral")), d_c_etr_calc$NumN)
}
# Get list of bad (removed) blocks
d_c_etr_b <- d_c_etr_calc %>%
  filter(NumP<9 | NumN<9) %>%
  select (file, ParticipantID, Group, NumP, NumN) %>%
  distinct()
# When <9 (and trial thresh = .5), no blocks removed

# Select only good blocks
d_c_etr_g<- d_c_etr_calc %>%
  filter(NumP>=9 & NumN>=9) %>%
  select (file, ParticipantID, Group, NumP, NumN) %>%
  distinct()
# Get data only from good blocks
d_c_etr_p<- d_c_etr %>%
  inner_join(d_c_etr_g)

# # Save tables
# write_csv(d_t_etr_p, "d_etr_verbtw.csv.gz")
# write_csv(d_c_etr_p, "d_etr_nountw.csv.gz")

# ## Sanity check --> Summarize data ##
# # Make target column a factor instead of logical
# d_t_etr_p$target_look <- as.integer(d_t_etr_p$target_aoi)
# # Average across condition and file
# summ <- describe_data(d_t_etr_p,
#                       describe_column = 'target_look',
#                       group_columns = c('cond_clean','ParticipantID'))
# # Plot
# plot(summ)

# Make time window data - aggregated by participant
  # target 
tw_t <- make_time_window_data(d_t_etr_p, 
                                 predictor_columns = c("Group", "cond_clean", "audstim", "setid"),
                                 summarize_by = c("trial", "file", "ParticipantID"))
  # cohort
d_c_etr_pc <- d_c_etr_p %>%
  select(-target_aoi)
tw_c <- make_time_window_data(d_c_etr_pc,
                                  aois = c("cohort_aoi", "urt_aoi", "urc_aoi"),
                                 predictor_columns = c("Group", "cond_clean", "audstim", "setid"),
                                 summarize_by = c("trial", "file", "ParticipantID"))

# # Plot data
# tw_t <- subset(tw, tw$AOI == "target_aoi")
# plot(tw_t, predictor_columns="cond_clean", dv = "LogitAdjusted")
# plot(tw_t, predictor_columns="cond_clean", dv = "Prop")
# plot(tw_t, predictor_columns="cond_clean", dv = "Elog")

# Save cleaned time window data to csv #
write_csv(tw_t, "tw_clean_Target.csv.gz")
write_csv(tw_c, "tw_clean_Cohort.csv.gz")

#############################

#### ---- Analysis ---- ####
# Read in data
tw_t <- read_csv("tw_clean_Target.csv.gz")
tw_c <- read_csv("tw_clean_Cohort.csv.gz") 
part<-read_csv("data/cblomquist/CI_SP_participants.csv")

# Filter participant table to info of interest (change AgeGroup to 0s (Y) and 1s (O))
part_var <- part %>%
  select(ParticipantID, Condition, Gender, Handedness, Age_mos, AGMatch, VocabMatch,
         EVT_SS, EVT_Perc, EVT_GSV, 
         PPVT_SS, PPVT_Perc, PPVT_GSV, 
         KBIT_Raw, KBIT_SS, KBIT_Perc) %>%
  mutate(Group = Condition) %>%
  select(-Condition)
# Join tables, assign numerical values to condition
tw_t_pv <- tw_t %>%
  inner_join(part_var) %>%
  mutate(cond_clean = ifelse(cond_clean == "Neutral", 0, cond_clean),
         cond_clean = ifelse(cond_clean == "Informative", 1, cond_clean))
tw_c_pv <- tw_c %>%
  inner_join(part_var) %>%
  mutate(cond_clean = ifelse(cond_clean == "Neutral", 0, cond_clean),
         cond_clean = ifelse(cond_clean == "Informative", 1, cond_clean))

#### Group comparison for TARGET - AGE Matches ####
# Filter 
tw_AG <- tw_t_pv %>%
  filter(AGMatch == 1 )

# Filter to just looks to target AOI
tw_AG_T<- tw_AG %>%
  filter(AOI == "target_aoi")

#Edit weights
weights_barr2012 <- function(y, N, offset = 0.5) {
  1 / (1/(y + offset) + 1/(N - y + offset))
}
tw_AG_T <- tw_AG_T %>%
  mutate(WeightsBarr = weights_barr2012(SamplesInAOI, SamplesTotal))

# Make vocab datasets (mean-centered, and accounting for if some pts missing certain test scores)
tw_AG_T_ID<- tw_AG_T %>%
  select(ParticipantID, PPVT_GSV, EVT_GSV, KBIT_SS) %>%
  unique()
tw_AG_T$EVT_GSV_C<- mc(tw_AG_T$EVT_GSV)
tw_AG_T$PPVT_GSV_C<- mc(tw_AG_T$PPVT_GSV)
tw_AG_T$PPVT_SS_C<- mc(tw_AG_T$PPVT_SS)
tw_AG_T$KBIT_SS_C<- mc(tw_AG_T$KBIT_SS)
tw_AG_T$AgeMos_C<- mc(tw_AG_T$Age_mos)
tw_AG_Tv<- tw_AG_T %>%
  filter(!is.na(PPVT_GSV))
tw_AG_Tev<- tw_AG_T %>%
  filter(!is.na(EVT_GSV))
tw_AG_Tk<- tw_AG_T %>%
  filter(!is.na(KBIT_SS))

### starting some class modeling here
cdata <- d_t_etr_p
head(cdata) %>% as.data.frame()
xtabs(~ cond, cdata)

cdata_sub <- mutate(cdata, 
                look_to_target = as.numeric(target_aoi)) %>%
  group_by(trial, ParticipantID) %>% 
  mutate(yprev = lag(look_to_target)) %>% ungroup() %>%
  filter(ParticipantID %in% unique(cdata$ParticipantID)[25:34]) %>%
  filter(!is.na(yprev))

glm0 <- glm(look_to_target ~ 1,
            family = "binomial", data = cdata_sub)
summary(glm0)

glm1 <- glm(look_to_target ~ 1 + yprev + cond,
            family = "binomial", data = cdata_sub)
summary(glm1)



# jump to collapsed data
colldata <- tw_AG_T %>% 
  mutate(PPVT_scaled = as.numeric(scale(PPVT_GSV_C)))

glm0 <- glm(cbind(SamplesInAOI, SamplesTotal - SamplesInAOI) ~
              1,
            family = "binomial", data = colldata)
summary(glm0)

glm1 <- glm(cbind(SamplesInAOI, SamplesTotal - SamplesInAOI) ~
               1 + cond_clean * Group,
             family = "binomial", data = colldata)
summary(glm1)

glm2 <- glm(cbind(SamplesInAOI, SamplesTotal - SamplesInAOI) ~
              1 + cond_clean * Group * PPVT_scaled,
            family = "binomial", data = colldata)
summary(glm2)
summary(colldata$PPVT_GSV_C)

glmer2_a <- 
  glmer(cbind(SamplesInAOI, SamplesTotal - SamplesInAOI) ~
          1 + cond_clean * Group * PPVT_scaled + 
          (1 | ParticipantID) +
          (1 | setid),
          family = "binomial", data = colldata)
summary(glmer2_a, corr = FALSE)

AIC(glm2)
AIC(glmer2_a)

glmer2_a <- 
  glmer(cbind(SamplesInAOI, SamplesTotal - SamplesInAOI) ~
          1 + cond_clean * Group * PPVT_scaled + 
          (1 | ParticipantID) +
          (1 | setid),
        family = "binomial", data = colldata)
summary(glmer2_a, corr = FALSE)

glmer2_b <- 
  glmer(cbind(SamplesInAOI, SamplesTotal - SamplesInAOI) ~
          1 + cond_clean * Group * PPVT_scaled + 
          (1 + cond_clean | ParticipantID) +
          (1 | setid),
        family = "binomial", data = colldata)
summary(glmer2_b, corr = FALSE)

glmer2_c <- 
  glmer(cbind(SamplesInAOI, SamplesTotal - SamplesInAOI) ~
          1 + cond_clean * Group * PPVT_scaled + 
          (1 | ParticipantID) +
          (1 + cond_clean * Group * PPVT_scaled | setid),
        family = "binomial", data = colldata,
        control = glmerControl(optimizer = "bobyqa",
                               optCtrl = list(maxfun = 50000)))
summary(glmer2_c, corr = FALSE)

glmer2_d <- 
  glmer(cbind(SamplesInAOI, SamplesTotal - SamplesInAOI) ~
          1 + cond_clean * Group * PPVT_scaled + 
          (1 + cond_clean | ParticipantID) +
          (1 + cond_clean * Group * PPVT_scaled | setid),
        family = "binomial", data = colldata,
        control = glmerControl(optimizer = "bobyqa",
                               optCtrl = list(maxfun = 50000)))
summary(glmer2_d, corr = FALSE)


# Base model
lmb <- lmer(Elog ~ cond_clean +  (1 | setid) + (1 | ParticipantID),
            data = tw_AG_T, weights = WeightsBarr,
            control=lmerControl(optimizer="bobyqa", optCtrl = list(maxfun=4000000)), REML = FALSE)
summary(lmb)
# With random slopes
lmb_rs <- lmer(Elog ~ cond_clean +  (cond_clean | setid) + (cond_clean | ParticipantID),
            data = tw_AG_T, weights = WeightsBarr,
            control=lmerControl(optimizer="bobyqa", optCtrl = list(maxfun=4000000)), REML = FALSE)
summary(lmb_rs) # *Fit is singular

# Add group
lm1 <- lmer(Elog ~ cond_clean * Group + (1 | setid) + (1 | ParticipantID),
            data = tw_AG_T, weights = WeightsBarr,
            control=lmerControl(optimizer="bobyqa", optCtrl = list(maxfun=4000000)), REML = FALSE)
summary(lm1)
anova(lmb, lm1)
# Does adding group improve model? - Yes
# With random slopes
lm1_rs <- lmer(Elog ~ cond_clean * Group +  (cond_clean | setid) + (cond_clean | ParticipantID),
               data = tw_AG_T, weights = WeightsBarr,
               control=lmerControl(optimizer="bobyqa", optCtrl = list(maxfun=4000000)), REML = FALSE)
summary(lm1_rs) # *Fit is singular
anova(lmb_rs, lm1_rs) 

# Add PPVT
lm2a <- lmer(Elog ~ cond_clean*Group*PPVT_GSV_C + (1 | setid) + (1 | ParticipantID),
             data = tw_AG_T, weights = WeightsBarr, 
             control=lmerControl(optimizer="bobyqa", optCtrl = list(maxfun=4000000)), REML = FALSE)
summary(lm2a)
anova(lm1, lm2a)
# Does adding PPVT improve model? - Yes
# With random slopes
lm2a_rs <- lmer(Elog ~ cond_clean * Group * PPVT_GSV_C +  (cond_clean | setid) + (cond_clean | ParticipantID),
               data = tw_AG_T, weights = WeightsBarr,
               control=lmerControl(optimizer="bobyqa", optCtrl = list(maxfun=4000000)), REML = FALSE)
summary(lm2a_rs) # *Fit is singular
anova(lm1_rs, lm2a_rs)

# Add EVT
lm2b <- lmer(Elog ~ cond_clean*Group*EVT_GSV_C + (1 | setid) + (1 | ParticipantID),
             data = tw_AG_Tev, weights = WeightsBarr, 
             control=lmerControl(optimizer="bobyqa", optCtrl = list(maxfun=4000000)), REML = FALSE)
lm1b <- lmer(Elog ~ cond_clean*Group + (1 | setid) + (1 | ParticipantID),
             data = tw_AG_Tev, weights = WeightsBarr, 
             control=lmerControl(optimizer="bobyqa", optCtrl = list(maxfun=4000000)), REML = FALSE)
summary(lm2b)
anova(lm1b, lm2b)
# Does adding EVT improve model? - No



#### Group comparison for COHORT - AGE Matches ####
# Filter
twc_AG <- tw_c_pv %>%
  filter(AGMatch == 1 )

# Filter to just looks to target AOI
tw_AG_C<- tw_AG %>%
  filter(AOI == "cohort_aoi")

#Edit weights
weights_barr2012 <- function(y, N, offset = 0.5) {
  1 / (1/(y + offset) + 1/(N - y + offset))
}
tw_AG_C <- tw_AG_C %>%
  mutate(WeightsBarr = weights_barr2012(SamplesInAOI, SamplesTotal))

# Make vocab datasets (mean-centered, and accounting for if some pts missing vocab)
# Add vocab -- mean centered
tw_AG_C_ID<- tw_AG_C %>%
  select(ParticipantID, PPVT_GSV, EVT_GSV, KBIT_SS) %>%
  unique()
tw_AG_C$EVT_GSV_C<- mc(tw_AG_C$EVT_GSV)
tw_AG_C$PPVT_GSV_C<- mc(tw_AG_C$PPVT_GSV)
tw_AG_C$PPVT_SS_C<- mc(tw_AG_C$PPVT_SS)
tw_AG_C$KBIT_SS_C<- mc(tw_AG_C$KBIT_SS)
tw_AG_C$AgeMos_C<- mc(tw_AG_C$Age_mos)
tw_AG_Cv<- tw_AG_C %>%
  filter(!is.na(PPVT_GSV))
tw_AG_Cev<- tw_AG_C %>%
  filter(!is.na(EVT_GSV))
tw_AG_Ck<- tw_AG_C %>%
  filter(!is.na(KBIT_SS))

# Base model
# Base model
lmbc <- lmer(Elog ~ cond_clean +  (1 | setid) + (1 | ParticipantID),
            data = tw_AG_C, weights = WeightsBarr,
            control=lmerControl(optimizer="bobyqa", optCtrl = list(maxfun=4000000)), REML = FALSE)
summary(lmb)
# With random slopes
lmbc_rs <- lmer(Elog ~ cond_clean +  (cond_clean | setid) + (cond_clean | ParticipantID),
               data = tw_AG_C, weights = WeightsBarr,
               control=lmerControl(optimizer="bobyqa", optCtrl = list(maxfun=4000000)), REML = FALSE)
summary(lmbc_rs) 

# Add group
lm1c <- lmer(Elog ~ cond_clean * Group + (1 | setid) + (1 | ParticipantID),
            data = tw_AG_C, weights = WeightsBarr,
            control=lmerControl(optimizer="bobyqa", optCtrl = list(maxfun=4000000)), REML = FALSE)
summary(lm1c)
anova(lmbc, lm1c)
# Does adding group improve model? - Yes
# With random slopes
lm1c_rs <- lmer(Elog ~ cond_clean * Group +  (cond_clean | setid) + (cond_clean | ParticipantID),
               data = tw_AG_C, weights = WeightsBarr,
               control=lmerControl(optimizer="bobyqa", optCtrl = list(maxfun=4000000)), REML = FALSE)
summary(lm1c_rs) 
anova(lmbc_rs, lm1c_rs) 

# Add PPVT
lm2ac <- lmer(Elog ~ cond_clean*Group*PPVT_GSV_C + (1 | setid) + (1 | ParticipantID),
             data = tw_AG_C, weights = WeightsBarr, 
             control=lmerControl(optimizer="bobyqa", optCtrl = list(maxfun=4000000)), REML = FALSE)
summary(lm2ac)
anova(lm1c, lm2ac)
# Does adding PPVT improve model? - Yes
# With random slopes
lm2ac_rs <- lmer(Elog ~ cond_clean * Group * PPVT_GSV_C +  (cond_clean | setid) + (cond_clean | ParticipantID),
                data = tw_AG_C, weights = WeightsBarr,
                control=lmerControl(optimizer="bobyqa", optCtrl = list(maxfun=4000000)), REML = FALSE)
summary(lm2ac_rs) 
anova(lm1c_rs, lm2ac_rs)

# Add EVT
lm2bc <- lmer(Elog ~ cond_clean*Group*EVT_GSV_C + (1 | setid) + (1 | ParticipantID),
             data = tw_AG_Cev, weights = WeightsBarr, 
             control=lmerControl(optimizer="bobyqa", optCtrl = list(maxfun=4000000)), REML = FALSE)
lm1bc <- lmer(Elog ~ cond_clean*Group + (1 | setid) + (1 | ParticipantID),
             data = tw_AG_Cev, weights = WeightsBarr, 
             control=lmerControl(optimizer="bobyqa", optCtrl = list(maxfun=4000000)), REML = FALSE)
summary(lm2bc)
anova(lm1bc, lm2bc)
# Does adding EVT improve model? - Yes
# With random slopes
lm2bc_rs <- lmer(Elog ~ cond_clean * Group * EVT_GSV_C +  (cond_clean | setid) + (cond_clean | ParticipantID),
                data = tw_AG_Cev, weights = WeightsBarr,
                control=lmerControl(optimizer="bobyqa", optCtrl = list(maxfun=4000000)), REML = FALSE)
lm1bc_rs <- lmer(Elog ~ cond_clean * Group +  (cond_clean | setid) + (cond_clean | ParticipantID),
                data = tw_AG_Cev, weights = WeightsBarr,
                control=lmerControl(optimizer="bobyqa", optCtrl = list(maxfun=4000000)), REML = FALSE)
summary(lm2bc_rs) 
anova(lm1bc_rs, lm2bc_rs)


#### Group comparison for TARGET - VOCAB Matches ####
# Filter 
tw_V <- tw_t_pv %>%
  filter(VocabMatch == 1 )

# Filter to just looks to target AOI
tw_V_T<- tw_V %>%
  filter(AOI == "target_aoi")

#Edit weights
weights_barr2012 <- function(y, N, offset = 0.5) {
  1 / (1/(y + offset) + 1/(N - y + offset))
}
tw_V_T <- tw_V_T %>%
  mutate(WeightsBarr = weights_barr2012(SamplesInAOI, SamplesTotal))

# Make vocab datasets (mean-centered, and accounting for if some pts missing vocab)
# Add vocab -- mean centered
tw_V_T_ID<- tw_V_T %>%
  select(ParticipantID, PPVT_GSV, EVT_GSV, KBIT_SS) %>%
  unique()
tw_V_T$EVT_GSV_C<- mc(tw_V_T$EVT_GSV)
tw_V_T$PPVT_GSV_C<- mc(tw_V_T$PPVT_GSV)
tw_V_T$PPVT_SS_C<- mc(tw_V_T$PPVT_SS)
tw_V_T$KBIT_SS_C<- mc(tw_V_T$KBIT_SS)
tw_V_T$AgeMos_C<- mc(tw_V_T$Age_mos)
tw_V_Tv<- tw_V_T %>%
  filter(!is.na(PPVT_GSV))
tw_V_Tev<- tw_V_T %>%
  filter(!is.na(EVT_GSV))
tw_V_Tk<- tw_V_T %>%
  filter(!is.na(KBIT_SS))

# Base model
lmbv <- lmer(Elog ~ cond_clean +  (1 | setid) + (1 | ParticipantID),
            data = tw_V_T, weights = WeightsBarr,
            control=lmerControl(optimizer="bobyqa", optCtrl = list(maxfun=4000000)), REML = FALSE)
summary(lmbv)
# With random slopes
lmbv_rs <- lmer(Elog ~ cond_clean +  (cond_clean | setid) + (cond_clean | ParticipantID),
               data = tw_V_T, weights = WeightsBarr,
               control=lmerControl(optimizer="bobyqa", optCtrl = list(maxfun=4000000)), REML = FALSE)
summary(lmbv_rs) # *Fit is singular

# Add group
lm1v <- lmer(Elog ~ cond_clean * Group + (1 | setid) + (1 | ParticipantID),
            data = tw_V_T, weights = WeightsBarr,
            control=lmerControl(optimizer="bobyqa", optCtrl = list(maxfun=4000000)), REML = FALSE)
summary(lm1v)
anova(lmbv, lm1v)
# Does adding group improve model? - Yes
# With random slopes
lm1v_rs <- lmer(Elog ~ cond_clean * Group +  (cond_clean | setid) + (cond_clean | ParticipantID),
               data = tw_V_T, weights = WeightsBarr,
               control=lmerControl(optimizer="bobyqa", optCtrl = list(maxfun=4000000)), REML = FALSE)
summary(lm1v_rs) # *Fit is singular
anova(lmbv_rs, lm1v_rs) 

# Add PPVT
lm2av <- lmer(Elog ~ cond_clean*Group*PPVT_GSV_C + (1 | setid) + (1 | ParticipantID),
             data = tw_V_T, weights = WeightsBarr, 
             control=lmerControl(optimizer="bobyqa", optCtrl = list(maxfun=4000000)), REML = FALSE)
summary(lm2av)
anova(lm1v, lm2av)
# Does adding PPVT improve model? - Yes
# With random slopes
lm2av_rs <- lmer(Elog ~ cond_clean * Group * PPVT_GSV_C +  (cond_clean | setid) + (cond_clean | ParticipantID),
                data = tw_V_T, weights = WeightsBarr,
                control=lmerControl(optimizer="bobyqa", optCtrl = list(maxfun=4000000)), REML = FALSE)
summary(lm2av_rs) # *Fit is singular
anova(lm1v_rs, lm2av_rs) # Does not improve model when there are random slopes

# Add EVT
lm2bv <- lmer(Elog ~ cond_clean*Group*EVT_GSV_C + (1 | setid) + (1 | ParticipantID),
             data = tw_V_Tev, weights = WeightsBarr, 
             control=lmerControl(optimizer="bobyqa", optCtrl = list(maxfun=4000000)), REML = FALSE)
lm1bv <- lmer(Elog ~ cond_clean*Group + (1 | setid) + (1 | ParticipantID),
             data = tw_V_Tev, weights = WeightsBarr, 
             control=lmerControl(optimizer="bobyqa", optCtrl = list(maxfun=4000000)), REML = FALSE)
summary(lm2bv)
anova(lm1bv, lm2bv)
# Does adding EVT improve model? - No

#### Group comparison for COHORT - VOCAB Matches ####
# Filter 
twc_V <- tw_c_pv %>%
  filter(VocabMatch == 1 )

# Filter to just looks to target AOI
tw_V_C<- twc_V %>%
  filter(AOI == "cohort_aoi")

#Edit weights
weights_barr2012 <- function(y, N, offset = 0.5) {
  1 / (1/(y + offset) + 1/(N - y + offset))
}
tw_V_C <- tw_V_C %>%
  mutate(WeightsBarr = weights_barr2012(SamplesInAOI, SamplesTotal))

# Make vocab datasets (mean-centered, and accounting for if some pts missing vocab)
# Add vocab -- mean centered
tw_V_C_ID<- tw_V_C %>%
  select(ParticipantID, PPVT_GSV, EVT_GSV, KBIT_SS) %>%
  unique()
tw_V_C$EVT_GSV_C<- mc(tw_V_C$EVT_GSV)
tw_V_C$PPVT_GSV_C<- mc(tw_V_C$PPVT_GSV)
tw_V_C$PPVT_SS_C<- mc(tw_V_C$PPVT_SS)
tw_V_C$KBIT_SS_C<- mc(tw_V_C$KBIT_SS)
tw_V_C$AgeMos_C<- mc(tw_V_C$Age_mos)
tw_V_Cv<- tw_V_C %>%
  filter(!is.na(PPVT_GSV))
tw_V_Cev<- tw_V_C %>%
  filter(!is.na(EVT_GSV))
tw_V_Ck<- tw_V_C %>%
  filter(!is.na(KBIT_SS))

# Base model
lmbvc <- lmer(Elog ~ cond_clean +  (1 | setid) + (1 | ParticipantID),
             data = tw_V_C, weights = WeightsBarr,
             control=lmerControl(optimizer="bobyqa", optCtrl = list(maxfun=4000000)), REML = FALSE)
summary(lmbvc)
# With random slopes
lmbvc_rs <- lmer(Elog ~ cond_clean +  (cond_clean | setid) + (cond_clean | ParticipantID),
                data = tw_V_C, weights = WeightsBarr,
                control=lmerControl(optimizer="bobyqa", optCtrl = list(maxfun=4000000)), REML = FALSE)
summary(lmbvc_rs) # *Fit is singular

# Add group
lm1vc <- lmer(Elog ~ cond_clean * Group + (1 | setid) + (1 | ParticipantID),
             data = tw_V_C, weights = WeightsBarr,
             control=lmerControl(optimizer="bobyqa", optCtrl = list(maxfun=4000000)), REML = FALSE)
summary(lm1vc)
anova(lmbvc, lm1vc)
# Does adding group improve model? - Yes
# With random slopes
lm1vc_rs <- lmer(Elog ~ cond_clean * Group +  (cond_clean | setid) + (cond_clean | ParticipantID),
                data = tw_V_C, weights = WeightsBarr,
                control=lmerControl(optimizer="bobyqa", optCtrl = list(maxfun=4000000)), REML = FALSE)
summary(lm1vc_rs) # *Fit is singular
anova(lmbvc_rs, lm1vc_rs) 

# Add PPVT
lm2avc <- lmer(Elog ~ cond_clean*Group*PPVT_GSV_C + (1 | setid) + (1 | ParticipantID),
              data = tw_V_C, weights = WeightsBarr, 
              control=lmerControl(optimizer="bobyqa", optCtrl = list(maxfun=4000000)), REML = FALSE)
summary(lm2avc)
anova(lm1vc, lm2avc)
# Does adding PPVT improve model? - Yes
# With random slopes
lm2avc_rs <- lmer(Elog ~ cond_clean * Group * PPVT_GSV_C +  (cond_clean | setid) + (cond_clean | ParticipantID),
                 data = tw_V_C, weights = WeightsBarr,
                 control=lmerControl(optimizer="bobyqa", optCtrl = list(maxfun=4000000)), REML = FALSE)
summary(lm2avc_rs) # *Fit is singular
anova(lm1vc_rs, lm2avc_rs) # Does not improve model when there are random slopes

# Add EVT
lm2bvc <- lmer(Elog ~ cond_clean*Group*EVT_GSV_C + (1 | setid) + (1 | ParticipantID),
              data = tw_V_Cev, weights = WeightsBarr, 
              control=lmerControl(optimizer="bobyqa", optCtrl = list(maxfun=4000000)), REML = FALSE)
lm1bvc <- lmer(Elog ~ cond_clean*Group + (1 | setid) + (1 | ParticipantID),
              data = tw_V_Cev, weights = WeightsBarr, 
              control=lmerControl(optimizer="bobyqa", optCtrl = list(maxfun=4000000)), REML = FALSE)
summary(lm2bvc)
anova(lm1bvc, lm2bvc)
# Does adding EVT improve model? - No

