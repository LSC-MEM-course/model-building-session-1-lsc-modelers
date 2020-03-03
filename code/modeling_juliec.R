### Read me: This contains data cleaning and ggplots for SU and ID plots with SE for YN and ON.
### USE ME!

#Last updated: 1/15/2020

library(arm)
library(readxl)
library(lme4)
library(lmerTest)
# library(lmtest)
library(tidyverse)
library(reshape2)

# library(ez)

# setwd("U:/sgsalant/Lab Share/Julies Folder/Dissertation/Results/Experiment 2")

#Functions
std.err <- function(x) {
  sd(x)/sqrt(length(x))
}


#Import Dataset (number values only)
s1 <- read_xlsx("data/Summary_Exp2.xlsx", sheet = 2)
names(s1)<-make.names(names(s1),unique = TRUE) 
age.data <- read_xlsx("data/Demographics.xlsx", sheet = 2)
names(age.data)<-make.names(names(age.data),unique = TRUE) 
s1<- left_join(s1, age.data, by = "Subject.Number")

#Add NIH Uncorrected scores only
nih.score <- read_xlsx("data/Assessment Scores.xlsx")
names(nih.score)<-make.names(names(nih.score),unique = TRUE) 
nih_short <- dplyr::select(nih.score, Subject.Number, Inst, RawScore, Uncorrected.Standard.Score)

#Spread data
nih_wide <- dcast(nih_short, Subject.Number ~ Inst, value.var="Uncorrected.Standard.Score")
s1<- left_join(s1, nih_wide, by = "Subject.Number")

#Make new column with combined Score
s1$Score <- rowSums(s1[,c("Score_1", "Score_2")], na.rm=TRUE)

#Rename columns

s1$TMR <- s1$Target.to.Masker.Ratio
s1$TrialNum <- s1$Overall.Trial.Number.for.Subject
s1$Dim.Change <- s1$`NIH Toolbox Dimensional Change Card Sort Test Age 12+ v2.1`
s1$Flanker <- s1$`NIH Toolbox Flanker Inhibitory Control and Attention Test Age 12+ v2.1`
s1$WM <- s1$`NIH Toolbox List Sorting Working Memory Test Age 7+ v2.1`
s1$ListSort <- s1$`NIH Toolbox List Sorting Working Memory Test Age 7+ v2.1`


#Subset with relevant columns
s1_short <- select(s1, Subject.Number:TrialNum, Talker.ID.Response, Presentation, TMR, Condition, Response.Time, Score, Familiar.Talker, Age, Group, Sex, MoCA, LSPAN, Dim.Change, Flanker, WM, ListSort)


#Modify Condition Numbers
s1_short$Condition[s1_short$Condition==8] <- 1
s1_short$Condition[s1_short$Condition==9] <- 2
s1_short$Condition[s1_short$Condition==10] <- 3
s1_short$Condition[s1_short$Condition==18] <- 11
s1_short$Condition[s1_short$Condition==19] <- 12
s1_short$Condition[s1_short$Condition==20] <- 13
s1_short$Condition[s1_short$Condition==31] <- 31
s1_short$Condition[s1_short$Condition==32] <- 31
s1_short$Condition[s1_short$Condition==41] <- 41
s1_short$Condition[s1_short$Condition==42] <- 41

##Subset for Proposal Data
#s1_short <- s1_short %>%
 # filter(!(Subject.Number == c(427, 428)))

################## MONOTIC Conditions - Speech Understanding #####################


#Subset for Monotic
mono_mean <- s1_short %>%
  group_by(Condition, TMR, Group) %>%
  summarize(mean = mean(Score), se = std.err(Score)) %>%
  filter(!(Condition >3)) %>%
  filter(!(TMR == 3))

#Factor and Leveling
mono_mean$TMR <-as.factor(mono_mean$TMR)
mono_mean$Group <- as.factor(mono_mean$Group)
mono_mean$Condition <- as.factor(mono_mean$Condition)
mono_mean$TMR <- recode(mono_mean$TMR, "1" = "0", "2" = "-5", "4" = "Quiet")
mono_mean$TMR <- relevel(mono_mean$TMR, "-5", "0", "Quiet")
mono_mean$Group <- relevel(mono_mean$Group, "YN", "ON")

meanQ <- data.frame(Group = c("YN", "ON"), quiet = c(0.9982143, 0.97875))

#Barplot of Age/TMR for Monotic Condtions - Quiet as a reference
ggplot(mono_mean, aes(TMR, mean, fill = Condition)) +
  geom_bar(stat = "identity", color = "black",
           position = "dodge") +
  geom_hline(aes(linetype = Group, yintercept = quiet), meanQ, show.legend = FALSE) +
  scale_linetype_manual(values=c("twodash", "twodash")) +
  facet_wrap(~ Group) +
  geom_errorbar(aes(ymin = mean - (se), ymax = mean + (se),
                    group = Condition),
                position = position_dodge(width = .9),
                width = .3) +
  labs(x = "SNR") +
  scale_fill_manual(name = "Target Talker + \nSame-Ear Masker",
                    labels=c("FAM + UnF",
                             "UnF + FAM",
                             "UnF + UnF"),
                    values=c("#8856a7","#8c96c6","#edf8fb")) +
  scale_y_continuous(name = "Proportion Corrent", expand = c(0,0), 
                     breaks = c(0, .2, .4, .6, .8, 1.0)) +
  expand_limits(y = c(0, 1.1)) +
  theme_bw(base_size = 16)

ggsave("Monotic_SpeechScore.jpeg")

#######################################
# Start of 2/18/20 class

# 1: build small models to make sense to you as you go

# 2: start with a model that (at least initially) makes
#    sense as a model of the dependent variable

# Here, the DV is called "Score"

#Subset for Monotic
mono_all <- s1_short %>%
  filter(!(Condition >3)) %>%
  filter(!(TMR == 3))

#Factor and Leveling
mono_all$TMR <-as.factor(mono_all$TMR)
mono_all$Group <- as.factor(mono_all$Group)
mono_all$Condition <- as.factor(mono_all$Condition)
mono_all$TMR <- recode(mono_all$TMR, "1" = "0", "2" = "-5", "4" = "Quiet")
mono_all$TMR <- relevel(mono_all$TMR, "-5", "0", "Quiet")
mono_all$Group <- relevel(mono_all$Group, "YN", "ON")

mono_all_unf <- filter(mono_all, Condition == "3")

ggplot(mono_all_unf, aes(Group, Score)) + 
  geom_boxplot(aes(fill = TMR))

unique(mono_all_unf$Score)

ggplot(mono_all_unf, aes(Score)) + geom_histogram()

# looks pretty logistic!
mono_all_unf.glm0 <- glm(cbind(Score*4, 4 - Score*4) ~ 1,
                         data = mono_all_unf,
                         family = "binomial")
summary(mono_all_unf.glm0)                        

invlogit(0.41937)

# 3: sensible model has your experimental variables/design in it!

mono_all_unf.glm1 <- glm(cbind(Score*4, 4 - Score*4) ~ 
                           1 + Group,
                         data = mono_all_unf,
                         family = "binomial")
summary(mono_all_unf.glm1)                        

mono_all_unf.glm2 <- glm(cbind(Score*4, 4 - Score*4) ~ 
                           1 + Group * TMR,
                         data = mono_all_unf,
                         family = "binomial")
summary(mono_all_unf.glm2)                        

# 4: we have a fixed effect model, great! Now let's add
#    sensible random effects

# Some rules of thumb:
#   - "keep it maximal"
#   - don't go nuts

mono_all_unf.glmer2a <- glmer(cbind(Score*4, 4 - Score*4) ~ 
                                1 + Group * TMR +
                                (1 | Subject.Number) +
                                (1 | Overall.Trial.Number.for.Subject),
                              data = mono_all_unf,
                              family = "binomial")
summary(mono_all_unf.glmer2a) 

mono_all_unf.glmer2a <- glmer(cbind(Score*4, 4 - Score*4) ~ 
                                1 + Group * TMR +
                                (1 | Subject.Number) +
                                (1 | Overall.Trial.Number.for.Subject),
                              data = mono_all_unf,
                              family = "binomial")
summary(mono_all_unf.glmer2a)

# 4a: what random effects even make sense?

# here: Group is between subjects, so it doesn't make sense as a subject-level random slope

mono_all_unf.glmer2b <- glmer(cbind(Score*4, 4 - Score*4) ~ 
                                1 + Group * TMR +
                                scale(Overall.Trial.Number.for.Subject) +
                                (1 | Subject.Number),
                              data = mono_all_unf,
                              family = "binomial")
summary(mono_all_unf.glmer2b)


mono_all_unf.glmer2c <- glmer(cbind(Score*4, 4 - Score*4) ~ 
                                1 + Group * TMR +
                                scale(Overall.Trial.Number.for.Subject) +
                                (1 | Subject.Number) +
                                (1 | Target.Talker) +
                                (1 | Masker1.Talker),
                              data = mono_all_unf,
                              family = "binomial")
summary(mono_all_unf.glmer2c)

mono_all_unf.glmer2d <- glmer(cbind(Score*4, 4 - Score*4) ~ 
                                1 + Group * TMR +
                                scale(Overall.Trial.Number.for.Subject) +
                                (1 | Subject.Number) +
                                (1 | Target.Talker),
                              data = mono_all_unf,
                              family = "binomial")
summary(mono_all_unf.glmer2d)

AIC(mono_all_unf.glmer2c)
AIC(mono_all_unf.glmer2d)

anova(mono_all_unf.glmer2c, mono_all_unf.glmer2d)

mono_all_unf.glmer2e <- glmer(cbind(Score*4, 4 - Score*4) ~ 
                                1 + Group * TMR +
                                scale(Overall.Trial.Number.for.Subject) +
                                (1 + TMR | Subject.Number) +
                                (1 | Target.Talker) +
                                (1 | Masker1.Talker),
                              data = mono_all_unf,
                              family = "binomial")
summary(mono_all_unf.glmer2e)

mono_all_unf.glmer3c <- glmer(cbind(Score*4, 4 - Score*4) ~ 
                                1 + Group * TMR +
                                scale(Overall.Trial.Number.for.Subject) +
                                + LSPAN +
                                (1 | Subject.Number) +
                                (1 | Target.Talker) +
                                (1 | Masker1.Talker),
                              data = mono_all_unf,
                              family = "binomial",
                              glmerControl(optimizer = "bobyqa"))
summary(mono_all_unf.glmer3c)

mono_all_unf.glmer3e <- glmer(cbind(Score*4, 4 - Score*4) ~ 
                                1 + Group * TMR +
                                scale(Overall.Trial.Number.for.Subject) +
                                + LSPAN +
                                (1 + TMR | Subject.Number) +
                                (1 | Target.Talker) +
                                (1 | Masker1.Talker),
                              data = mono_all_unf,
                              family = "binomial",
                              glmerControl(optimizer = "bobyqa"))
summary(mono_all_unf.glmer3e)

mono_all <- mono_all %>% mutate(Condition = relevel(Condition, "3"))

mono_all.glmer3e <- glmer(cbind(Score*4, 4 - Score*4) ~ 
                                1 + Group * TMR +
                                scale(Overall.Trial.Number.for.Subject) +
                                LSPAN + 
                                (1 + TMR | Subject.Number) +
                                (1 | Target.Talker) +
                                (1 | Masker1.Talker),
                              data = mono_all,
                              family = "binomial",
                              glmerControl(optimizer = "bobyqa"))
summary(mono_all.glmer3e)

mono_all.glmer4e <- glmer(cbind(Score*4, 4 - Score*4) ~ 
                            1 + Group * TMR * Condition +
                            scale(Overall.Trial.Number.for.Subject) +
                            LSPAN + 
                            (1 + TMR | Subject.Number) +
                            (1 | Target.Talker) +
                            (1 | Masker1.Talker),
                          data = mono_all,
                          family = "binomial",
                          glmerControl(optimizer = "bobyqa"))
summary(mono_all.glmer4e)

mono_all.glmer4e <- glmer(cbind(Score*4, 4 - Score*4) ~ 
                            1 + Group * TMR * Condition +
                            scale(Overall.Trial.Number.for.Subject) +
                            LSPAN + 
                            (1 + TMR * Condition | Subject.Number) +
                            (1 | Target.Talker) +
                            (1 | Masker1.Talker),
                          data = mono_all,
                          family = "binomial",
                          glmerControl(optimizer = "bobyqa"))
summary(mono_all.glmer4e)
################## MONOTIC Conditions - Talker ID #####################


#Subset for Monotic
mono_mean_ID <- s1_short %>%
  group_by(Condition, TMR, Group) %>%
  summarize(mean = mean(Score), se = std.err(Score)) %>%
  filter(Condition %in% c(11, 12, 13)) %>%
  filter(!(TMR == 3))

#Factor and Leveling
mono_mean_ID$TMR <-as.factor(mono_mean_ID$TMR)
mono_mean_ID$Group <- as.factor(mono_mean_ID$Group)
mono_mean_ID$Condition <- as.factor(mono_mean_ID$Condition)

mono_mean_ID$TMR <- recode(mono_mean_ID$TMR, "1" = "0", "2" = "-5")
mono_mean_ID$TMR <- relevel(mono_mean_ID$TMR, "-5", "0")
mono_mean_ID$Group <- relevel(mono_mean_ID$Group, "YN", "ON")

meanQ_ID <- data.frame(Group = c("YN", "ON"), quiet = c(0.9892857, 0.945))

#

#Barplot of Age/TMR for Monotic Condtions - Quiet as a reference
ggplot(mono_mean_ID, aes(TMR, mean, fill = Condition)) +
  geom_bar(stat = "identity", color = "black",
           position = "dodge") +
  geom_hline(aes(linetype = Group, yintercept = quiet), meanQ_ID, show.legend = FALSE) +
  scale_linetype_manual(values=c("twodash", "twodash")) +
  facet_wrap(~ Group) +
  geom_errorbar(aes(ymin = mean - (se), ymax = mean + (se),
                    group = Condition),
                position = position_dodge(width = .9),
                width = .3) +
  labs(x = "SNR", title = "Monotic Talker Identification") +
  scale_fill_manual(name = "Target Talker + \nSame-Ear Masker",
                   labels=c("FAM + UnF",
                           "UnF + FAM",
                            "UnF + UnF"),
                   values=c("#8856a7","#8c96c6","#edf8fb")) +
  scale_y_continuous(name = "Proportion Corrent", expand = c(0,0), 
                     breaks = c(0, .2, .4, .6, .8, 1.0)) +
  expand_limits(y = c(0, 1.1)) +
  theme_bw()


ggsave("MonoticID_SpeechScore.jpeg")


################## ALL Dichotic Conditions - Speech Understanding #####################

#Subset for Dichotic
di_mean <- s1_short %>%
  group_by(Presentation, Condition, Group) %>%
  summarize(mean = mean(Score), se = std.err(Score)) %>%
  filter(Presentation == 2) %>%
  filter(Condition < 10)

#Factor and Leveling
di_mean$Group <- as.factor(di_mean$Group)
di_mean$Condition <- as.factor(di_mean$Condition)
di_mean$Group <- relevel(di_mean$Group, "YN", "ON")
di_mean$Condition <- factor(di_mean$Condition, levels = c("4", "6", "7", "5"))
di_mean$Condition <- recode(di_mean$Condition, "4" = "Familiar Target", "5" = "All Unfamiliar", "6" = "Familiar Masker - Target Ear", "7" = "Familiar Masker - Non-Target Ear")


ggplot(di_mean, aes(Condition, mean, fill = Condition)) +
  geom_bar(stat = "identity", color = "black",
           position = "dodge") +
  geom_errorbar(aes(ymin = mean - (se), ymax = mean + (se),
                    group = Condition),
                position = position_dodge(width = .9),
                width = .3) +
  facet_wrap(~Group) +
  labs(x = "Dichotic Condition") +
  scale_fill_manual(values=c("#006d2c","#b2e2e2", "#253494","#41b6c4"))+
  scale_y_continuous(name = "Proportion Corrent", expand = c(0,0), limits = 0:1.2) +
  theme_bw(base_size = 16) +  
  theme(axis.text.x=element_blank(), axis.ticks.x = element_blank(), legend.position = "none")
  


ggsave("DichoticALL_SpeechScore_OnYn_SE.jpeg")

################## FTalker Dichotic Conditions - Speech Understanding ###################

#Subset for Dichotic FT
diFT_mean <- s1_short %>%
  group_by(Presentation, Condition, Group) %>%
  summarize(mean = mean(Score), se = std.err(Score)) %>%
  filter(Presentation == 2) %>%
  filter(Condition %in% c(4, 5))

#Factor and Leveling
diFT_mean$Group <- as.factor(diFT_mean$Group)
diFT_mean$Condition <- as.factor(diFT_mean$Condition)
diFT_mean$Group <- relevel(diFT_mean$Group, "YN", "ON")
# diFT_mean$Condition <- recode(diFT_mean$Condition, "4" = "Familiar Target", "5" = "All Unfamiliar")


ggplot(diFT_mean, aes(Group, mean, fill = Condition)) +
  geom_bar(stat = "identity", color = "black",
           position = "dodge", width = .7) +
  geom_errorbar(aes(ymin = mean - (se), ymax = mean + (se),
                    group = Condition),
                position = position_dodge(width = .7),
                width = .3) +
  labs(x = "Age Group") +
  scale_fill_manual(name = "Condition",
                    labels=c("Familiar Talker","Unfamiliar Talker"),
                    values=c("#006d2c","#b2e2e2"))+
  scale_y_continuous(name = "Proportion Corrent", expand = c(0,0), limits = 0:1.2) +
  theme_bw(base_size = 16) +
  theme(legend.position = c(.98, .98),
        legend.justification = c("right", "top"),
        legend.title.align = .47)



ggsave("DichoticFT_SpeechScore_OnYn_SE.jpeg", width = 4, height = 4.5)
?ggsave

################## FMasker Dichotic Conditions - Speech Understanding ###################

#Subset for Dichotic FM
diFM_mean <- s1_short %>%
  group_by(Presentation, Condition, Group) %>%
  summarize(mean = mean(Score), se = std.err(Score)) %>%
  filter(Presentation == 2) %>%
  filter(Condition %in% c(5, 6, 7))

#Factor and Leveling
diFM_mean$Group <- as.factor(diFM_mean$Group)
diFM_mean$Condition <- as.factor(diFM_mean$Condition)
diFM_mean$Group <- relevel(diFM_mean$Group, "YN", "ON")
diFM_mean$Condition <- factor(diFM_mean$Condition, levels = c("6", "7", "5"))
# diFM_mean$Condition <- recode(diFM_mean$Condition, "5" = "All Unfamiliar", "6" = "Familiar -\nTarget Ear", "7" = "Familiar -\nNon-Target Ear")


ggplot(diFM_mean, aes(Group, mean, fill = Condition)) +
  geom_bar(stat = "identity", color = "black",
           position = "dodge", width = .7) +
  geom_errorbar(aes(ymin = mean - (se), ymax = mean + (se),
                    group = Condition),
                position = position_dodge(width = .7),
                width = .3) +
  labs(x = "Age Group") +
  scale_fill_manual(name = "Condition",
                    labels=c("Familiar - Target Ear",
                             "Familiar - Opposite Ear", "Both Unfamiliar"),
                    values=c("#253494","#41b6c4","#b2e2e2")) +
  scale_y_continuous(name = "Proportion Corrent", expand = c(0,0), limits = 0:1.2) +
  theme_bw(base_size = 16) +
  theme(legend.position = c(.98, .98),
        legend.justification = c("right", "top"),
        legend.title.align = .47)


ggsave("DichoticFM_SpeechScore_OnYn_SE.jpeg", width = 5, height = 5.3)
?ggsave


#################################################

################## FTalker Dichotic Conditions - Talker ID ###################

#Subset for Dichotic FT
diFT_mean_ID <- s1_short %>%
  group_by(Presentation, Condition, Group) %>%
  summarize(mean = mean(Score), se = std.err(Score)) %>%
  filter(Presentation == 2) %>%
  filter(Condition %in% c(14, 15, 41))

#Factor and Leveling
diFT_mean_ID$Group <- as.factor(diFT_mean_ID$Group)
diFT_mean_ID$Condition <- as.factor(diFT_mean_ID$Condition)
diFT_mean_ID$Group <- relevel(diFT_mean_ID$Group, "YN", "ON")
# diFT_mean_ID$Condition <- recode(diFT_mean_ID$Condition, "4" = "Familiar Target", "5" = "All Unfamiliar")


ggplot(diFT_mean_ID, aes(Group, mean, fill = Condition)) +
  geom_bar(stat = "identity", color = "black",
           position = "dodge", width = .7) +
  geom_errorbar(aes(ymin = mean - (se), ymax = mean + (se),
                    group = Condition),
                position = position_dodge(width = .7),
                width = .3) +
  labs(x = "Age Group") +
  scale_fill_manual(name = "Condition",
                    labels=c("Familiar Target","All Unfamiliar"),
                    values=c("#006d2c","#b2e2e2"))+
  scale_y_continuous(name = "Proportion Corrent", expand = c(0,0), 
                     breaks = c(0, .2, .4, .6, .8, 1.0)) +
  expand_limits(y = c(0, 1.1)) +
  theme_bw(base_size = 16)



ggsave("DichoticFT_ID_SE.jpeg", width = 6.4, height = 5)

################## FMasker Dichotic Conditions - Talker ID ###################

#Subset for Dichotic FM
diFM_mean_ID <- s1_short %>%
  group_by(Presentation, Condition, Group) %>%
  summarize(mean = mean(Score), se = std.err(Score)) %>%
  filter(Presentation == 2) %>%
  filter(Condition %in% c(15, 16, 17))

#Factor and Leveling
diFM_mean_ID$Group <- as.factor(diFM_mean_ID$Group)
diFM_mean_ID$Condition <- as.factor(diFM_mean_ID$Condition)
diFM_mean_ID$Group <- relevel(diFM_mean_ID$Group, "YN", "ON")
diFM_mean_ID$Condition <- factor(diFM_mean_ID$Condition, levels = c("16", "17", "15"))
# diFM_mean_ID$Condition <- recode(diFM_mean_ID$Condition, "5" = "All Unfamiliar", "6" = "Familiar -\nTarget Ear", "7" = "Familiar -\nNon-Target Ear")


ggplot(diFM_mean_ID, aes(Group, mean, fill = Condition)) +
  geom_bar(stat = "identity", color = "black",
           position = "dodge", width = .7) +
  geom_errorbar(aes(ymin = mean - (se), ymax = mean + (se),
                    group = Condition),
                position = position_dodge(width = .7),
                width = .3) +
  labs(x = "Age Group") +
  scale_fill_manual(name = "Condition",
                    labels=c("Familiar - Target Ear",
                             "Familiar - Opposite Ear", "Both Unfamiliar"),
                    values=c("#253494","#41b6c4","#b2e2e2")) +
  scale_y_continuous(name = "Proportion Corrent", expand = c(0,0), 
                     breaks = c(0, .2, .4, .6, .8, 1.0)) +
  expand_limits(y = c(0, 1.1)) +
  theme_bw(base_size = 16)



ggsave("DichoticFM_ID_SE.jpeg")



################################################################

#Statistics #

############################################################


#Monotic Speech Understanding

#ANOVA
mono <- s1_short %>%
  group_by(Subject.Number, Condition, TMR, Group) %>%
  summarize(mean = mean(Score), se = std.err(Score)) %>%
  filter(!(Condition >3)) %>%
  filter(!(TMR == 3))

mono$TMR <-as.factor(mono$TMR)
mono$Group <- as.factor(mono$Group)
mono$Condition <- as.factor(mono$Condition)
mono$Subject.Number <-as.factor(mono$Subject.Number)
mono$TMR <- recode(mono$TMR, "1" = "0", "2" = "-5", "3" = "+5")
mono$TMR <- relevel(mono$TMR, "0", "-5")
mono$Group <- relevel(mono$Group, "YN", "ON")

mono$Condition <- relevel(mono$Condition, ref = "3")


mono1.aov <- aov(mean ~ Condition * TMR * Group, data = mono)
summary(mono1.aov)
anova(mono1.aov)

TukeyHSD(mono1.aov, "Condition", ordered = TRUE)
TukeyHSD(mono1.aov, "TMR", ordered = TRUE)
TukeyHSD(mono1.aov, "Group", ordered = TRUE)

#Dichotic speech understanding - all conditions

di_su <- s1_short %>%
  group_by(Subject.Number, Condition, TMR, Group) %>%
  summarize(mean = mean(Score), se = std.err(Score)) %>%
  filter(Condition >3 & Condition < 8)

di_su$Group <- as.factor(di_su$Group)
di_su$Condition <- as.factor(di_su$Condition)
di_su$Subject.Number <-as.factor(di_su$Subject.Number)
di_su$Group <- relevel(di_su$Group, "YN", "ON")

di1.aov <- aov(mean ~ Condition * Group, data = di_su)
summary(di1.aov)
anova(di1.aov)

TukeyHSD(di1.aov, "Condition", ordered = TRUE)
TukeyHSD(di1.aov, "Group", ordered = TRUE)

#Dichotic conditions - effect of familiar masker


di_su_m <- s1_short %>%
  group_by(Subject.Number, Condition, TMR, Group) %>%
  summarize(mean = mean(Score), se = std.err(Score)) %>%
  filter(Condition >4 & Condition < 8)

di_su_m$Group <- as.factor(di_su_m$Group)
di_su_m$Condition <- as.factor(di_su_m$Condition)
di_su_m$Subject.Number <-as.factor(di_su_m$Subject.Number)
di_su_m$Group <- relevel(di_su_m$Group, "YN", "ON")

di2.aov <- aov(mean ~ Condition * Group, data = di_su_m)
summary(di2.aov)
anova(di2.aov)

TukeyHSD(di2.aov, "Condition", ordered = TRUE)
TukeyHSD(di2.aov, "Group", ordered = TRUE)

############################################

#Monotic Talker ID

#ANOVA
mono_ID <- s1_short %>%
  group_by(Subject.Number, Condition, TMR, Group) %>%
  summarize(mean = mean(Score), se = std.err(Score)) %>%
  filter(Condition > 10 & Condition < 14) %>%
  filter(!(TMR == 3))

mono_ID$TMR <-as.factor(mono_ID$TMR)
mono_ID$Group <- as.factor(mono_ID$Group)
mono_ID$Condition <- as.factor(mono_ID$Condition)
mono_ID$Subject.Number <-as.factor(mono_ID$Subject.Number)
mono_ID$Group <- relevel(mono_ID$Group, "YN", "ON")

mono_ID1.aov <- aov(mean ~ Condition * TMR * Group, data = mono_ID)
summary(mono_ID1.aov)
anova(mono_ID1.aov)

TukeyHSD(mono_ID1.aov, "Condition", ordered = TRUE)
TukeyHSD(mono_ID1.aov, "Group", ordered = TRUE)

TukeyHSD(mono_ID1.aov, "Condition:Group", ordered = TRUE)

#Dichotic Talker ID

di_ID <- s1_short %>%
  group_by(Subject.Number, Condition, Group) %>%
  summarize(mean = mean(Score), se = std.err(Score)) %>%
  filter(Condition >13 & Condition < 18)

di_ID$Group <- as.factor(di_ID$Group)
di_ID$Condition <- as.factor(di_ID$Condition)
di_ID$Subject.Number <-as.factor(di_ID$Subject.Number)
di_ID$Group <- relevel(di_ID$Group, "YN", "ON")

di_ID1.aov <- aov(mean ~ Condition * Group, data = di_ID)
summary(di_ID1.aov)
anova(di_ID1.aov)

TukeyHSD(di_ID1.aov, "Condition", ordered = TRUE)
TukeyHSD(di_ID1.aov, "Group", ordered = TRUE)

TukeyHSD(di_ID1.aov, "Condition:Group", ordered = TRUE)

#Dichotic conditions - effect of familiar masker


di_su_m <- s1_short %>%
  group_by(Subject.Number, Condition, TMR, Group) %>%
  summarize(mean = mean(Score), se = std.err(Score)) %>%
  filter(Condition >4 & Condition < 8)

di_su_m$Group <- as.factor(di_su_m$Group)
di_su_m$Condition <- as.factor(di_su_m$Condition)
di_su_m$Subject.Number <-as.factor(di_su_m$Subject.Number)
di_su_m$Group <- relevel(di_su_m$Group, "YN", "ON")

di2.aov <- aov(mean ~ Condition * Group, data = di_su_m)
summary(di2.aov)
anova(di2.aov)

TukeyHSD(di2.aov, "Condition", ordered = TRUE)
TukeyHSD(di2.aov, "Group", ordered = TRUE)


#################### Monotic vs. Dichotic Conditions  ################################


#Subset with relevant columns
s2_short <- select(s1, Subject.Number:Overall.Trial.Number.for.Subject, Talker.ID.Response, Presentation, Target.to.Masker.Ratio, Condition, Response.Time, Score, Familiar.Talker, Age, Group, Sex, MoCA, LSPAN)

#Modify Condition Numbers
s2_short$Cond_Merge <- s2_short$Condition

s2_short$Cond_Merge[s2_short$Cond_Merge == 4] <- 1
s2_short$Cond_Merge[s2_short$Cond_Merge == 6] <- 2
s2_short$Cond_Merge[s2_short$Cond_Merge == 5] <- 3


#Rename columns
s2_short$TMR <- s2_short$Target.to.Masker.Ratio
s2_short$TrialNum <- s2_short$Overall.Trial.Number.for.Subject

#Subset
m_d <- s2_short %>%
  group_by(Presentation, TMR, Cond_Merge, Condition, Group) %>%
  summarize(mean = mean(Score), se = std.err(Score)) %>%
  filter(!(Condition >6)) %>%
  filter((TMR == 1))

#Factor and Leveling
m_d$TMR <-as.factor(m_d$TMR)
m_d$Group <- as.factor(m_d$Group)
m_d$Condition <- as.factor(m_d$Condition)
m_d$Presentation <- as.factor(m_d$Presentation)
m_d$Cond_Merge <- as.factor(m_d$Cond_Merge)
m_d$TMR <- recode(m_d$TMR, "1" = "0", "2" = "-5", "4" = "Quiet")
m_d$Cond_Merge <- recode(m_d$Cond_Merge, "1" = "Familiar Target", 
                         "2" = "Familiar Masker", "3" = "Both Unfamiliar")
m_d$Presentation <- recode(m_d$Presentation, "1" = "Monotic", "2" = "Dichotic")
m_d$Group <- relevel(m_d$Group, "YN", "ON")

#Barplot of Age/Mono-Dichotic _ Condition on x axis
ggplot(m_d, aes(Presentation, mean)) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.3, color = "black") +
  geom_point(aes(color = Cond_Merge), size = 4.5) +
  geom_point(shape = 1, size = 4.5, color = "black") +
  geom_line(aes(group = Cond_Merge)) +
  facet_wrap(~Group) +
  scale_y_continuous(name = "Proportion Corrent", expand = c(0,0), 
                     breaks = c(0, .2, .4, .6, .8, 1.0)) +
  scale_colour_manual(values=c("#8856a7","#8c96c6","#edf8fb")) +
  expand_limits(y = c(0, 1.1)) +
  theme_bw(base_size = 16) 
  #theme(legend.position = "none")

ggsave("MonoDichotic_SpeechScore.jpeg")

m_d2 <- m_d
m_d2$TMR <-as.factor(m_d2$TMR)
m_d2$Group <- as.factor(m_d2$Group)
m_d2$Condition <- as.factor(m_d2$Condition)
m_d2$Presentation <- as.factor(m_d2$Presentation)
m_d2$Cond_Merge <- as.factor(m_d2$Cond_Merge)
m_d2$TMR <- recode(m_d2$TMR, "1" = "0", "2" = "-5", "4" = "Quiet")
m_d2$Cond_Merge <- recode(m_d2$Cond_Merge, "1" = "Familiar\nTarget", 
                         "2" = "Familiar\nMasker", "3" = "Both\nUnfamiliar")
m_d2$Presentation <- recode(m_d2$Presentation, "1" = "Monotic", "2" = "Dichotic")
m_d2$Group <- relevel(m_d2$Group, "YN", "ON")


ggplot(m_d2, aes(Cond_Merge, mean)) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.3, color = "black") +
  geom_point(aes(color = Presentation), size = 4.5) +
  geom_point(shape = 1, size = 4.5, color = "black") +
  geom_line(aes(group = Presentation)) +
  facet_wrap(~Group) +
  scale_y_continuous(name = "Proportion Corrent", expand = c(0,0), 
                     breaks = c(0, .2, .4, .6, .8, 1.0)) +
  scale_colour_manual(name = "Presentation",
                      labels=c("Monotic",
                               "Dichotic"),
                      values=c("blue","green")) +
  labs(x = "Condition") + 
  expand_limits(y = c(0, 1.1)) +
  theme_bw()

ggsave("MonoDichotic_Flip_SpeechScore.jpeg")


