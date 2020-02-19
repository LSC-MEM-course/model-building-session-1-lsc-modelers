
library(dplyr)
library(tidyverse)
library(ggplot2)
library(openxlsx)
library(lmerTest)
library(emmeans)

library(MuMIn)
library(AICcmodavg)

# load adaptation data
# setwd("U:\\sgsalant\\Lab Share\\RB folder\\Dissertation\\Data") #set working directory

adapt = read.csv('data/adapt.csv')

adapt$PresentationOrder = factor(adapt$PresentationOrder, levels = c("First", "Second", "Third", "Fourth", "Fifth"))

summary(adapt)

#plotting, for sanity
# ggplot(adapt, aes(Block, PercentCorrect, colour = Condition)) + geom_smooth(method = 'loess') + facet_wrap(~Group2)

condition_names = as_labeller(c(`C1` = "C1: Single ENG", `C2` = "C2: Single SPA", `C3` = "C3: Single JAP", `C4` = "C4: Multiple SPA", `C5` = "C5: Multiple L1s"))

#Plot by condition
ggplot(adapt, aes(Block, PercentCorrect,
                  #fill = Group2,
                  colour = Group2)) +
  geom_smooth(method = 'loess', size=1.1) +
  facet_wrap(.~Condition, nrow=1, labeller = condition_names) +
  xlab("Sentence number") +
  theme_bw() +
  theme(#panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank(),
    #panel.background = element_blank(),
    legend.title = element_text(),
    legend.position = c(0.07, 0.15),
    legend.background = element_rect(fill="white",
                                     size=0.5, linetype="solid",
                                     colour ="black"),
    plot.title = element_text(color="black", size=14, face="bold.italic"),
    axis.title.x = element_text(color="black", size=12, face="bold"),
    axis.title.y = element_text(color="black", size=12, face="bold"),
    axis.text.x = element_text(color="black", size=12),
    axis.text.y = element_text(color="black", size=12),
    strip.text = element_text(size = 11, face = "bold")) +
  labs(fill = "Listener Group", colour = "Listener Group") +
  scale_colour_brewer(palette = "Dark2") +
  #scale_fill_brewer(palette = "Dark2") +
  labs(title = "Rapid adaptation", y = "Proportion Correct", x = "Sentence number")


#Plot by presentation order
# ggplot(adapt, aes(Block, PercentCorrect,
#                   #fill = Group2,
#                   colour = Group2)) +
#   geom_smooth(method = 'loess', size=1.1) +
#   facet_wrap(.~PresentationOrder, nrow=1) +
#   xlab("Sentence number") +
#   theme_bw() +
#   theme(#panel.grid.major = element_blank(),
#     #panel.grid.minor = element_blank(),
#     #panel.background = element_blank(),
#     legend.title = element_text(),
#     legend.position = c(0.07, 0.15),
#     legend.background = element_rect(fill="white",
#                                      size=0.5, linetype="solid",
#                                      colour ="black"),
#     plot.title = element_text(color="black", size=14, face="bold.italic"),
#     axis.title.x = element_text(color="black", size=12, face="bold"),
#     axis.title.y = element_text(color="black", size=12, face="bold"),
#     axis.text.x = element_text(color="black", size=12),
#     axis.text.y = element_text(color="black", size=12),
#     strip.text = element_text(size = 11, face = "bold")) +
#   labs(fill = "Listener Group", colour = "Listener Group") +
#   scale_colour_brewer(palette = "Dark2") +
#   #scale_fill_brewer(palette = "Dark2") +
#   labs(title = "Rapid adaptation, by presentation order", y = "Proportion Correct", x = "Sentence number")


######creating polynomial time terms
t = poly(unique(adapt$Block), 4)
adapt[, paste("ot", 1:4, sep="")] = t[adapt$Block, 1:4]


####GCM Models - DO NOT RUN THESE! Just here to see the structure. Load in the objects instead###########

# model.1217 = glmer(PercentCorrect ~ (ot1+ot2+ot3)*Group*Condition
#                    + (ot1+ot2+ot3|Subject)
#                    + (ot1+ot2|Subject:Condition),
#                    #+ (Subject|Sentence),
#                    data = adapt, family = binomial, weights = Mean.PossibleKeywords,
#                    glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun=20000)), REML = F)
#
# model.1217a = glmer(PercentCorrect ~ (ot1+ot2+ot3)*Group*Condition + PresentationOrder*(ot1+ot2+ot3)
#                     + (ot1+ot2+ot3|Subject)
#                     + (ot1+ot2|Subject:Condition),
#                     #+ (Subject|Sentence),
#                     data = adapt, family = binomial, weights = Mean.PossibleKeywords,
#                     glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun=25000)), REML = F)
#
# model.1217a2 = glmer(PercentCorrect ~ (ot1+ot2+ot3)*Group*Condition + PresentationOrder*(ot1+ot2+ot3)
#                      + PresentationOrder*Condition
#                      + (ot1+ot2+ot3|Subject)
#                      + (ot1+ot2|Subject:Condition),
#                      #+ (Subject|Sentence),
#                      data = adapt, family = binomial, weights = Mean.PossibleKeywords,
#                      glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun=50000)), REML = F)
ggplot(adapt, aes(PercentCorrect)) + geom_histogram()

windows()
ggplot(adapt, aes(x = Block)) +
    geom_line(aes(y = ot1), linetype = 1) +
    geom_line(aes(y = ot2), linetype = 2) +
    geom_line(aes(y = ot3), linetype = 3)



model.1 = lmer(PercentCorrect ~ (ot1+ot2+ot3) + Group + Condition +
                   (ot1+ot2+ot3|Subject),
                   #+ (Subject|Sentence),
                   data = adapt, #family = binomial, # weights = Mean.PossibleKeywords,
                   REML = FALSE)
                   #glmerControl(optimizer = "bobyqa"), optCtrl = list(maxfun=20000)), REML = F)

summary(model.1)

adapt %>% group_by(Subject, Condition) %>% summarize(N = n()) %>%
    summarize(N.cond = n()) %>% summary()

model.1a = lmer(PercentCorrect ~ (ot1+ot2+ot3) + Condition + Group +
                   (ot1+ot2+ot3 + Condition| Subject),
                   #+ (Subject|Sentence),
                   data = adapt, #family = binomial, # weights = Mean.PossibleKeywords,
                   REML = FALSE)
summary(model.1a)

anova(model.1, model.1a)

model.1 = lmer(PercentCorrect ~ (ot1+ot2+ot3) + Group + Condition +
                   (ot1+ot2+ot3|Subject),
                   #+ (Subject|Sentence),
                   data = adapt, #family = binomial, # weights = Mean.PossibleKeywords,
                   REML = FALSE)
                   #glmerControl(optimizer = "bobyqa"), optCtrl = list(maxfun=20000)), REML = F)

summary(model.1)

glm_fit1 <- glmer(cbind(Mean.Response.RESP, Mean.PossibleKeywords - Mean.Response.RESP) ~
                      (ot1 + ot2 + ot3) * Condition + (1 + ot1 + ot2 + ot3|Subject),
                  data = adapt,
                  family = "binomial",
                  control = glmerControl(optimizer = "bobyqa"))
summary(glm_fit1, corr = FALSE)

load('saved_models/model1217.RData')
load('saved_models/model1217a.RData')
load('saved_models/model1217a2.RData')

anova(model.1217, model.1217a, model.1217a2) #model.1217a is preferable

model.sel(model.1217, model.1217a, model.1217a2)
aictab(list(model.1217, model.1217a, model.1217a2))

summary(model.1217a)
summary(model.1217)
summary(model.1217a2)

######QUESTION 1: How to interpret this? Where to go from here? ##########



##########################Part 2: Using slopes as DV##########################
adapt.ran = as.data.frame(ranef(model.1217a, condVar = TRUE, whichel = "Subject:Condition"))

summary(adapt.ran)
head(adapt.ran)


adapt.ran = separate(adapt.ran, grp, into = c("Subject", "Condition"), sep = ":", extra = "merge")

#add group column
adapt.ran$Subject = as.numeric(adapt.ran$Subject)
adapt.ran$group = factor(ifelse(adapt.ran$Subject<400, "OHI", ifelse(adapt.ran$Subject < 7200, "YNH", "ONH")), levels = c("YNH", "ONH", "OHI"))

#####load in and clean audcog data (to use as predictors for random slopes)

#read in
audcog = read.xlsx('Auditory and Cognitive Info_RBP1v2.xlsx', 1, cols = 1:39)
summary(audcog)

#clean out old subjects
audcog$version = as.factor(ifelse(audcog$Subject<300, "old", ifelse(audcog$Subject>1000, "new", ifelse(audcog$Subject > 300 & audcog$Subject < 900, "new", "old"))))

audcog = audcog[!(audcog$version == "old"), ]
audcog=audcog[!(audcog$Subject == "7202"), ]
audcog=audcog[!(audcog$Subject == "7120"), ]

audcog$Subject = as.factor(audcog$Subject)
audcog$group = as.factor(audcog$Group)

audcog.min = audcog[ , c("Subject","Age", "Sex", "MOCA", "TrailA", "TrailB", "TrailABA", "Stroop-C", "Stroop-I", "Stroop-N", "Stroop-I-C", "ListSortingStandard", "ListSortingAgeCorr", "VocabStandard", "VocabAgeCorr", "CardSortStandard", "CardSortAgeCorr", "PTA", "HFPTA")]

summary(audcog.min)

##merge audcog data with random slopes
ranef.audcog = merge(adapt.ran, audcog.min, by = "Subject")
ranef.audcog$Condition = factor(ranef.audcog$Condition)
ranef.audcog$Sex = factor(ranef.audcog$Sex)
ranef.audcog$Subject = factor(ranef.audcog$Subject)

ranef.audcog$StroopEffect = ranef.audcog$`Stroop-I-C`

summary(ranef.audcog)

##########z-transforming predictors of interest
ranef.audcog$StroopEffect_z = scale(ranef.audcog$StroopEffect)
ranef.audcog$VocabStandard_z = scale(ranef.audcog$VocabStandard)
ranef.audcog$TrailABA_z = scale(ranef.audcog$TrailABA)
ranef.audcog$ListSortingStandard_z = scale(ranef.audcog$ListSortingStandard)


#Modelling
linearmod.0 = lmer(condval ~ term + (1|Subject), data = ranef.audcog)

linearmod.1 = lmer(condval ~ Condition + term + (1|Subject), data = ranef.audcog)

linearmod.2 = lmer(condval ~ Condition + group + term + (1|Subject), data = ranef.audcog)

############QUESTION 2: What to do about singular fits? How to do this correctly? ##########




##############################Part 3: other data to look at for simpler modelling###########################



#creating a table where I can calculate magnitude of adaptation
adapt$block.fifths = as.factor(ifelse(adapt$Block <=5, "Block1", ifelse(adapt$Block >5 & adapt$Block <= 10, "Block2", ifelse(adapt$Block >10 & adapt$Block <= 15, "Block3", ifelse(adapt$Block >15 & adapt$Block <= 20, "Block4", ifelse(adapt$Block >20 & adapt$Block <= 25, "Block5", "Block6" ))))))

adapt.mag = as.data.frame(summarize(group_by(adapt, Condition, Subject, Group2, block.fifths), mean = mean(PercentCorrect, na.rm = TRUE)))

summary(adapt.mag)

adapt.mag_wide = pivot_wider(adapt.mag, names_from = block.fifths, values_from = mean)

summary(adapt.mag_wide)

#replacing 0 values with minimum value (0.05)
adapt.mag_wide$Block1.adj = ifelse(adapt.mag_wide$Block1 == 0, 0.05, adapt.mag_wide$Block1)

adapt.mag_wide$RelativeChange = (adapt.mag_wide$Block6 - adapt.mag_wide$Block1.adj)/adapt.mag_wide$Block1.adj


##plotting magnitude of change
ggplot(adapt.mag_wide, aes(Condition, RelativeChange)) + geom_boxplot(aes(fill = Condition)) + facet_grid(~Group2) + geom_hline(yintercept=0, color = "black")

##adding in the auditory/cognitive data
mag.with.cog = merge(adapt.mag_wide, audcog.min, by = "Subject")
mag.with.cog$StroopEffect = mag.with.cog$`Stroop-I-C`

summary(mag.with.cog)

##########z-transforming predictors of interest
mag.with.cog$StroopEffect_z = scale(mag.with.cog$StroopEffect)
mag.with.cog$VocabStandard_z = scale(mag.with.cog$VocabStandard)
mag.with.cog$TrailABA_z = scale(mag.with.cog$TrailABA)
mag.with.cog$ListSortingStandard_z = scale(mag.with.cog$ListSortingStandard)


#modelling
######QUESTION 3: How do you decide what to leave in and take out? p-values or model comparisons? etc.#####

magcog.0 = lmer(RelativeChange ~ Condition + (1|Subject), data = mag.with.cog)
magcog.1 = lmer(RelativeChange ~ Condition*Group2 + (1|Subject), data = mag.with.cog)
magcog.2 = lmer(RelativeChange ~ Condition*Group2 + Condition*StroopEffect_z + (1|Subject), data = mag.with.cog)
magcog.3 = lmer(RelativeChange ~ Condition*Group2 + Condition*StroopEffect_z + Condition*VocabStandard_z + (1|Subject), data = mag.with.cog)
magcog.4 = lmer(RelativeChange ~ Condition*Group2 + Condition*StroopEffect_z + Condition*VocabStandard_z + Condition*TrailABA_z + (1|Subject), data = mag.with.cog)

magcog.5 = lmer(RelativeChange ~ Condition*StroopEffect_z + Condition*VocabStandard_z + Condition*TrailABA_z + (1|Subject), data = mag.with.cog)

magcog.6 = lmer(RelativeChange ~ Condition*StroopEffect_z + (1|Subject), data = mag.with.cog)

#...you could do this forever.
summary(magcog.6)
