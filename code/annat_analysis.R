library(MASS)
library(tidyverse)
library(lme4)
library(lmerTest)
library(car)
library(arm)
library(readxl)
library(fastDummies)

data_avg <- read_excel("Averages.xlsx", sheet = 1)
data_trials <- read_excel("Trial by Trial (1).xlsx", sheet = 1)

# data_trials %>% group_by(Pt) %>% summarize(PC = mean(`AO PC`),
#                                            mean_ACC = mean(ST_ACC, na.rm = TRUE))

data_trials <- data_trials %>%
    filter(!is.na(Keywords)) %>%
    mutate(n_responses = 5,
           Condition_cat = paste0("Condition_", Condition)) %>%
    dummy_columns("Condition")

head(data_trials) %>% as.data.frame()

cell_means <- data_trials %>%
    group_by(Condition_cat, Group) %>% summarize(mean_PC = mean(PC))

ggplot(data_trials, aes(Condition_cat, PC)) + geom_boxplot() + facet_wrap(~ Group)
windows()
ggplot(cell_means, aes(Condition_cat, mean_PC)) + geom_bar(stat = "identity") + facet_wrap(~ Group)
ggplot(data_trials, aes(PC)) + geom_bar() +
    facet_grid(Group ~ Condition_cat)

glm_fit1 <- glm(cbind(Keywords, n_responses - Keywords) ~
                    Condition_cat,
                data = data_trials,
                family = "binomial")
summary(glm_fit1)

summary(glm_fit1)$coefficients["(Intercept)", "Estimate"]
invlogit(summary(glm_fit1)$coefficients["(Intercept)", "Estimate"])
data_trials %>% filter(!is.na(Keywords)) %>%
    group_by(Condition_cat) %>% summarize(mean_ACC = mean(PC))

glmer_fit1a <- glmer(cbind(Keywords, n_responses - Keywords) ~
                         Condition_cat + (1|Pt),
                     data = data_trials,
                     family = "binomial")
summary(glmer_fit1a, corr = FALSE)

head(ranef(glmer_fit1a)$Pt)
head(ranef(glmer_fit1a)$Pt$`(Intercept)`)

pt_avgs <- data_trials %>%
    group_by(Pt) %>% summarize(mean_ACC_bypt = mean(PC))

pt_avgs <- pt_avgs %>%
    mutate(blups = invlogit(ranef(glmer_fit1a)$Pt$`(Intercept)`)) %>%
    pivot_longer(mean_ACC_bypt:blups, names_to = "type", values_to = "values")

ggplot(pt_avgs, aes(values)) + geom_density(aes(fill = type), alpha = .2)

glmer_fit1b <- glmer(cbind(Keywords, n_responses - Keywords) ~
                         Condition_cat + (1|Pt) + (1|Item),
                     data = data_trials,
                     family = "binomial")
summary(glmer_fit1b, corr = FALSE)

glmer_fit2b <- glmer(cbind(Keywords, n_responses - Keywords) ~
                         Condition_cat + Group + (1|Pt) + (1|Item),
                     data = data_trials,
                     family = "binomial")
summary(glmer_fit2b, corr = FALSE)

glmer_fit3b <- glmer(cbind(Keywords, n_responses - Keywords) ~
                         Condition_cat * Group + (1|Pt) + (1|Item),
                     data = data_trials,
                     family = "binomial")
summary(glmer_fit3b, corr = FALSE)

anova(glmer_fit2b, glmer_fit3b)

glmer_fit3c <- glmer(cbind(Keywords, n_responses - Keywords) ~
                         Condition_cat * Group + (1 + Condition_cat|Pt) + (1|Item),
                     data = data_trials,
                     family = "binomial")
summary(glmer_fit3c, corr = FALSE)

anova(glmer_fit3b, glmer_fit3c)

glmer_fit3d <- glmer(cbind(Keywords, n_responses - Keywords) ~
                         (Condition_2 + Condition_3 + Condition_4 +
                          Condition_5 + Condition_6) * Group +
                          (1 + Condition_2 + Condition_3 + Condition_4 +
                          Condition_5 + Condition_6|Pt) + (1|Item),
                     data = data_trials,
                     family = "binomial")

summary(glmer_fit3d, corr = FALSE)

anova(glmer_fit3c, glmer_fit3d) # identical

glmer_fit3e <- glmer(cbind(Keywords, n_responses - Keywords) ~
                         (Condition_2 + Condition_3 + Condition_4 +
                          Condition_5 + Condition_6) * Group +
                          (1 + Condition_2 + Condition_3 + Condition_4 +
                          Condition_5 + Condition_6||Pt) + (1|Item),
                     data = data_trials,
                     family = "binomial")
summary(glmer_fit3e, corr = FALSE)

anova(glmer_fit3d, glmer_fit3e)

allFit(glmer_fit3d)

glmer_fit3d <- glmer(cbind(Keywords, n_responses - Keywords) ~
                         (Condition_2 + Condition_3 + Condition_4 +
                          Condition_5 + Condition_6) * Group +
                         (1 + Condition_2 + Condition_3 + Condition_4 +
                          Condition_5 + Condition_6|Pt) + (1|Item),
                     data = data_trials,
                     family = "binomial",
                     control = glmerControl(optimizer = "bobyqa"))

summary(glmer_fit3d, corr = FALSE)
pca_fit3d <- rePCA(glmer_fit3d)
summary(pca_fit3d)
screeplot(pca_fit3d$Pt, type = "lines")
