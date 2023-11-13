# Scratch assay data analysis
# The effect of treatment (Group) and treatment protein concentration (CTHRC1) on cell migration 
# measured as the size of cell gap (Scratch_FOV) over time (Day 0, 1, 2)

#
library(dplyr)
library(tidyr)
library(ggplot2)
library(lme4)
library(lmerTest)
library(RColorBrewer)
library(ciTools)
##


# load dataset
dat1 <- read.csv("C:/Users/aggiej/Documents/Fellowship/Tissue culture/MSC_CM_in_vitro/Scratch_assay/Scratch_results_merged.csv") %>% 
  arrange(Recipient_ID)


head(dat1)
colnames(dat1)

# convert Group variable to factor
dat1$Group <- as.factor(dat1$Group)

# create baseline variable = size of scratch on Day 0 (before treatment)
dat_baseline <- dat1 %>%
  arrange(Recipient_ID) %>%
  dplyr::select(Recipient_ID, Well, Day, Scratch_FOV) %>%
  filter(Day == 0) %>%
  dplyr::select(-Day) %>%
  rename(Scratch_FOV_baseline = Scratch_FOV)

# add baseline to the dataset
dat2 <- dat1 %>%
  dplyr::select(-Scratch_area, -FOV_area, -X) %>%
  left_join(dat_baseline, by = join_by(Well, Recipient_ID)) %>%
  filter(Day!=0) %>%
  mutate(Day = as.factor(Day))

## MIXED-EFFECT MODEL
dat2 <- dat2 %>% mutate(CTHRC2 = (CTHRC1/sd(CTHRC1))) # scale CTHRC1
dat3 <- dat2 %>% mutate(CTHRC2 = ifelse(Group == 0, 0 , CTHRC2)) # code CTHRC2 as interaction term with Group

# visualise distribution of the response variable
dat2 %>% ggplot(aes(x = Scratch_FOV, fill = Day)) + geom_histogram(position = 'identity', alpha = 0.5, bins = 15)

# MODEL 1
m1 <- lmer(Scratch_FOV ~ Scratch_FOV_baseline + Day + Recipient_age*CTHRC2 + (1 | Recipient_ID),
          data = dat3)
summary(m1)

# save model coefficients to a data frame
coefs <- as.data.frame(coef(summary(m1)))
write.csv(coefs, 'Scratch_M1_results.csv')

# visualisation of data model prediction
dat4 <- dat3 %>% mutate(Scratch_FOV_baseline = mean(Scratch_FOV_baseline)) %>% # fix value of predictors not shown on graph
  mutate(Recipient_age = mean(Recipient_age))

dat5 <- dat4 %>% mutate(fit = predict(m1, newdata = dat4, re.form = NA))

# alternative prediction with confidence intervals
dat5a <- add_ci(df = dat4,
                fit = m1,
                type = 'parametric',
                includeRanef = F,
                names = c('LCB', 'UCB'),
                alpha = 0.05)

# plot data and predicted trendline
v1a <- dat5a %>% mutate(CTHRC2a = CTHRC1/sd(CTHRC1)) %>%
  filter(Day == 2) %>%
  ggplot(aes(x = CTHRC2a, y = Scratch_FOV)) +
  geom_jitter(aes(x = CTHRC2a, colour = Group), size = 1) + 
  geom_line(aes(x = CTHRC2a, y = pred, colour = Group)) +
  geom_ribbon(aes(x = CTHRC2a, ymin = LCB, ymax = UCB, group = Group), alpha = 0.2) +
  theme_classic(base_size = 9) +
  ggtitle('Effect of CTHRC1 concentration in MSC CM \n on tenocyte migration') +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab('Scratch area / FOV') +
  xlab('CTHRC1 [sd]') +
  scale_colour_manual(name = 'Treatment', labels = c('Control', 'CM'), values = palette('Dark2'))

v1a

ggsave('Fig1_scratch.tiff', v1a, units = 'mm', width = 100, height = 75)

# summarise scratch size by day/group
dat3sum <- dat3 %>% 
  mutate(CTHRC2a = CTHRC1/sd(CTHRC1)) %>%
  filter(Day == 2) %>% 
  group_by(Group, CTHRC2a) %>% summarise(mean = mean(Scratch_FOV))

dat3sum

dat3sum[10,3]/dat3sum[1,3] # mean ratio Day 1
dat3sum[18,3]/dat3sum[9,3] # mean ratio Day 2




