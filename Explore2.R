# ---- Clear Env ----
rm(list=ls())

# ---- Packages ----
library(reshape2)
library(lme4)
library(sjPlot)
library(glmmTMB)
library(performance)
library(tidyverse)

# ---- Load data ----
data <- read_tsv("donnees_centauree.csv")

# ---- Sort date formatting ----

data$Date_de_germination <- as.Date(data$Date_de_germination,
                                    format = "%d/%m/%Y")

gathered_data <- data %>%
  rename("2005-12-01" = "Taille_Dec_05",
         "2006-02-01" = "Taille_Fev_06",
         "2006-03-01" = "Taille_Mars_06",
         "2006-06-01" = "Taille_Juin_06",
         "2006-09-01" = "Taille_Sept_06") %>%
  gather(key = "Date", 
         value = "Rosette_size", 
         "2005-12-01":"2006-09-01")

gathered_data$Date <- as.Date(gathered_data$Date, format = "%Y-%m-%d")

gathered_data$traitement <- factor(gathered_data$traitement,      # Reordering group factor levels
                                   levels = c("temoin", "peu", "mi", "dense"))

gathered_data <- gathered_data %>%
  mutate(light = PAR/reference)

# ---- Raw data plots

#xParRef_yRosette_facSpecDate_GrpPop
ggplot(gathered_data, aes(x = light, y = Rosette_size, colour = pop)) +
  geom_jitter(alpha = 0.3, size = 2) +
  facet_grid(espece~Date) +
  labs(x = "PAR/ref (1 = max light recieved & min competition)",
       y = "Rosette size") +
  theme(axis.text.x = element_text(angle = 90)) +
  theme(axis.text=element_text(size=12),
        strip.text = element_text(size = 12),
        axis.title=element_text(size=14,face="bold"))

ggsave(filename = "xParRef_yRosette_facSpecDate_GrpPop.png",
       dpi = "retina",
       height = 6,
       width = 12)

# x1minusParRef_yRosette_facSpecDate_GrpPop
ggplot(gathered_data, aes(x = (1-light), y = Rosette_size, colour = pop)) +
  geom_jitter(alpha = 0.3, size = 2) +
  facet_grid(espece~Date) +
  labs(x = "1-PAR/ref (0 = max light recieved & min competition)",
       y = "Rosette size") +
  theme(axis.text.x = element_text(angle = 90)) +
  theme(axis.text=element_text(size=12),
        strip.text = element_text(size = 12),
        axis.title=element_text(size=14,face="bold"))

ggsave(filename = "x1minusParRef_yRosette_facSpecDate_GrpPop.png",
       dpi = "retina",
       height = 6,
       width = 12)


# xDate_yRosette_facSpecTrait_GrpPop
ggplot(gathered_data, aes(x = Date, y = Rosette_size, colour = pop)) +
  geom_jitter(alpha = 0.3, size = 2) +
  facet_grid(espece~traitement) +
  scale_x_continuous(breaks = as.numeric(gathered_data$Date),
                     labels = format(gathered_data$Date, format = "%b")) +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(y = "Rosette size") +
  theme(axis.text=element_text(size=12),
        strip.text = element_text(size = 12),
        axis.title=element_text(size=14,face="bold"))

ggsave(filename = "xDate_yRosette_facSpecTrait_GrpPop.png",
       dpi = "retina",
       height = 6,
       width = 12)

# xTraitement_yRosette_facSpecDate_GrpPop
ggplot(gathered_data, aes(x = traitement, y = Rosette_size, colour = pop)) +
  geom_jitter(alpha = 0.3, size = 2) +
  facet_grid(espece~Date) +
  #scale_x_continuous(breaks = as.numeric(gathered_data$Date),
  #                   labels = format(gathered_data$Date, format = "%b")) +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(y = "Rosette size") +
  theme(axis.text=element_text(size=12),
        strip.text = element_text(size = 12),
        axis.title=element_text(size=14,face="bold"))

ggsave(filename = "xTraitement_yRosette_facSpecDate_GrpPop.png",
       dpi = "retina",
       height = 6,
       width = 12)

# ---- Time after germinations ----
# Maybe it is more informative to look at time after germination, not measurement date

d1 <- gathered_data %>%
  filter(!is.na(Date_de_germination)) %>%
  mutate(time_after_germ = difftime(Date, Date_de_germination, units = "days")) %>%
  mutate(time_after_germ = replace(time_after_germ, which(time_after_germ<0), NA))

d2 <- gathered_data %>%
  filter(is.na(Date_de_germination)) %>%
  mutate(time_after_germ = NA)

data_after_germ <- full_join(d1, d2)

ggplot(data = data_after_germ, aes(x = time_after_germ, y = Rosette_size, group = Plante, colour = pop)) +
  geom_line(alpha = 0.3) +
  geom_point(alpha = 0.1, size = 2) +
  facet_grid(espece~traitement) +
  labs(y = "Rosette size",
       x = "Days after germination") +
  theme(axis.text=element_text(size=12),
        strip.text = element_text(size = 12),
        axis.title=element_text(size=14,face="bold"))

ggsave(filename = "xtimeaftergerm_yRosette_facSpecTrait_GrpPop.png",
       dpi = "retina",
       height = 6,
       width = 12)

# ---- Modelling with 0 set to NA ----

nozero <- data_after_germ %>%
  mutate(Rosette_size = replace(Rosette_size, which(Rosette_size==0), NA))

lmer1 <- lmer(Rosette_size~light*time_after_germ*espece+Cotyledons+(1|Plante)+(1|pop),
              data = nozero)

summary(lmer1)
plot_model(lmer1)
plot_model(lmer1, type = "diag")

# pop is telling us nothing, so will remove

lmer2 <- lmer(Rosette_size~light*time_after_germ*espece+Cotyledons+(1|Plante),
              data = nozero)

summary(lmer2)
plot_model(lmer2)
plot_model(lmer2, type = "diag")

# will try log transform to rosette_size to try improve fit to assumption of norm dist res

lmer3 <- lmer(log(Rosette_size)~light*time_after_germ*espece+Cotyledons+(1|Plante),
              data = nozero)

summary(lmer3)
plot_model(lmer3)
plot_model(lmer3, type = "diag")

# not very sure, but will keep log atm
# will check if interactions are kept

lmer4 <- update(lmer3, .~. -light:time_after_germ:espece)
summary(lmer4)
plot_model(lmer4)
plot_model(lmer4, type = "diag")

compare_performance(lmer4, lmer3) # lmer4 is much better

# interaction between time and species isn't having much effect, will remove

lmer5 <- update(lmer4, .~. -time_after_germ:espece)
summary(lmer5)
plot_model(lmer5)
plot_model(lmer5, type = "diag")

# interaction between light and time is tiny

lmer6 <- update(lmer5, .~. -time_after_germ:light)
summary(lmer6)
plot_model(lmer6)
plot_model(lmer6, type = "diag")

compare_performance(lmer6, lmer5, lmer4, lmer3) # but aparently should be kept in, will come back to this

# try removeing time after germ

lmer7 <- update(lmer6, .~. -time_after_germ)
summary(lmer7)
plot_model(lmer7)
plot_model(lmer7, type = "diag")

compare_performance(lmer7, lmer6, lmer5, lmer4, lmer3) # definitly worse

# will try remove cot

lmer8 <- update(lmer6, .~. -Cotyledons)
summary(lmer8)
plot_model(lmer8)
plot_model(lmer8, type = "diag")

compare_performance(lmer8, lmer7, lmer6, lmer5, lmer4, lmer3) # also worse, back to model 5

summary(lmer5)
plot_model(lmer5, type = "int")
plot_model(lmer5, type = "pred")

ggplot()



# ---- GLMM (section not neat) ----
glm1 <- glmer(Rosette_size~light+(1|Plante),
              data = data_after_germ,
      family = "poisson")

withzeros <- data_after_germ %>%
  mutate(Rosette_size = replace(Rosette_size, which(Rosette_size==NA), 0))

glmm0 <- glmmTMB(Rosette_size~light+espece+time_after_germ+(1|Plante),
        data = withzeros,
        family = ziGamma(link = "log"),
        ziformula = ~1)

glmm1 <- glmmTMB(Rosette_size~light*espece+time_after_germ+(1|Plante),
        data = withzeros,
        family = ziGamma(link = "log"),
        ziformula = ~1)

glmm2 <- glmmTMB(Rosette_size~light*espece*time_after_germ+(1|Plante),
        data = withzeros,
        family = ziGamma(link = "log"),
        ziformula = ~1)

AIC(glmm0,glmm1,glmm2)

plot_model(glmm2, type = "pred")

glmm00 <- glmmTMB(Rosette_size~traitement+espece+time_after_germ+(1|Plante),
                 data = withzeros,
                 family = ziGamma(link = "log"),
                 ziformula = ~1)

glmm11 <- glmmTMB(Rosette_size~traitement*espece+time_after_germ+(1|Plante),
                 data = withzeros,
                 family = ziGamma(link = "log"),
                 ziformula = ~1)

glmm22 <- glmmTMB(Rosette_size~traitement*espece*time_after_germ+(1|Plante),
                 data = withzeros,
                 family = ziGamma(link = "log"),
                 ziformula = ~1)



AIC(glmm00,glmm11,glmm22)

check_model(glmm22)

plot_model(glmm22, type = "int")

# ---- TESTING ----

nozero <- data_after_germ %>%
  mutate(Rosette_size = replace(Rosette_size, which(Rosette_size==0), NA))

lmer1 <- lmmm(log(Rosette_size)~light+time_after_germ+espece+pop+Cotyledons+(1|Plante),
              data = nozero)

summary(lmer1)

plot_model(lmer1)
plot_model(lmer1, type = "diag")


