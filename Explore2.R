# ---- Clear Env ----
rm(list=ls())

# ---- Packages ----
library(reshape2)
library(lme4)
library(sjPlot)
library(glmmTMB)
library(performance)
library(equatiomatic)
library(ggeffects)
library(cowplot)
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

rm(d1,d2)

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

# ---- Modelling with max size ---- #

maxdata <- data_after_germ %>%
  group_by(Plante) %>%
  filter(Rosette_size == max(Rosette_size)) %>%
  filter(time_after_germ == max(time_after_germ)) # gets rid of duplicates and keeps oldest value

maxdata<-maxdata %>%
  mutate(comp = 1-light)

maxdata <- maxdata %>%
  filter(Rosette_size > 0)

maxdata <- maxdata %>%
  rename("Competition" = "comp",
         "Light" = "light",
         "Species" = "espece",
         "Population" = "pop")

ggplot(data = maxdata, aes(x = Competition, y = Rosette_size, colour = Population)) +
  geom_point(alpha = 0.3, size = 4) +
  #geom_smooth(method = "lm") +
  facet_grid(~Species) +
  labs(y = "Max rosette size achieved",
       x = "Competition (0 - low, 1 - high)") +
  theme(axis.text=element_text(size=12),
        strip.text = element_text(size = 12),
        axis.title=element_text(size=14,face="bold"))

ggsave(filename = "maxgrowth.png",
       dpi = "retina",
       height = 6,
       width = 12)

maxdata %>%
  group_by(espece, pop) %>%
  count() # very unbalanced

pairs(maxdata)

hist(maxdata$Rosette_size)

ggplot()

# messing around with aov
aov1<-aov(Rosette_size~espece:pop, data = maxdata)
summary(aov1)
TukeyHSD(aov1)

lmer1 <- lmer(Rosette_size~Light*Species+Cotyledons+(1+Light|Species:Population), data = maxdata)
summary(lmer1)



deviance(lmer1)

plot_model(lmer1, type = "diag")
plot_model(lmer1, type = "re")
plot_model(lmer1, type = "int")

extract_eq(lmer1.5)


maxdata$Species <- recode(maxdata$Species,mac = "a", cor = "b")
lmer1.5 <- lmer(Rosette_size~Light*Species+Cotyledons+(1|Species:Population), data = maxdata)

summary(lmer1.5)
plot_model(lmer1.5, type = "resid")
plot_model(lmer1.5, type = "int")

AIC(lmer1, lmer1.5, lmer1.6, lmer2, lmer3, lmer4, lmer5, lmer6, lmer7)

lmer1.6 <- lm(Rosette_size~Light*Species+Cotyledons, data = maxdata)
summary(lmer1.6)
plot_model(lmer1.6, type = "diag")
plot_model(lmer1.6, type = "int")

extract_eq(lmer1.5, use_coefs = FALSE)


anova(lmer1.5, lmer1, dist = "Chi")
anova(lmer1.5, lmer1.6, dist = "Chi")
drop1(lmer1, test = "Chisq")

# random effect can be removed

lmer2 <- update(lmer1.6, .~. -Cotyledons)
summary(lmer2)
plot_model(lmer2)
plot_model(lmer2, type = "diag")

anova(lmer2, lmer1.6) # cotyledons should be kept in

lmer3 <- update(lmer1.6, .~. -Light:Species)
summary(lmer3)
plot_model(lmer3)

anova(lmer3, lmer1.6) # interaction between light and species can be removed

lmer4 <- update(lmer3, .~. -Species)
summary(lmer4)
plot_model(lmer4)
plot_model(lmer4, type = "pred")
anova(lmer4, lmer3) # species is supported to be removed

extract_eq(lmer4,show_distribution = TRUE, use_coefs = TRUE)

lmer5 <- update(lmer4, .~. -Cotyledons)
summary(lmer5)
plot_model(lmer5)
anova(lmer5, lmer4) # cotyledons should be kept still

lmer6 <- update(lmer4, .~. -Light)
summary(lmer6)
plot_model(lmer6)
anova(lmer6, lmer4) # light should be kept in

lmer7 <- update(lmer5, .~. -Light)
summary(lmer7)
anova(lmer7, lmer4)



lmtest <- lm(Rosette_size~Light*Cotyledons, data = maxdata)
plot_model(lmtest, type = "int")
anova(lmtest)
summary(lmtest)

lmtest2 <- lm(Rosette_size~Competition, data = maxdata)
plot_model(lmtest2, type = "diag")
anova(lmtest2)
summary(lmtest2)


plot_model(lmer4, type = "pred",terms = c("Light"),show.data = TRUE, dot.size = 2) +
  xlim(0,1.08)
plot_model(lmer4, type = "pred",terms = c("Cotyledons"),show.data = TRUE, dot.size = 2)

lmer4pred <- ggpredict(lmer4, terms = "Light")
lmerpred1 <- ggpredict(lmer4, terms = "Cotyledons")

g1 <- ggplot(data = maxdata, aes(x = Light, y = Rosette_size)) +
  geom_point(shape=21, fill = "black", size = 3, alpha = 0.3) +
  geom_line(data = lmer4pred, aes(y = predicted, x = x), size = 1.5) +
  geom_ribbon(data = lmer4pred, aes(y = predicted, x = x, ymin = conf.low, ymax = conf.high), alpha = 0.1) +
  labs(y = "Maximum rosette size (mm)",
       x = "Light (received PAR / reference PAR)") +
  scale_x_continuous(expand = c(0, 0), limits = c(0,1.09)) + scale_y_continuous(expand = c(0, 0), limits = c(0,150)) +
  annotate(geom="text", x=0.80, y=5, label="*Adjusted for Cotyledon = 9.40",
             color="black") +
  theme_cowplot()

g2 <- ggplot(data = maxdata, aes(x = Cotyledons, y = Rosette_size)) +
  geom_point(shape=21, fill = "black", size = 3, alpha = 0.3) +
  geom_line(data = lmerpred1, aes(y = predicted, x = x), size = 1.5) +
  geom_ribbon(data = lmerpred1, aes(y = predicted, x = x, ymin = conf.low, ymax = conf.high), alpha = 0.1) +
  labs(y = "Maximum rosette size (mm)",
       x = "Cotyledon size (mm)") +
  scale_x_continuous(expand = c(0, 0), limits = c(0,26)) + scale_y_continuous(expand = c(0, 0), limits = c(0,150)) +
  annotate(geom="text", x=20, y=5, label="*Adjusted for Light = 0.73",
           color="black") +
  theme_cowplot() +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank())

g1g2 <- grid.arrange(g1,g2, ncol = 2, nrow = 1)

ggsave(filename = "test.png",
       plot = g1g2,
       dpi = "retina",
       height = 5,
       width = 10)

lmer15pred <- ggpredict(lmer1.5, terms = c("Light", "Species"),type = "fixed")
lmer15pred <- lmer15pred %>%
  rename("Species" = "group")

lmer15pred$Species <- recode_factor(lmer15pred$Species, mac = "C. maculosa", cor = "C. corymbosa")

lmer151pred <- ggpredict(lmer1.5, terms = c("Cotyledons", "Species"),type = "fixed")
lmer151pred <- lmer151pred %>%
  rename("Species" = "group")

lmer151pred$Species <- recode_factor(lmer151pred$Species, mac = "C. maculosa", cor = "C. corymbosa")

maxdata1 <- maxdata
maxdata1$Species <- recode_factor(maxdata1$Species, mac = "C. maculosa", cor = "C. corymbosa")


#temp1 <- lmer15pred %>%
 # filter(conf.low < 0) %>%
  #mutate(conf.low = 0)

#temp2 <- lmer15pred %>%
 # filter(conf.low > 0)

#lmer15pred <- full_join(temp1, temp2)

#rm(temp1, temp2)

g1 <- ggplot(data = maxdata1, aes(x = Light, y = Rosette_size)) +
  geom_point(data = maxdata1 %>%
               filter(Species == "a"), aes(x = Light, y = Rosette_size), size = 3, alpha = 0.3, colour = "#F8766D") +
  geom_point(data = maxdata1 %>%
               filter(Species == "b"), aes(x = Light, y = Rosette_size), size = 3, alpha = 0.3, colour = "#00BFC4") +
  geom_line(data = lmer15pred, aes(y = predicted, x = x, colour = Species), size = 1.5) +
  geom_ribbon(data = lmer15pred, aes(y = predicted, x = x, ymin = conf.low, ymax = conf.high, colour = Species, fill = Species), alpha = 0.1) +
  labs(y = "Maximum rosette size (mm)",
       x = "Light (received PAR / reference PAR)",
       title = "A") +
  scale_x_continuous(expand = c(0, 0), limits = c(0,1.09)) + scale_y_continuous(expand = c(0, 0), limits = c(0,150)) +
  annotate(geom="text", x=0.80, y=5, label="*Adjusted for Cotyledon = 9.40",
           color="black") +
  theme_cowplot()+
  theme(legend.position = "none")

g2 <- ggplot(data = maxdata1, aes(x = Cotyledons, y = Rosette_size)) +
  geom_point(data = maxdata1 %>%
               filter(Species == "a"), aes(x = Cotyledons, y = Rosette_size), size = 3, alpha = 0.3, colour = "#F8766D") +
  geom_point(data = maxdata1 %>%
               filter(Species == "b"), aes(x = Cotyledons, y = Rosette_size), size = 3, alpha = 0.3, colour = "#00BFC4") +
  geom_line(data = lmer151pred, aes(y = predicted, x = x, colour = Species), size = 1.5) +
  geom_ribbon(data = lmer151pred, aes(y = predicted, x = x, ymin = conf.low, ymax = conf.high, colour = Species, fill = Species), alpha = 0.1) +
  labs(y = "Maximum rosette size (mm)",
       x = "Cotyledon size (mm)",
       title = "B") +
  scale_x_continuous(expand = c(0, 0), limits = c(0,26)) + scale_y_continuous(expand = c(0, 0), limits = c(0,150)) +
  annotate(geom="text", x=21, y=5, label="*Adjusted for Light = 0.73",
           color="black") +
  theme_cowplot() +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        legend.position = c(0.7,0.9),
        legend.title = element_blank())

g1g2 <- grid.arrange(g1,g2, ncol = 2, nrow = 1)

ggsave(filename = "Figure2.png",
       plot = g1g2,
       dpi = "retina",
       height = 5,
       width = 10)


# ---- Germ Data ----

didnotgerm <- data %>%
  filter(is.na(Date_de_germination)) %>%
  mutate(Date_de_germination = 0)
didgerm <- data %>%
  filter(!is.na(Date_de_germination)) %>%
  mutate(Date_de_germination = 1)

binarydata <- full_join(didnotgerm, didgerm)

binarydata<-binarydata %>%
  mutate(light = PAR/reference)

binarydata<-binarydata %>%
  mutate(comp = 1-light)

binarydata <- binarydata %>%
  rename("Germination" = "Date_de_germination",
         "Competition" = "comp",
         "Light" = "light",
         "Species" = "espece",
         "Population" = "pop")

binarydata %>%
  group_by(Species) %>%
  count()


ggplot(data = binarydata, aes(y = Germination, x = Light, colour = Species)) +
  geom_point(alpha = 0.1, size = 4)

binarydata$Species <- recode(binarydata$Species,mac = "a", cor = "b")


glmer1 <- glmer(Germination~Species*Light+(1|Species:Population), data = binarydata, family = "binomial")
summary(glmer1)
plot_model(glmer1, type = "re")

model_performance(glmer1)


pred <- data.frame(Species = levels(binarydata$Species))
predict(glmer1, newdat = pred, type = "response")

glmer1sim <- simulate(glmer1, 100000) %>%
  colSums()

hist(glmer1sim, breaks = 100)
abline(v = 488, col = "red")

glmer1sim %>%
  colSums()

binarydata$Germination %>%
  sum()

anova(glmer1)

glmer1.5 <- glm(Germination~Species*Light, data = binarydata, family = "binomial")
glmer1.5 <- lmer(Germination~Species*Light+(1|Species:Population), data = binarydata)

glmer2 <- update(glmer1, .~. -Species:Light)
summary(glmer2)

glmer3 <- update(glmer2, .~. -Light)
summary(glmer3)

anova(glmer1, glmer2)

model_performance(glmer1)
model_performance(glmer1.5)
model_performance(glmer2)
model_performance(glmer3)

extract_eq(glmer1.5, show_distribution = TRUE)

plot_model(glmer1, type = "pred", terms = c("Light", "Species"), show.data = TRUE, dot.size = 4, ci.style = "bar") +
  xlim(0,1.03)



glmerpred <- ggpredict(glmer1, terms = c("Light [all]", "Species"),type = "fixed")

glmerpred <- glmerpred %>%
  rename("Species" = "group")

glmerpred$Species <- recode_factor(glmerpred$Species, a = "C. maculosa", b = "C. corymbosa")

ggplot(data = binarydata, aes(x = Light, y = Germination)) +
  geom_point(data = binarydata %>%
               filter(Species == "a"), aes(x = Light, y = Germination), size = 3, alpha = 0.1, colour = "#F8766D") +
  geom_point(data = binarydata %>%
               filter(Species == "b"), aes(x = Light, y = Germination), size = 3, alpha = 0.1, colour = "#00BFC4") +
  geom_line(data = glmerpred, aes(y = predicted, x = x, colour = Species), size = 1.5) +
  geom_ribbon(data = glmerpred, aes(y = predicted, x = x, ymin = conf.low, ymax = conf.high, colour = Species, fill = Species), alpha = 0.1) +
  labs(y = "P(Germination)",
       x = "Light (received PAR / reference PAR)") +
 # scale_x_continuous(expand = c(0, 0), limits = c(0,1.09)) + scale_y_continuous(expand = c(0, 0), limits = c(0,150)) +
  #annotate(geom="text", x=0.80, y=5, label="*Adjusted for Cotyledon = 9.40",
   #        color="black") +
  scale_x_continuous(breaks = c(0,0.25,0.5,0.75,1)) +
  theme_cowplot() +
  theme(legend.position = c(0.65,0.25),
        legend.title = element_blank())

ggsave(filename = "Figure1.png",
       dpi = "retina",
       height = 5,
       width = 5)

# ---- Modelling with 0 set to NA ----

nozero <- data_after_germ %>%
  mutate(Rosette_size = replace(Rosette_size, which(Rosette_size==0), NA))

ggplot(data = nozero, aes(x = time_after_germ, y = Rosette_size, group = Plante, colour = pop)) +
  geom_line(alpha = 0.3) +
  geom_point(alpha = 0.1, size = 2) +
  facet_grid(espece~traitement) +
  labs(y = "Rosette size",
       x = "Days after germination") +
  theme(axis.text=element_text(size=12),
        strip.text = element_text(size = 12),
        axis.title=element_text(size=14,face="bold"))

nozero %>%
  group_by(espece) %>%
  count()


lmer1 <- lmer(Rosette_size~light*espece*time_after_germ+Cotyledons+(1|Plante)+(1|pop),
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

lmer3 <- lmer(log(Rosette_size)~light*espece*time_after_germ+Cotyledons+(1|Plante),
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

ggsave(filename = "model5homo.png",
       width = 7,
       height = 5,
       dpi = "retina")
plot_model(lmer5, type = "diag")

lmer5nolog <- lmer(Rosette_size ~ light + espece + time_after_germ + Cotyledons +  (1 | Plante) + light:espece + light:time_after_germ,
     data = nozero)

lmer5diag <- plot_model(lmer5, type = "diag")

ggsave(filename = "logmodel5diag1.png",
       plot = lmer5diag[[1]],
       width = 7,
       height = 5,
       dpi = "retina")
ggsave(filename = "logmodel5diag2.png",
       plot = lmer5diag[[2]]$Plante,
       width = 7,
       height = 5,
       dpi = "retina")
ggsave(filename = "logmodel5diag3.png",
       plot = lmer5diag[[3]],
       width = 7,
       height = 5,
       dpi = "retina")
ggsave(filename = "logmodel5diag4.png",
       plot = lmer5diag[[4]],
       width = 7,
       height = 5,
       dpi = "retina")

lmer5nologdiag <- plot_model(lmer5nolog, type = "diag")

ggsave(filename = "model5diag1.png",
       plot = lmer5nologdiag[[1]],
       width = 7,
       height = 5,
       dpi = "retina")
ggsave(filename = "model5diag2.png",
       plot = lmer5nologdiag[[2]]$Plante,
       width = 7,
       height = 5,
       dpi = "retina")
ggsave(filename = "model5diag3.png",
       plot = lmer5nologdiag[[3]],
       width = 7,
       height = 5,
       dpi = "retina")
ggsave(filename = "model5diag4.png",
       plot = lmer5nologdiag[[4]],
       width = 7,
       height = 5,
       dpi = "retina")

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
plot_model(lmer5)
predint <- plot_model(lmer5, type = "int")

predint[[1]]
ggsave(filename = "predint.png",
       width = 5,
       height = 3,
       dpi = "retina")


plot_model(lmer5, type = "pred")



testMOD <- lmer(log(Rosette_size)~time_after_germ*light*espece+Cotyledons+(1|Plante),
     data = nozero)

plot_model(testMOD)

str(nozero)

nozero$time_after_germ <- as.numeric(nozero$time_after_germ) # important

hist(nozero$time_after_germ)

lmerTEST <- lmer(log(Rosette_size)~espece*light*time_after_germ+Cotyledons+(1|Plante), data = nozero)
plot_model(lmerTEST)
plot_model(lmerTEST, type = "int")
plot_model(lmerTEST, type = "pred")
equatiomatic::extract_eq(lmTEST, wrap = TRUE)

lmTEST <- lmer(Rosette_size~light*espece*time_after_germ+Cotyledons+(1|Plante), data = nozero)

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

summary(glmm2)

AIC(glmm0,glmm1,glmm2)

plot_model(glmm2, type = "pred")
plot_model(glmm2, type = "diag")

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


