# Data Exploration Scrip
# IR Moodie
# 16/02/21
# Last modified: 16/02/21

# ---- Clear Environment ----
rm(list=ls())

# ---- Load Packages ----
library(tidyverse)
library(reshape2)

# ---- Load Data ----
data <- read_tsv("donnees_centauree.csv")

# ---- Clean Data ----
# Date format
data$Date_de_germination <- as.Date(data$Date_de_germination, format = "%d/%m/%Y")

head(data)

# If plant didn't germinate, it should not have a rosette size
data_nogerm <- data %>%
  filter(is.na(Date_de_germination)) %>%
  mutate(Taille_Dec_05 = NA,
         Taille_Fev_06 = NA,
         Taille_Mars_06 = NA,
         Taille_Juin_06 = NA,
         Taille_Sept_06 = NA)

data_germ <- data %>%
  filter(!is.na(Date_de_germination))

fulldata <- full_join(data_nogerm, data_germ)

fulldata %>%
  select(Taille_Dec_05) %>%
  rename("Dec 2005" = "Taille_Dec_05")

gathered_data <- fulldata %>%
  rename("2005-12-01" = "Taille_Dec_05",
         "2006-02-01" = "Taille_Fev_06",
         "2006-03-01" = "Taille_Mars_06",
         "2006-06-01" = "Taille_Juin_06",
         "2006-09-01" = "Taille_Sept_06") %>%
  gather(key = "Date", value = "Rosette_size", "2005-12-01":"2006-09-01")

gathered_data$Date <- as.Date(gathered_data$Date, format = "%Y-%m-%d")

ggplot(light, aes(x = light, y = Rosette_size)) +
  geom_jitter(alpha = 0.3) +
  facet_grid(~Date)

ggsave(filename = "testplot2.png")

# ---- Explore ----

light <- gathered_data %>%
  mutate(light = PAR/reference)

lm1 <- lm(light~traitement, data = light)

anova(lm1)






# Dates of germination
ggplot(data = data, aes(x = Date_de_germination)) +
  geom_histogram(binwidth = 14) +
  scale_x_date(date_breaks = "months", date_labels = "%b-%y") +
  theme_minimal()

ggplot(data = data, aes)







plot(traitement~PAR, data = data)


ggplot(data %>%
         filter(traitement != "temoin"),
       aes(y = PAR, x = traitement)) +
  geom_jitter(width = 0.2)

ggplot(data %>%
         filter(PAR != "100"),
       aes(x = PAR, y = Cotyledons)) +
  geom_point()

ggplot(data,
       aes(x = pop, y = Cotyledons)) +
  geom_jitter(width = 0.2) +
  

ggplot(data %>%
         filter(traitement != "temoin"),
       aes(y = PAR, x = traitement)) +
  geom_jitter(width = 0.2)
  


ggplot(data %>%
         filter(Date_de_germination != NA),
       aes(x = Date_de_germination)) +
  geom_histogram(stat = "count")



data$Date_de_germination <- strptime(data$Date_de_germination, "%d/%m/%Y")       

str(data)

data$Date_de_germination <- as.Date(data$Date_de_germination, format = "%d/%m/%Y")
