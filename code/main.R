# ========== Librairies ========== #
library(survival)
library(lattice)
library(Formula)
library(ggplot2)
library(Hmisc)
library(foreign)
library(nlme) 
library(splines) 
library(MASS)
library(JM) 
library(tidyverse) 
library(tidyr)
library(magrittr) 
library(ggkm) 
library(parallel)
library(mvtnorm)
library(lcmm)
library(nnet)
library(epiDisplay)
library(epiR)
library(lme4)
library(lmerTest)
library(moments)
library(ggpubr)
library(bestNormalize)

# ========== Importation et construction de la base  ========== #

# ____ Définir le répertoire ____ #
setwd("/home/tuvu/Downloads/projet_STA302/projet")

# ____ Importation ____ #

df <- read.csv("projet.csv", sep=";", dec = ".", na.strings = NA)

# ____ Types de variables ____ #

df$APOE4 <- as.factor(df$APOE4)
df$DIPNIV0 <- as.factor(df$DIPNIV0)
df$DEM0_8 <- as.factor(df$DEM0_8)
df$DC8 <- as.factor(df$DC8)
df$STATUT5 <- as.factor(df$STATUT5)
df$STATUT6 <- as.factor(df$STATUT6)
df$STATUT7 <- as.factor(df$STATUT7)
df$STATUT8 <- as.factor(df$STATUT8)
df$AGE5 <- as.numeric(df$AGE5)
df$AGE6 <- as.numeric(df$AGE6)
df$AGE7 <- as.numeric(df$AGE7)
df$AGE8 <- as.numeric(df$AGE8)
df$AGEDEM8 <- as.numeric(df$AGEDEM8)
df$AGEFIN8 <- as.numeric(df$AGEFIN8)
df$MMSTT5 <- as.numeric(df$MMSTT5)
df$MMSTT6 <- as.numeric(df$MMSTT6)
df$MMSTT7 <- as.numeric(df$MMSTT7)
df$MMSTT8 <- as.numeric(df$MMSTT8)
df$ISA5_60 <- as.numeric(df$ISA5_60)
df$ISA6_60 <- as.numeric(df$ISA6_60)
df$ISA7_60 <- as.numeric(df$ISA7_60)
df$ISA8_60 <- as.numeric(df$ISA8_60)
df$BENTON5 <- as.numeric(df$BENTON5)
df$BENTON6 <- as.numeric(df$BENTON6)
df$BENTON7 <- as.numeric(df$BENTON7)
df$BENTON8 <- as.numeric(df$BENTON8)
df$LG_AX_OD <- as.numeric(df$LG_AX_OD)
df$LG_AX_OG <- as.numeric(df$LG_AX_OG)

# ____ Création de nouvelles variables ____ #

df$RNFLG <- ifelse(!is.na(df$RNFLGD_1), df$RNFLGD_1, df$RNFLGG_1) # RNFLG droit si pas NA, RNGLG gauche sinon
df$AX <- ifelse(!is.na(df$RNFLGD_1), df$LG_AX_OD, df$LG_AX_OG) # AXIALE droit si RNFLG droit pas NA, AXIALE gauche sinon

# ____ Création du jeu de données de travail "base" ____ #

df$AGEINT <- df$AGE5

df1 <- gather(df, MMSTTX, MMSE, MMSTT5, MMSTT6, MMSTT7, MMSTT8)
df1 <- df1[order(df1$ID, df1$MMSTTX), ]

df2 <- gather(df, ISAX_60, ISAAC, ISA5_60, ISA6_60, ISA7_60, ISA8_60)
df2 <- df2[order(df2$ID, df2$ISAX_60), ]

df3 <- gather(df, BENTONX, BENTON, BENTON5, BENTON6, BENTON7, BENTON8)
df3 <- df3[order(df3$ID, df3$BENTONX), ]

df$STATUT5 <- as.character(df$STATUT5)
df$STATUT6 <- as.character(df$STATUT6)
df$STATUT7 <- as.character(df$STATUT7)
df$STATUT8 <- as.character(df$STATUT8)

df4 <- gather(df, STATUTX, STATUT, STATUT5, STATUT6, STATUT7, STATUT8)
df4 <- df4[order(df4$ID, df4$STATUTX), ]
df4$STATUTX <- as.factor(df4$STATUTX)
df4$STATUT <- as.factor(df4$STATUT)

df5 <- gather(df, AGEX, AGE, AGE5, AGE6, AGE7, AGE8)
df5 <- df5[order(df5$ID, df5$AGEX), ]

df6 <- cbind(df1, df2[c("ISAX_60", "ISAAC")])
df7 <- cbind(df6, df3[c("BENTONX", "BENTON")])
df8 <- cbind(df7, df4[c("STATUTX", "STATUT")])
df9 <- cbind(df8, df5[c("AGEX", "AGE")])

df9$DIPNIV <- df9$DIPNIV0

df9$DIPNIV <- as.character(df9$DIPNIV)

df9$DIPNIV[df9$DIPNIV == "1" | df9$DIPNIV == "2"] <- "0"
df9$DIPNIV[df9$DIPNIV == "3" | df9$DIPNIV == "4"] <- "1"
df9$DIPNIV[df9$DIPNIV == "5"] <- "2"
df9$DIPNIV <- as.factor(df9$DIPNIV)

df9$SEXE[df9$SEXE == "1"] <- "0"
df9$SEXE[df9$SEXE == "2"] <- "1"
df9$SEXE <- as.factor(df9$SEXE)

data <- df9[,c("ID", "SEXE", "DIPNIV", "AGEX", "AGE", "AGEINT", "AGEDEM8", "AGEFIN8", "DEM0_8", "DC8", "APOE4", "RNFLG", "RNFLGD_1", "RNFLGG_1", "AX", "MMSTTX", "MMSE", "ISAX_60", "ISAAC", "BENTONX", "BENTON", "STATUTX", "STATUT")]

# Supprimer les objets inutiles de l'environnement R
rm(df1, df2, df3, df4, df5, df6, df7, df8, df9)

# Création de la variable de délai
data$OBSTIME <- (data$AGE - data$AGEINT)/10

# Création de la variable de AGE_MEAN
#data$AGE_MEAN <- (data$AGEINT-mean(data$AGE, na.rm=TRUE))/10
data$AGE75 <- (data$AGEINT-75)/10
data$AGE65 <- (data$AGEINT-65)/10

# ========== Vérifications ========== #

# Données manquantes dans "data"
sapply(data, function(x) sum(is.na(x)))
sum(is.na(data))

# ========== Analyse descriptive ========== #

summary(data)

# Nombre d'observations
df$ID %>% unique() %>% length()

str(data)

# Vérification de la distribution de ISAAC
hist(data$ISAAC) # ISAAC suit bien la loi normale

# ========== Graphiques ========== #

p1 = ggplot(data, aes(x = OBSTIME, y = ISAAC, group = ID)) +
  geom_line(aes(color=ID))+
  geom_point()+
  scale_x_continuous("OBSTIME", breaks=seq(0, 14, 1)) +
  ggtitle("ISAAC scores by Indivisuals")

p2 = ggplot(data, aes(x = OBSTIME, y = ISAAC)) +
  geom_smooth(method = 'loess')+
  geom_point()+
  scale_x_continuous("OBSTIME", breaks=seq(0, 14, 1)) +
  ggtitle("Overall ISAAC scores Trend")

gridExtra::grid.arrange(p1,p2,ncol=2)

ggplot(data, aes(x=OBSTIME, y=ISAAC, group=ID)) +
  geom_line(aes(color=ID))+
  geom_point(aes(color=ID))+
  facet_grid(. ~ SEXE)

ggplot(data, aes(x=OBSTIME, y=ISAAC)) +
  geom_smooth(method = 'loess')+
  geom_point()+
  facet_grid(. ~ SEXE)

ggplot(data, aes(x=OBSTIME, y=ISAAC, group=ID)) +
  geom_line(aes(color=ID))+
  geom_point(aes(color=ID))+
  facet_grid(. ~ APOE4)

ggplot(data, aes(x=OBSTIME, y=ISAAC)) +
  geom_smooth(method = 'loess')+
  geom_point()+
  facet_grid(. ~ APOE4)

ggplot(data, aes(x=OBSTIME, y=ISAAC, group=ID)) +
  geom_line(aes(color=ID))+
  geom_point(aes(color=ID))+
  facet_grid(. ~ DIPNIV)

ggplot(data, aes(x=OBSTIME, y=ISAAC)) +
  geom_smooth(method = 'loess')+
  geom_point()+
  facet_grid(. ~ DIPNIV)

### =========== TEST LA NORMALITE DES SCORES COGNITIFS ============= ###
hist(data$ISAAC) 
hist(data$MMSE) 
hist(data$BENTON)

# ISAAC suit bien la loi normale.
# Par contre, MMSE et BENTON ne suit pas la loi normale. Alors, il est important de faire la transformation de données pour MMSE et BENTON.

# ======================= Partie 1: Modèles linéaires mixtes ============================= #

######## ============ ETUDE DE LA RELATION ENTRE RNFL ET SCORE ISAAC #######################

# === Comparaison d'AIC entre modèle à intercept aléatoire et modèles à intercept et pente aléatoire === #

# Modèle à intercept aléatoire
mi1ISAAC <- lme(fixed = ISAAC ~ OBSTIME + RNFLG + AGE75 + SEXE + DIPNIV + APOE4, data = data, random = ~ 1 | ID, method = "ML", na.action = na.omit)  
summary(mi1ISAAC) # AIC: 12452.91

# Modèle à intercept et pente aléatoires
mip2ISAAC <- lme(fixed = ISAAC ~ OBSTIME + RNFLG + AGE75 + SEXE + DIPNIV + APOE4, data = data, random = ~ OBSTIME | ID, method = "ML", na.action = na.omit)  
summary(mip2ISAAC) # AIC: 12344.76

### --- On garde finalement le modèle à intercept et pente aléatoires --- #

### ============ TEST UNIVARIE DES INTERACTIONS ENTRE CHAQUE VARIABLE ET LE DELAI ============== ###

### === RNFLG * OBSTIME === ###
mu1ISAAC <- lme(fixed = ISAAC ~ OBSTIME + RNFLG + AGE75 + SEXE + DIPNIV + APOE4 + RNFLG*OBSTIME, data = data, random = ~ OBSTIME | ID, method = "ML", na.action = na.omit)  
summary(mu1ISAAC) 

### === AGE75 * OBSTIME === ###
mu2ISAAC <- lme(fixed = ISAAC ~ OBSTIME + RNFLG + AGE75 + SEXE + DIPNIV + APOE4 + AGE75*OBSTIME, data = data, random = ~ OBSTIME | ID, method = "ML", na.action = na.omit)  
summary(mu2ISAAC) 

### === SEXE * OBSTIME === ###
mu3ISAAC <- lme(fixed = ISAAC ~ OBSTIME + RNFLG + AGE75 + SEXE + DIPNIV + APOE4 + SEXE*OBSTIME, data = data, random = ~ OBSTIME | ID, method = "ML", na.action = na.omit)  
summary(mu3ISAAC) 

### === DIPNIV * OBSTIME === ###
mu4ISAAC <- lme(fixed = ISAAC ~ OBSTIME + RNFLG + AGE75 + SEXE + DIPNIV + APOE4 + DIPNIV*OBSTIME, data = data, random = ~ OBSTIME | ID, method = "ML", na.action = na.omit)  
summary(mu4ISAAC) 

### === APOE4 * OBSTIME === ###
mu5ISAAC <- lme(fixed = ISAAC ~ OBSTIME + RNFLG + AGE75 + SEXE + DIPNIV + APOE4 + APOE4*OBSTIME, data = data, random = ~ OBSTIME | ID, method = "ML", na.action = na.omit)  
summary(mu5ISAAC) 

### === MODELE MULTIVARIE === ###
mfinal_ISAAC <- lme(fixed = ISAAC ~ OBSTIME + RNFLG + AGE75 + SEXE + APOE4 + DIPNIV + OBSTIME*(RNFLG + AGE75 + APOE4), data = data, random = ~ OBSTIME | ID, method = "ML", na.action = na.omit)
summary(mfinal_ISAAC)

####### ====== VALIDATION DU MODELE ======== #######
q1ISAAC <- qqnorm(mfinal_ISAAC)
pl1ISAAC <- plot(mfinal_ISAAC)
gridExtra::grid.arrange(q1ISAAC,pl1ISAAC,ncol=2)

# -------------------- TRANSFORMATION DE SCORE BENTON ------------------------ #

## racine carrée pour biais modéré: sqrt(max(x+1) - x) pour les données asymétriques négatives
data$BENTONsqrtT <- sqrt(max(data$BENTON + 1,na.rm = TRUE) - data$BENTON)

######## ============ ETUDE DE LA RELATION ENTRE RNFL ET SCORE BENTON #######################

# === Comparaison d'AIC entre modèle à intercept aléatoire et modèles à intercept et pente aléatoire === #

# Modèle à intercept aléatoire
mi1BENTONsqrtT <- lme(fixed = BENTONsqrtT ~ OBSTIME + RNFLG + AGE75 + SEXE + DIPNIV + APOE4, data = data, random = ~ 1 | ID, method = "ML", na.action = na.omit)  
summary(mi1BENTONsqrtT) # AIC: 1779.72

# Modèle à intercept et pente aléatoires
mip2BENTONsqrtT <- lme(fixed = BENTONsqrtT ~ OBSTIME + RNFLG + AGE75 + SEXE + DIPNIV + APOE4, data = data, random = ~ OBSTIME | ID, method = "ML", na.action = na.omit)  
summary(mip2BENTONsqrtT) # AIC: 1780.538

### --- On garde finalement le modèle à intercept et pente aléatoires --- #

### ============ TEST UNIVARIE DES INTERACTIONS ENTRE CHAQUE VARIABLE ET LE DELAI ============== ###

### === RNFLG * OBSTIME === ###
mu1BENTONsqrtT <- lme(fixed = BENTONsqrtT ~ OBSTIME + RNFLG + AGE75 + SEXE + DIPNIV + APOE4 + RNFLG*OBSTIME, data = data, random = ~ OBSTIME | ID, method = "ML", na.action = na.omit)  
summary(mu1BENTONsqrtT) # 0.9274

### === AGE75 * OBSTIME === ###
mu2BENTONsqrtT <- lme(fixed = BENTONsqrtT ~ OBSTIME + RNFLG + AGE75 + SEXE + DIPNIV + APOE4 + AGE75*OBSTIME, data = data, random = ~ OBSTIME | ID, method = "ML", na.action = na.omit)  
summary(mu2BENTONsqrtT) # 0.0339

### === SEXE * OBSTIME === ###
mu3BENTONsqrtT <- lme(fixed = BENTONsqrtT ~ OBSTIME + RNFLG + AGE75 + SEXE + DIPNIV + APOE4 + SEXE*OBSTIME, data = data, random = ~ OBSTIME | ID, method = "ML", na.action = na.omit)  
summary(mu3BENTONsqrtT) # 0.0989

### === DIPNIV * OBSTIME === ###
mu4BENTONsqrtT <- lme(fixed = BENTONsqrtT ~ OBSTIME + RNFLG + AGE75 + SEXE + DIPNIV + APOE4 + DIPNIV*OBSTIME, data = data, random = ~ OBSTIME | ID, method = "ML", na.action = na.omit)  
summary(mu4BENTONsqrtT) # 0.4094, 0.9236

### === APOE4 * OBSTIME === ###
mu5BENTONsqrtT <- lme(fixed = BENTONsqrtT ~ OBSTIME + RNFLG + AGE75 + SEXE + DIPNIV + APOE4 + APOE4*OBSTIME, data = data, random = ~ OBSTIME | ID, method = "ML", na.action = na.omit)  
summary(mu5BENTONsqrtT) # 0.0290

### === MODELE MULTIVARIE === ###
mfinal_BENTONsqrtT <- lme(fixed = BENTONsqrtT ~ OBSTIME + RNFLG + AGE75 + SEXE + APOE4 + DIPNIV + OBSTIME*(RNFLG + AGE75 + APOE4), data = data, random = ~ OBSTIME | ID, method = "ML", na.action = na.omit)
summary(mfinal_BENTONsqrtT)

####### ====== VALIDATION DU MODELE ======== #######
q1BENTONsqrtT <- qqnorm(mfinal_BENTONsqrtT)
pl1BENTONsqrtT <- plot(mfinal_BENTONsqrtT)
gridExtra::grid.arrange(q1BENTONsqrtT,pl1BENTONsqrtT,ncol=2)


######## ============ ETUDE DE LA RELATION ENTRE RNFL ET SCORE MMSE #######################
## log pour une plus grande asymétrie: log10(max(x+1) - x) pour les données asymétriques négatives
data$logMMSE <- log10(max(data$MMSE + 1,na.rm = TRUE) - data$MMSE)

# === Comparaison d'AIC entre modèle à intercept aléatoire et modèles à intercept et pente aléatoire === #

# Modèle à intercept aléatoire
mi1logMMSE <- lme(fixed = logMMSE ~ OBSTIME + RNFLG + AGE75 + SEXE + DIPNIV + APOE4, data = data, random = ~ 1 | ID, method = "ML", na.action = na.omit)  
summary(mi1logMMSE) # AIC: 125.5908

# Modèle à intercept et pente aléatoires
mip2logMMSE <- lme(fixed = logMMSE ~ OBSTIME + RNFLG + AGE75 + SEXE + DIPNIV + APOE4, data = data, random = ~ OBSTIME | ID, method = "ML", na.action = na.omit)  
summary(mip2logMMSE) # AIC: 109.391

### --- On garde finalement le modèle à intercept et pente aléatoires --- #

### ============ TEST UNIVARIE DES INTERACTIONS ENTRE CHAQUE VARIABLE ET LE DELAI ============== ###

### === RNFLG * OBSTIME === ###
mu1logMMSE <- lme(fixed = logMMSE ~ OBSTIME + RNFLG + AGE75 + SEXE + DIPNIV + APOE4 + RNFLG*OBSTIME, data = data, random = ~ OBSTIME | ID, method = "ML", na.action = na.omit)  
summary(mu1logMMSE) # 0.1430

### === AGE75 * OBSTIME === ###
mu2logMMSE <- lme(fixed = logMMSE ~ OBSTIME + RNFLG + AGE75 + SEXE + DIPNIV + APOE4 + AGE75*OBSTIME, data = data, random = ~ OBSTIME | ID, method = "ML", na.action = na.omit)  
summary(mu2logMMSE) # 0.0058

### === SEXE * OBSTIME === ###
mu3logMMSE <- lme(fixed = logMMSE ~ OBSTIME + RNFLG + AGE75 + SEXE + DIPNIV + APOE4 + SEXE*OBSTIME, data = data, random = ~ OBSTIME | ID, method = "ML", na.action = na.omit)  
summary(mu3logMMSE) # 0.6737

### === DIPNIV * OBSTIME === ###
mu4logMMSE <- lme(fixed = logMMSE ~ OBSTIME + RNFLG + AGE75 + SEXE + DIPNIV + APOE4 + DIPNIV*OBSTIME, data = data, random = ~ OBSTIME | ID, method = "ML", na.action = na.omit)  
summary(mu4logMMSE) # 0.6945, 0.4400

### === APOE4 * OBSTIME === ###
mu5logMMSE <- lme(fixed = logMMSE ~ OBSTIME + RNFLG + AGE75 + SEXE + DIPNIV + APOE4 + APOE4*OBSTIME, data = data, random = ~ OBSTIME | ID, method = "ML", na.action = na.omit)  
summary(mu5logMMSE) # 0.7483

### === MODELE MULTIVARIE === ###
mfinal_logMMSE <- lme(fixed = logMMSE ~ OBSTIME + RNFLG + AGE75 + SEXE + APOE4 + DIPNIV + OBSTIME*(RNFLG + AGE75), data = data, random = ~ OBSTIME | ID, method = "ML", na.action = na.omit)
summary(mfinal_logMMSE)

####### ====== VALIDATION DU MODELE ======== #######
q1logMMSE <- qqnorm(mfinal_logMMSE)
pl1logMMSE <- plot(mfinal_logMMSE)
gridExtra::grid.arrange(q1logMMSE,pl1logMMSE,ncol=2)
