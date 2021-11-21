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

# ========== Importation et construction de la base  ========== #

# ____ Définir le répertoire ____ #
setwd("/home/tuvu/Downloads/projet_STA302/projet")

# ____ Importation ____ #

df <- read.csv("projet.csv", sep=";", dec = ".", na.strings = NA)

# ____ Types de variables ____ #

#df$SEXE <- as.factor(df$SEXE)
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

# ____ Cr?ation de nouvelles variables ____ #

df$RNFLG <- ifelse(!is.na(df$RNFLGD_1), df$RNFLGD_1, df$RNFLGG_1) # RNFLG droit si pas NA, RNGLG gauche sinon
df$AX <- ifelse(!is.na(df$RNFLGD_1), df$LG_AX_OD, df$LG_AX_OG) # AXIALE droit si RNFLG droit pas NA, AXIALE gauche sinon

# for (i in c(1:nrow(df))){
#   if(!is.na(df$RNFLGD_1)){
#     df$RNFLG <- df$RNFLGD_1
#     df$AXE <- df$LG_AX_OD
#   }
#   else{
#     df$RNFLG <- df$RNFLGG_1
#     df$AXE <- df$LG_AX_OG
#   }

# ____ Cr?ation du jeu de donn?es de travail "base" ____ #

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

# Cr?ation de la variable de d?lai
data$OBSTIME <- (data$AGE - data$AGEINT)/10


# ========== V?rifications ========== #

# Donnees manquantes dans "data"
sapply(data, function(x) sum(is.na(x)))
sum(is.na(data))

# ========== Analyse descriptive ========== #

summary(data)

# Nombre d'observations
df$ID %>% unique() %>% length()

str(data)

# ========== Graphiques ========== #

# pMMSE <- (ggplot(data)
#           + geom_line(aes(x = OBSTIME, y = MMSE, group = ID), color="grey30", alpha = 0.8)
#           + stat_smooth(aes(x = OBSTIME, y = MMSE), method = "loess", size = 0.75)
#           + theme_bw()
#           + xlab("Temps depuis l'entrée dans l'étude (en dizaines d'années)")
#           + ylab("MMSE")
# )
# pMMSE

pISAAC <- (ggplot(data)
           + geom_line(aes(x = OBSTIME, y = ISAAC, group = ID), color="grey30", alpha = 0.8)
           + stat_smooth(aes(x = OBSTIME, y = ISAAC), method = "loess", size = 0.75)
           + theme_bw()
           + xlab("Temps depuis l'entrée dans l'étude (en dizaines d'années)")
           + ylab("ISAAC")
)
pISAAC

# pBENTON <- (ggplot(data)
#             + geom_line(aes(x = OBSTIME, y = BENTON, group = ID), color="grey30", alpha = 0.8)
#             + stat_smooth(aes(x = OBSTIME, y = BENTON), method = "loess", size = 0.75)
#             + theme_bw()
#             + xlab("Temps depuis l'entrée dans l'étude (en dizaines d'années)")
#             + ylab("BENTON")
# )
# pBENTON


# ========== Mod?les ========== #

# Comparaison des AIC
# c(AIC(m2), AIC(m1))

# Test de l'ajout de la pente
# devm1m2 <- 2*logLik(m1) - 2*logLik(m2)
# pm1m2 <- 0.5*(1-pchisq(devm1m2,df=2)) + 0.5*(1-pchisq(devm1m2,df=1))
# pm1m2

m1 <- lme(fixed = ISAAC ~ OBSTIME + SEXE + DIPNIV + AGEINT + APOE4 + RNFLG, data = data, random = ~ OBSTIME | ID, method = "ML", na.action = na.omit)
summary(m1)

m2 <- lme(fixed = ISAAC ~ OBSTIME + SEXE, data = data, random = ~ OBSTIME | ID, method = "ML", na.action = na.omit)
summary(m2)

m3 <- lme(fixed = ISAAC ~ OBSTIME + SEXE + AGEINT + I(RNFLG/10), data = data, random = ~ OBSTIME + I(RNFLG/10)| ID, method = "ML", na.action = na.omit)
summary(m3)

data$AGE_MEAN <- (data$AGE-mean(data$AGE, na.rm=TRUE))/10

m4 <- lme(fixed = ISAAC ~ SEXE + AGEINT + RNFLG + DIPNIV, data = data, random = ~ AGE_MEAN + I(AGE_MEAN^2) | ID, method = "ML", na.action = na.omit)
summary(m4)


#m4 <- lcmm::hlme(fixed = ISAAC ~ SEXE + AGEINT + RNFLG + DIPNIV, data = data, random = ~ AGE_MEAN + I(AGE_MEAN^2), subject = "ID", na.action = na.omit)
#summary(m4)
