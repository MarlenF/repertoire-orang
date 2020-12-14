## R code for analysis of individual repertoire sizes

rm(list=ls())
setwd("C:/Users/Marlen Fröhlich/Documents/R")
rep <- read.table ("orangutan_indrep.csv", header=TRUE, sep=",",stringsAsFactors=TRUE)
library(arm)
library(car)

test.data=rep
test.data$age.dep=as.numeric(test.data$age==levels(test.data$age)[2])
test.data$age.imm=as.numeric(test.data$age==levels(test.data$age)[3])
test.data$sex.code=as.numeric(test.data$sex==levels(test.data$sex)[2])
test.data$z.context=as.vector(scale(test.data$no_context))
test.data$z.sample=as.vector(scale(test.data$no_sample))

contr=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=10000000))

# collinearity: max vif = 2.85
vif(lm(rnorm(nrow(test.data)) ~  species + setting +  age.dep + age.imm + sex.code + no_context + log(no_sample), data = test.data))

# run full model (all individuals)
mod.rep = lmer(formula = log(no_rep) ~ species * setting +  age.dep + age.imm + sex.code + no_context + log(no_sample) +
                  (0+age.dep|group)  + (0+age.imm|group)+ (0+sex.code|group) + #(0+ no_context|group) + 
                  (1|group),  data = test.data)

# run reduced model without interaction term (all individuals)
red.rep = lmer(formula = log(no_rep) ~ species + setting +  age.dep + age.imm + sex.code + no_context + log(no_sample) +
                 (0+age.dep|group)  + (0+age.imm|group)+ (0+sex.code|group) + #(0+ no_context|group) + 
                 (1|group),  data = test.data)

# run null model (all individuals)
null.rep = lmer(formula = log(no_rep) ~  age.dep + age.imm + sex.code + no_context + log(no_sample) +
                 (0+age.dep|group)  + (0+age.imm|group)+ (0+sex.code|group) + #(0+ no_context|group)  + 
                 (1|group),  data = test.data)
###################################################################################################3
# run full model (restricted dataset)
mod.rep.asy = lmer(formula = log(no_rep_asy) ~ species * setting +  age.dep + age.imm + sex.code + no_context + log(no_sample) +
                 (0+age.dep|group)  + (0+age.imm|group)+ (0+sex.code|group) + #(0+ no_context|group) + 
                 (1|group),  data = test.data)

# run reduced model without interaction term (restricted dataset)
red.rep.asy = lmer(formula = log(no_rep_asy) ~ species + setting +  age.dep + age.imm + sex.code + no_context + log(no_sample) +
                     (0+age.dep|group)  + (0+age.imm|group)+ (0+sex.code|group) + #(0+ no_context|group) + 
                     (1|group),  data = test.data)

# run null model (restricted dataset)
null.rep.asy = lmer(formula = log(no_rep_asy) ~  age.dep + age.imm + sex.code + no_context + log(no_sample) +
                  (0+age.dep|group)  + (0+age.imm|group)+ (0+sex.code|group) + #(0+ no_context|group)  + 
                  (1|group),  data = test.data)

length(residuals(mod.rep)) #70
length(residuals(mod.rep.asy)) #44
length(residuals(null.rep)) #70
length(residuals(null.rep.asy)) #44

#Likelihood ratio test, all individuals
as.data.frame(anova(null.rep, mod.rep, test="Chisq"))

#Df      AIC      BIC   logLik deviance   Chisq Chi Df  Pr(>Chisq)
#null.rep 11 39.23088 63.96433 -8.61544 17.23088      NA     NA          NA
#mod.rep  14 31.17298 62.65191 -1.58649  3.17298 14.0579      3 0.002827391

#Likelihood ratio test, restricted
as.data.frame(anova(null.rep.asy, mod.rep.asy, test="Chisq"))

#Df        AIC      BIC   logLik  deviance    Chisq Chi Df  Pr(>Chisq)
#null.rep 11  -5.734699 13.89139 13.86735 -27.73470       NA     NA          NA
#mod.rep  14 -14.865271 10.11338 21.43264 -42.86527 15.13057      3 0.001708398

#get chi square and p-values (all)
drop1(red.rep,  test ="Chisq")
#species           1  28.505  1.325 0.2496437    
#setting           1  40.959 13.780 0.0002055 ***
#age.dep           1  30.872  3.692 0.0546636 .  
#age.imm           1  30.004  2.825 0.0928163 .  
#sex.code          1  27.782  0.603 0.4374396    
#no_context        1  37.752 10.573 0.0011475 ** 
#log(no_sample)    1 109.951 82.771 < 2.2e-16 ***

#get chi square and p-values (restricted)
drop1(red.rep.asy,  test ="Chisq")
#species           1 -12.0522  6.641 0.0099641 ** 
#setting           1  -5.9661 12.727 0.0003603 ***
#age.dep           1  -4.1226 14.571 0.0001350 ***
#age.imm           1 -15.2972  3.396 0.0653416 .  
#sex.code          1 -18.4706  0.223 0.6368329    
#no_context        1 -15.0392  3.654 0.0559256 .  
#log(no_sample)    1  20.7324 39.426 3.407e-10 ***


#get estimates and SE (all individuals)
round(summary(red.rep)$coefficients, 3)

#(Intercept)       -0.578      0.151  -3.827
#speciesSumatran    0.070      0.066   1.059
#settingwild       -0.353      0.067  -5.234
#age.dep            0.170      0.090   1.891
#age.imm            0.188      0.114   1.655
#sex.code          -0.059      0.074  -0.801
#no_context         0.114      0.036   3.212
#log(no_sample)     0.557      0.047  11.840

#get estimates and SE (restricted)
round(summary(red.rep.asy)$coefficients, 3)

#(Intercept)        0.453      0.227   1.998
#speciesSumatran    0.148      0.061   2.404
#settingwild       -0.290      0.069  -4.226
#age.dep            0.314      0.076   4.112
#age.imm            0.243      0.138   1.753
#sex.code          -0.039      0.064  -0.600
#no_context         0.058      0.030   1.941
#log(no_sample)     0.361      0.050   7.221

#dummy-code factors
test.data$age.dep.c= test.data$age.dep- mean(test.data$age.dep)
test.data$age.imm.c= test.data$age.imm- mean(test.data$age.imm)
test.data$sex.code.c= test.data$sex.code- mean(test.data$sex.code)
