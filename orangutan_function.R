rm(list=ls())
setwd("C:/Users/Marlen Fröhlich/Documents/R")
fun <- read.table ("orangutan_fun.csv", header=TRUE, sep=",", stringsAsFactors=TRUE)
xx=as.data.frame(na.omit(fun[, c("species","signal", "setting", "dominance_con", "context", "no_cases", "no_subjects")]))

library(arm)
library(car)

test.data=xx

test.data$z.subjects=as.vector(scale(test.data$no_subjects))
test.data$z.cases=as.vector(scale(test.data$no_cases))
test.data$context_pl=as.numeric(test.data$context==levels(test.data$context)[4])
contr=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=10000000))

# collinearity: max vif = 2.75 
vif(lm(rnorm(nrow(test.data)) ~  species + setting + log(no_cases) + log(no_subjects) , data = test.data))

#run the full model
mod.fun = lmer(formula = dominance_con ~ species * setting +  context_pl + log(no_cases) + log(no_subjects) +
               +(1|signal),  
                 data = test.data)
#run the null model
null.fun = lmer(formula = dominance_con ~ log(no_cases) + log(no_subjects) +
                 +(1|signal),  
               data = test.data)

length(residuals(mod.spec)) #114
length(residuals(null.spec)) #114

#Likelihood ratio test
as.data.frame(anova(null.fun, mod.fun, test="Chisq"))
#npar       AIC       BIC   logLik   deviance    Chisq Df   Pr(>Chisq)
#null.fun    5 -64.00926 -50.32827 37.00463  -74.00926       NA NA           NA
#mod.fun     9 -82.20468 -57.57890 50.10234 -100.20468 26.19542  4 2.889918e-05

# get chi square and p-values
drop1(mod.fun,  test ="Chisq")
#context_pl          1 -78.294  5.9106 0.015050 * 
#log(no_cases)       1 -84.162  0.0428 0.836111   
#log(no_subjects)    1 -79.627  4.5775 0.032395 * 
#species:setting     1 -73.810 10.3945 0.001264 **

#get estimates and SE
round(summary(mod.fun)$coefficients, 3)
#Estimate Std. Error t value
#(Intercept)               0.898      0.060  14.980
#speciesSum                0.002      0.039   0.060
#settingwild               0.059      0.043   1.367
#context_pl                0.095      0.038   2.499
#log(no_cases)             0.004      0.019   0.197
#log(no_subjects)         -0.072      0.034  -2.093
#speciesSum:settingwild   -0.175      0.055  -3.196


#####################plot site differences in monopol#################################3

# dummycoding for plot
context_pl.c= test.data$context_pl- mean(test.data$context_pl)

######################################################
#plot model
plot.fun = glmer(formula = dominance_con ~ species * setting +  context_pl.c + z.cases + z.subjects +
                   +(1|signal),  
                 data = test.data)


