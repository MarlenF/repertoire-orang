## R code for analysis of individual repertoire size (1)

rm(list=ls())
setwd("C:/Users/Marlen Fröhlich/Documents/R")
rep <- read.table ("dataset_ind.rep.csv", header=TRUE, sep=",",stringsAsFactors=TRUE)
library(arm)
library(car)

test.data=rep
test.data$age.dep=as.numeric(test.data$age==levels(test.data$age)[2])
test.data$age.imm=as.numeric(test.data$age==levels(test.data$age)[3])
test.data$sex.code=as.numeric(test.data$sex==levels(test.data$sex)[2])
obs.level.re=as.factor(1:nrow(test.data))
test.data$z.sample=as.vector(scale(test.data$no_sample))


contr=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=10000000))
contr2=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=10000000))


test.data2<- subset(test.data,no_sample>40)
obs.level.re2=as.factor(1:nrow(test.data2))


# collinearity: max vif = 1.3
vif(lm(rnorm(nrow(test.data)) ~  species + setting +  age.dep + age.imm + sex.code, data = test.data))

# run full model (all individuals)
mod.rep = glmer(formula = no_rep ~ species * setting +  age.dep + age.imm +  z.sample + 
                  (0+age.imm|group) +(1|obs.level.re), , 
                  ,family=poisson,  data = test.data, control=contr)

# run reduced model without interaction term (all individuals)
red.rep = glmer(formula = no_rep ~ species + setting +age.dep + age.imm +  z.sample + 
                  + (0+age.imm|group) +(1|obs.level.re), 
                ,family=poisson,  data = test.data, control=contr)

# run null model (all individuals)
null.rep = glmer(formula = no_rep ~   age.dep + age.imm +  z.sample + 
                   + (0+age.imm|group) +(1|obs.level.re), 
                 ,family=poisson,  data = test.data, control=contr)
###################################################################################################3
# run full model (restricted dataset)
mod.rep.asy = glmer(formula = no_rep ~ species * setting +  age.dep + age.imm +  z.sample + 
                      (0+age.imm|group), 
, 
                    ,family=poisson,  data = test.data2, control=contr)

# run reduced model without interaction term (restricted dataset)
red.rep.asy = glmer(formula = no_rep ~ species + setting +  age.dep + age.imm  +  z.sample + 
                      + (0+age.imm|group), 
                    ,family=poisson,  data = test.data2, control=contr)

# run null model (restricted dataset)
null.rep.asy = glmer(formula = no_rep ~  age.dep + age.imm  + z.sample  + 
                       + (0+age.imm|group), 
                     ,family=poisson,  data = test.data2, control=contr)


setwd("C:/Users/Marlen Fröhlich/Documents/R/roger/")
source("diagnostic_fcns.r")
# full dataset: dispersion parameter 0.88
overdisp.test(mod.rep)

# restricted dataset: dispersion parameter 0.577
overdisp.test(mod.rep.asy)



length(residuals(mod.rep)) #71
length(residuals(mod.rep.asy)) #51
length(residuals(null.rep)) #71
length(residuals(null.rep.asy)) #51

#Likelihood ratio test, all individuals
as.data.frame(anova(null.rep, mod.rep, test="Chisq"))


#Likelihood ratio test, restricted
as.data.frame(anova(null.rep.asy, mod.rep.asy, test="Chisq"))


#get chi square and p-values (all)
drop1(mod.rep,  test ="Chisq")
  

#get chi square and p-values (all)
drop1(red.rep,  test ="Chisq")


#get chi square and p-values (restricted)
drop1(mod.rep.asy,  test ="Chisq")


drop1(red.rep.asy,  test ="Chisq")



#get estimates and SE (all individuals)
round(summary(red.rep)$coefficients, 3)


#get estimates and SE (restricted)
round(summary(red.rep.asy)$coefficients, 3)



#dummy-code factors
test.data$age.dep.c= test.data$age.dep- mean(test.data$age.dep)
test.data$age.imm.c= test.data$age.imm- mean(test.data$age.imm)
test.data$sex.code.c= test.data$sex.code- mean(test.data$sex.code)
test.data$z.sample=as.vector(scale(test.data$no_sample))

test.data2$age.dep.c= test.data2$age.dep- mean(test.data2$age.dep)
test.data2$age.imm.c= test.data2$age.imm- mean(test.data2$age.imm)
#test.data2$sex.code.c= test.data2$sex.code- mean(test.data2$sex.code)
test.data2$z.sample=as.vector(scale(test.data2$no_sample))



plot.mod.rep = glmer(formula = no_rep ~ species + setting +  age.dep.c + age.imm.c +  z.sample +
                       + (0+age.imm|group) +(1|obs.level.re), 
                     ,family=poisson,  data = test.data, control=contr)

plot.mod.asy = glmer(formula = no_rep ~ species + setting +  age.dep.c + age.imm.c + z.sample +
                       + (0+age.imm|group), 
                     ,family=poisson,  data = test.data2, control=contr)


#####################plot setting and age differences#################################3
plot.age=aggregate(x=test.data$no_rep, by=test.data[, c("species", "setting","age")], FUN=mean)


library(ggplot2)
theme_marlen_ss <- theme(panel.background = element_blank(),
                         panel.border =element_rect(colour="black", fill=NA),
                         
                         plot.background = element_blank(),
                         panel.grid = element_blank(),
                         axis.line = element_line(colour ="black"),
                         axis.text.x = element_text (size = 12,colour= "black", family="sans"),
                         axis.text.y = element_text (size = 12,colour= "black", family="sans"),
                         axis.ticks.y = element_line(colour="black"),
                         axis.ticks.x = element_line(colour=NA),
                         axis.title.x = element_text(size = 12, vjust = -0.5, family="sans"),
                         axis.title.y = element_text(size = 12, vjust = 2, family="sans"),
                         legend.text=  element_text(size = 11, family="sans", margin = margin(t = 10)),
                         legend.key = element_blank(),
                         legend.position = "right",
                         legend.spacing.x = unit(0.2, 'cm'),
                         strip.text = element_text(size = 11))

levels(plot.age$species) <- c("Bornean", "Sumatran")

test.data2$age<- factor(test.data2$age, levels = c("Ad", "Im", "Dp"))



#modt$partner <- factor(suc$partner, levels = c("M", "K", "C"))
dodge.posn <- position_dodge(.9)
#svg("/Users/mfroehlich/R/self.svg", width=90/25.4, height=90/25.4, pointsize=9,family="arial")
mod.site <- ggplot(test.data2, aes(x = setting, y = no_rep))
mod.site + geom_boxplot(aes(fill = age), width = 0.9) +
  geom_point(aes(fill = age), width = 1,position= dodge.posn, shape = 1, colour = "black", alpha = 0.5) +
  theme_marlen_ss +
  scale_y_continuous("Individual repertoire size") +
  scale_x_discrete("Setting",
                   limits = c("captive", "wild"),
                   labels = c("captive", "wild"))+
  scale_fill_manual(values=c("blue", "orange", "grey"),name="Age class",
                    breaks=c("Ad","Im", "Dp"),
                    labels=c("Adult", "Older immature", "Younger immature"))+
  facet_wrap(~species)+
  stat_summary(fun=mean, geom="point",shape =23, fill ="black",aes(group=age), position=position_dodge(.9), 
               color="black", size=3)
