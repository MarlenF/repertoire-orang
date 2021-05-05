## R code for analysis of functional specificity (3)

rm(list=ls())
setwd("C:/Users/Marlen Fröhlich/Documents/R")
fun <- read.table ("dataset_func.spec.csv", header=TRUE, sep=",", stringsAsFactors=TRUE)
xx=as.data.frame(na.omit(fun[, c("species","signal", "setting","dominance_con", "max_cases", "context", "no_cases", "no_subjects")]))

library(arm)
library(car)

test.data=xx

test.data$z.subjects=as.vector(scale(test.data$no_subjects))
test.data$z.sample=as.vector(scale(test.data$no_cases))

test.data$context_pl=as.numeric(test.data$context==levels(test.data$context)[4])
contr=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=10000000))
dom_bin = cbind(test.data$max_cases, test.data$no_cases-test.data$max_cases)

contr=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=10000000))



# collinearity: max vif = 2.7 
vif(lm(rnorm(nrow(test.data)) ~  species + setting + log(no_cases) + log(no_subjects) , data = test.data))

#run the full model
mod.fun = glmer(formula = dom_bin ~ species * setting +  context_pl +  z.subjects + z.sample +
               +(1|signal),  family =binomial, control=contr,
                 data = test.data)
#run the null model
null.fun = glmer(formula = dom_bin ~ context_pl +  z.subjects + z.sample +
                  +(1|signal),  family =binomial, control=contr,
                data = test.data)
              

length(residuals(mod.fun)) #128
length(residuals(null.fun)) #128

#Likelihood ratio test
as.data.frame(anova(null.fun, mod.fun, test="Chisq"))


# get chi square and p-values
drop1(mod.fun,  test ="Chisq")

#post hoc sidak test of interaction effect
require(lsmeans)
lsm <- lsmeans(mod.fun, list(pairwise ~ setting|species, pairwise ~ species|setting))
lsm

  
#get estimates and SE
round(summary(mod.fun)$coefficients, 3)
####################3

# dummycoding for plot
context_pl.c= test.data$context_pl- mean(test.data$context_pl)

######################################################
path <- "C:/Users/Marlen Fröhlich/Documents/R/"


#plot model
plot.fun = glmer(formula = dom_bin ~ species * setting +  context_pl.c +  z.subjects + z.sample +
                   +(1|signal),  family =binomial, control=contr,
                 data = test.data)

setwd("C:/Users/Marlen Fröhlich/Documents/R/roger/")
source("local_jitter_AlexHausmann.R")
require(effects)

test.data$XPos <- ifelse(test.data$species=="Bor",1,2)
EF <- Effect(c("species","setting"),plot.fun,se=TRUE)
dat.EF <- as.data.frame(EF)

# Add colour column
test.data$colourBG <- ifelse(test.data$setting=="wild",rgb(255, 210, 128, maxColorValue=255),rgb(128, 128, 255, maxColorValue=255))
test.data$colourL <- ifelse(test.data$setting=="wild",rgb(255, 192, 77, maxColorValue=255),rgb(77, 77, 255, maxColorValue=255))

# Open empty plot (IMPORTANT: THE PLOT HAS TO BE OPEN BEFORE RUNNING THE FUNCTION)
path <- "C:\\Users\\Marlen Fröhlich\\Documents\\R\\"


svg(filename=paste0(path,"FunSpecN2.svg",sep=""), height=90/25.4, width=90/25.4, family="Arial", pointsize=9)
OF <- 0.1
par(mar=c(2.7, 3.2, 0.2, 0.2), mgp=c(1.3, 0.2, 0), tcl=-0.25, cex=1)
plot(c(0.5,2.5),c(0.35,1), type="n", axes=FALSE, xlab="Orang-utan species", ylab="") ; par(new=TRUE)
plot(c(0.5,2.5),c(0.35,1), type="n", axes=FALSE, xlab="", ylab="Functional specificity", mgp=c(2.2, 0.2, 0))

X0 <-local_jitter(fact_coord = test.data$XPos, gradual_coord = test.data$dominance_con, categories = as.character(test.data$setting), factorial_axis = 1, buffer = 0.45, sizes = sqrt(test.data$no_cases)/4, verbose=F, iterations=1000)
points(X0,test.data$dominance_con,cex=sqrt(test.data$no_cases)/4, pch=21, bg=test.data$colourBG, col=test.data$colourL)

arrows(1-OF,dat.EF$lower[1],1-OF,dat.EF$upper[1],code=3,length=0.1,angle=90)
points(x=1-OF,y=dat.EF$fit[1], pch=23, col="black", bg="blue", cex=3)

arrows(1+OF,dat.EF$lower[3],1+OF,dat.EF$upper[3],code=3,length=0.1,angle=90)
points(x=1+OF,y=dat.EF$fit[3], pch=23, col="black", bg="orange", cex=3)

arrows(2-OF,dat.EF$lower[2],2-OF,dat.EF$upper[2],code=3,length=0.1,angle=90)
points(x=2-OF,y=dat.EF$fit[2], pch=23, col="black", bg="blue", cex=3)

arrows(2+OF,dat.EF$lower[4],2+OF,dat.EF$upper[4],code=3,length=0.1,angle=90)
points(x=2+OF,y=dat.EF$fit[4], pch=23, col="black", bg="orange", cex=3)

axis(1,at=c(1,2), label=c("Bornean","Sumatran"), tcl=-0.25)
axis(2,at=seq(0,1,by=0.2), label=c("0.0","0.2","0.4","0.6","0.8","1.0"), tcl=-0.25, las=2, mgp=c(1.2, 0.4, 0))
legend("topright", pt.bg=c("blue","orange"), pch=23, legend=c("captive","wild"), bty="n", pt.cex=2)
box()
dev.off()
#border=c(rgb(77, 77, 255, maxColorValue=255),rgb(255, 192, 77, maxColorValue=255)), fill=c(rgb(128, 128, 255, maxColorValue=255), rgb(255, 210, 128, maxColorValue=255))



domrep=aggregate(x=test.data$dominance_con, by=test.data[, c("species", "setting")], FUN=mean)
domrep


