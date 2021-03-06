rm(list=ls())
library(lme4)
library(lmerTest)
library(car)


# step 0: define path and read data ---
path <- "C:\\Users\\uknief\\Documents\\CoAuthoredManuscripts\\OrangUtan_RepertoireSimilarity\\DiceCoefficient\\"
motoff <- read.table(paste0(path,"dataset_complete.csv"), header=TRUE, sep=",")



colnames(motoff)[which(colnames(motoff)=="ID_sign")] <- "animal_id"
colnames(motoff)[which(colnames(motoff)=="Signal")] <- "Behavior2"
motoff <- as.data.frame(na.omit(motoff[ ,c("Behavior2","animal_id","species","group","setting","ID_obs","class")]))


# step 1: omit all levels of animal_id that contributed fewer than 40 cases ---
NID <- data.frame(table(motoff$animal_id)) ; colnames(NID) <- c("animal_id","N")
NID$animal_id <- as.character(NID$animal_id)
NID40 <- subset(NID, N>=40)
motoff <- subset(motoff, animal_id %in% NID40$animal_id)


# step 2: subset data (mother/offspring repertoire) ---
#dat.ad <- subset(motoff, class=="Mot" & species=="Bor")
#dat.ad <- subset(motoff, class=="Mot" & species=="Sum")
#dat.ad <- subset(motoff, class=="Inf" & species=="Bor")
dat.ad <- subset(motoff, class=="Inf" & species=="Sum")


# step 3: create a table for the set of behaviors (Behavior2) used by every individual (animal_ID) AT LEAST TWICE, let's call this R ---
## create data frame in which every combination of individual and behaviour is counted
dat.ind.bhv <- data.frame(table(dat.ad$animal_id,dat.ad$Behavior2)) ; colnames(dat.ind.bhv) <- c("animal_id","Behavior2","N")
## subset dat.ind.bhv to include only rows where N>=2
dat.ind.bhv2 <- subset(dat.ind.bhv, N>=2)
## create data frame with the number of behaviours for each individual
dat.ind.bhv2.n <- data.frame(table(dat.ind.bhv2$animal_id)) ; colnames(dat.ind.bhv2.n) <- c("animal_id","N")


# step 4: calculate Dice coefficient based on formula: dc = (2 x number of behaviours two inds have in common)/(R of ind1 + R ind2) ---
## create data frame with animal_id as rows and columns
dat.ind.ind <- data.frame(matrix(rep(NA, nrow(dat.ind.bhv2.n) * nrow(dat.ind.bhv2.n)), ncol=nrow(dat.ind.bhv2.n)))
colnames(dat.ind.ind) <- dat.ind.bhv2.n$animal_id
rownames(dat.ind.ind) <- dat.ind.bhv2.n$animal_id
## loop over all pairs of individuals
for(i in 1:(nrow(dat.ind.bhv2.n))) {
	for(j in 1:(nrow(dat.ind.bhv2.n))) {
		Ind1 <- subset(dat.ind.bhv2, animal_id==dat.ind.bhv2.n$animal_id[i])
		Ind2 <- subset(dat.ind.bhv2, animal_id==dat.ind.bhv2.n$animal_id[j])
		ovrlp <- length(which(Ind1$Behavior2 %in% Ind2$Behavior2))
		RInd1 <- subset(dat.ind.bhv2.n, animal_id==dat.ind.bhv2.n$animal_id[i])$N	# which is the same as nrow(Ind1)
		RInd2 <- subset(dat.ind.bhv2.n, animal_id==dat.ind.bhv2.n$animal_id[j])$N	# which is the same as nrow(Ind2)
		dat.ind.ind[which(rownames(dat.ind.ind)==dat.ind.bhv2.n$animal_id[i]),which(colnames(dat.ind.ind)==dat.ind.bhv2.n$animal_id[j])] <- (2 * ovrlp) / (RInd1 + RInd2)
	}
}
## set diagonal to NA
Dice <- as.matrix(dat.ind.ind)
diag(Dice) <- NA
#Dice <- as.data.frame(Dice)


# step 5: define within and between setting diads ---
## create data frame with animal_id as rows and columns
dat.wth.btw <- data.frame(matrix(rep(NA, nrow(dat.ind.bhv2.n) * nrow(dat.ind.bhv2.n)), ncol=nrow(dat.ind.bhv2.n)))
colnames(dat.wth.btw) <- dat.ind.bhv2.n$animal_id
rownames(dat.wth.btw) <- dat.ind.bhv2.n$animal_id
## loop over all pairs of individuals
for(i in 1:(nrow(dat.ind.bhv2.n))) {
	for(j in 1:(nrow(dat.ind.bhv2.n))) {
		Ind1 <- subset(motoff, animal_id==dat.ind.bhv2.n$animal_id[i])
		Ind2 <- subset(motoff, animal_id==dat.ind.bhv2.n$animal_id[j])
		if(length(table(Ind1$setting))>1) { print("Ohoh, better check my data again!") }
		if(length(table(Ind2$setting))>1) { print("Ohoh, better check my data again!") }
		if(Ind1$setting[1]==Ind2$setting[1]) { wth.btw <- "W" } else { wth.btw <- "B" }
		dat.wth.btw[which(rownames(dat.wth.btw)==dat.ind.bhv2.n$animal_id[i]),which(colnames(dat.wth.btw)==dat.ind.bhv2.n$animal_id[j])] <- wth.btw
	}
}
## set diagonal to NA
wth.btw <- as.matrix(dat.wth.btw)
diag(wth.btw) <- NA
#wth.btw <- as.data.frame(wth.btw)


# step 6: select within and between group diads from the Dice data frame ---
Dice.Within <- ifelse(wth.btw=="W",Dice,NA)
Dice.Between <- ifelse(wth.btw=="B",Dice,NA)

(Emp.MDiceW <- mean(Dice.Within[lower.tri(Dice.Within, diag = FALSE)],na.rm=TRUE))
(Emp.MDiceB <- mean(Dice.Between[lower.tri(Dice.Between, diag = FALSE)],na.rm=TRUE))

# Bornean moms within: 0.7312891
# Bornean moms between: 0.5723846
# Bornean infants within: 0.7635561
# Bornean infants between: 0.6107028
# Sumatran moms within: 0.6217871
# Sumatran moms between: 0.5301574
# Sumatran infants within: 0.7277381
# Sumatran infants between: 0.6118929

## create data frame dd with columns Ind1, Ind2, Dice, Within/Between
yW <- expand.grid(rownames(Dice.Within), colnames(Dice.Within))
labs <- yW[as.vector(upper.tri(Dice.Within, diag = FALSE)), ]
yW <- cbind(labs, Dice.Within[upper.tri(Dice.Within,diag=FALSE)])
colnames(yW) <- c("Ind1","Ind2","Dice")
yW <- yW[!is.na(yW$Dice), ]
yW$WB <- rep("Within",nrow(yW))

yB <- expand.grid(rownames(Dice.Between), colnames(Dice.Between))
labs <- yB[as.vector(upper.tri(Dice.Between, diag = FALSE)), ]
yB <- cbind(labs, Dice.Between[upper.tri(Dice.Between,diag=FALSE)])
colnames(yB) <- c("Ind1","Ind2","Dice")
yB <- yB[!is.na(yB$Dice), ]
yB$WB <- rep("Between",nrow(yB))
dd <- rbind(yW,yB)
motoff.sub <- motoff[,c("animal_id","species","group","setting")]
motoff.sub <- motoff.sub[!duplicated(motoff.sub), ]
dd <- merge(dd, motoff.sub, by.x="Ind1", by.y="animal_id", sort=FALSE, all.y=FALSE)
dd <- merge(dd, motoff.sub, by.x="Ind2", by.y="animal_id", sort=FALSE)
colnames(dd)[which(colnames(dd)=="species.x")] <- "species.ind1"
colnames(dd)[which(colnames(dd)=="group.x")] <- "group.ind1"
colnames(dd)[which(colnames(dd)=="setting.x")] <- "setting.ind1"
colnames(dd)[which(colnames(dd)=="species.y")] <- "species.ind2"
colnames(dd)[which(colnames(dd)=="group.y")] <- "group.ind2"
colnames(dd)[which(colnames(dd)=="setting.y")] <- "setting.ind2"


# step 7: matrix permutation test ---
motoff.inc <- subset(motoff, animal_id %in% dat.ind.bhv2.n$animal_id)
ind.setting <- unique(motoff.inc[,c("animal_id","group","species","setting")])
ind.setting$species_setting <- paste(ind.setting$species, ind.setting$setting, sep="_")

nPerm <- 1000
dat.Perm <- c()
for(k in 1:nPerm) {
	## permute setting in motoff
	ind.setting.perm <- ind.setting
	ind.setting.perm$setting <- sample(ind.setting.perm$setting, nrow(ind.setting.perm), replace=FALSE)
	motoff.inc.perm <- motoff.inc[,-which(colnames(motoff.inc)=="setting")]
	motoff.inc.perm <- merge(motoff.inc.perm, ind.setting.perm, by.x="animal_id", by.y="animal_id", sort=FALSE)
	
	## permute setting blockwise in species_setting
#	ind.setting.perm <- ind.setting
#	for(l in 1:length(unique(ind.setting$species_setting))) {
#		pos <- which(ind.setting$species_setting %in% unique(ind.setting$species_setting)[l])
#		ind.setting.perm$setting[pos] <- sample(ind.setting.perm$setting[pos], length(pos), replace=FALSE)
#	}
#	motoff.inc.perm <- motoff.inc[,-which(colnames(motoff.inc) %in% c("group","species","setting"))]
#	motoff.inc.perm <- merge(motoff.inc.perm, ind.setting.perm, by.x="animal_id", by.y="animal_id", sort=FALSE)	
	
	## repeat step 5 & 6
	## create data frame with animal_id as rows and columns
	dat.wth.btw <- data.frame(matrix(rep(NA, nrow(dat.ind.bhv2.n) * nrow(dat.ind.bhv2.n)), ncol=nrow(dat.ind.bhv2.n)))
	colnames(dat.wth.btw) <- dat.ind.bhv2.n$animal_id
	rownames(dat.wth.btw) <- dat.ind.bhv2.n$animal_id
	## loop over all pairs of individuals
	for(i in 1:(nrow(dat.ind.bhv2.n))) {
		for(j in 1:(nrow(dat.ind.bhv2.n))) {
			Ind1 <- subset(motoff.inc.perm, animal_id==dat.ind.bhv2.n$animal_id[i])
			Ind2 <- subset(motoff.inc.perm, animal_id==dat.ind.bhv2.n$animal_id[j])
			if(length(table(Ind1$setting))>1) { print("Ohoh, better check my data again!") }
			if(length(table(Ind2$setting))>1) { print("Ohoh, better check my data again!") }
			if(Ind1$setting[1]==Ind2$setting[1]) { wth.btw <- "W" } else { wth.btw <- "B" }
			dat.wth.btw[which(rownames(dat.wth.btw)==dat.ind.bhv2.n$animal_id[i]),which(colnames(dat.wth.btw)==dat.ind.bhv2.n$animal_id[j])] <- wth.btw
		}
	}
	## set diagonal to NA
	wth.btw <- as.matrix(dat.wth.btw)
	diag(wth.btw) <- NA
	Dice.Within <- ifelse(wth.btw=="W",Dice,NA)
	Dice.Between <- ifelse(wth.btw=="B",Dice,NA)
	MDiceW <- mean(Dice.Within[lower.tri(Dice.Within, diag = FALSE)],na.rm=TRUE)
	MDiceB <- mean(Dice.Between[lower.tri(Dice.Between, diag = FALSE)],na.rm=TRUE)
	dat.Perm[k] <- MDiceW - MDiceB
	flush.console()
	if(k %% 10 == 0) { print(paste0("Finished ", k, " out of ", nPerm, " simulations")) }
}

hist(dat.Perm)
abline(v=Emp.MDiceW - Emp.MDiceB, col="red")
1 - sum(dat.Perm<=(Emp.MDiceW - Emp.MDiceB)) / length(dat.Perm)	# This is the P-value

# Setting Bornean moms Within/Between: P = 0
# Setting Bornean infants Within/Between: P = 0
# Setting Sumatran moms Within/Between: P = 0.018
# Setting Sumatran infants Within/Between: P = 0

# Boxplot Dice within between
path <- "C:\\Users\\Marlen Fr�hlich\\Documents\\R\\"
dice <- read.table(paste0(path,"data_dice.coeff.csv"), header=TRUE, sep=",")

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

levels(dice$species.ind1) <- c("Bornean", "Sumatran")

dodge.posn <- position_dodge(.9)
#svg("/Users/mfroehlich/R/self.svg", width=90/25.4, height=90/25.4, pointsize=9,family="arial")
mod.site <- ggplot(dice, aes(x = class, y = Dice))
mod.site + geom_boxplot(aes(fill = WB), width = 0.9) +
  geom_point(aes(fill = WB),position= dodge.posn, shape = 1, colour = "black", alpha = 0.5) +
  theme_marlen_ss +
  scale_y_continuous("Repertoire similarity (Dice coefficients)") +
  scale_x_discrete("Age class",
                   limits = c("Infs", "Moms"),
                   labels = c("Infants", "Mothers"))+
  scale_fill_manual(values=c("blue", "orange"),name="Setting",
                    breaks=c("Between","Within"),
                    labels=c("Between","Within"))+
  facet_wrap(~species.ind1)+
  stat_summary(fun=mean, geom="point",shape =23, fill ="black",aes(group=WB), position=position_dodge(.9), 
               color="black", size=3)



