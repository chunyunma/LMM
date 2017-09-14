library(pacman)
p_load(tcltk2,lme4,arm,tidyr,plyr,ggplot2,reshape2, lattice, GGally)
tclTaskSchedule(1000, {options(prompt=paste(Sys.time(),"> "))}, id = "ticktock", redo = TRUE)

DT3 <- read.csv(file.choose())
# pseudo_rawdata.csv


# data prep -------------------------
## rename columns, make them the same as Liu97
names(DT3)[2] <- 'sub'
names(DT3)[10] <- 'Lvalue'
names(DT3)[11] <- 'Rvalue'
names(DT3)[13] <- 'psize'  # this psize is product for mul and sum for add
names(DT3)[14] <- 'accuracy_obsolet'  # this column seems have been replaced by acc.recoded, and I would like to name acc.recoded as accuracy
names(DT3)[19] <- 'correct'
DT3$psizex <- DT3$Lvalue * DT3$Rvalue

### column orders
# 4 operation, 5 format, 10 left_digit, 11 right_digit, 12 response, 13 actual.answer, 16 RT, 17 strategy, 18 standard.set, 19 acc.recoded

### manually scale psize
# deleted multiple lines here


### manually scale psizex
DT3$c_psizex <- (DT3$psizex - 42)/20

# only need to run this once
# psizex <- sort(unique(DT3$psizex))
# c_psizex <- sort(unique(DT3$c_psizex))
# 
# psizex <- as.data.frame(cbind(psizex, c_psizex), col.names=c('psizex', 'c_psizex'))
# write.csv(psizex,"~/Documents/R/Pseudohomophone/psizex.csv")

### accuracy, etc.
DT3$correct <- factor(DT3$correct, levels = c(0,1,2,3)) # 1 = correct
# labels=c('Error', 'CorrectValid', 'InvalidCorrect', 'InvalidIncorrect'))
DT3$format <- factor(DT3$format, levels = c(1,2,3), labels = c("DIGIT","WORD", "PSEUDO"))

### format, dummy coding, with different reference group
# DT3$fmtD <- C(DT3$format, treatment, base = 1)
DT3$fmtW <- C(DT3$format, treatment, base = 2)
DT3$fmtP <- C(DT3$format, treatment, base = 3)

### deleted multiple lines here, old coding of fmt ...

# for running anova later
DT3$problemsize <- factor(DT3$problemsize, levels = c(1,2), labels = c('S', 'L')) #dichotomize psize

# have type coded both ways so I don't have to worry later
# nontie as the ref group
DT3$typexy <- C(factor(DT3$tie, levels = c(0,1), labels=c('NTIE','TIE')), treatment, base = 1)
# tie as the ref group
DT3$typexx <- C(factor(DT3$tie, levels = c(0,1), labels=c('NTIE','TIE')), treatment, base = 2)


DT3$strategy[DT3$strategy>5]<-NA
DT3$strategy <- factor(DT3$strategy, levels = c(1,2,3,4,5), 
                       labels=c('remember','count','transform','other','unsure'))
# retrieval as the ref group
DT3$retr <- C(factor(with(DT3, ifelse(DT3$strategy == 'remember', 0, 1)), levels = c(0,1), labels = c("RETR","NRETR")), treatment, base=1)
DT3$nretr <- C(factor(with(DT3, ifelse(DT3$strategy == 'remember', 0, 1)), levels = c(0,1), labels = c("RETR","NRETR")), treatment, base = 2)  # non-retrieval as the ref group
DT3$operation <- factor(DT3$operation, levels = c(1,2), labels = c("ADD","MUL"))
contrasts(DT3$operation) <- contr.treatment(2)
# base = 1

# add random effects of item, need to create item index first. 
# name items, 2x3 is different from 3x2, two x three is also diff from 2x3
DT3$item <- NULL
for (i in 1:length(DT3$sub)) {
  if (DT3$operation[i] == 'ADD' & DT3$format[i] == 'DIGIT') {
    DT3$item[i] <- paste(DT3$Lvalue[i], DT3$Rvalue[i], 'a','d', sep='')
  } else if (DT3$operation[i] == 'ADD' & DT3$format[i] == 'DIGIT') {
    DT3$item[i] <- paste(DT3$Lvalue[i], DT3$Rvalue[i], 'a','w', sep='')
  } else if (DT3$operation[i] == 'ADD' & DT3$format[i] == 'PSEUDO') {
    DT3$item[i] <- paste(DT3$Lvalue[i], DT3$Rvalue[i], 'a','p', sep='')
  } else if (DT3$operation[i] == 'MUL' & DT3$format[i] == 'DIGIT') {
    DT3$item[i] <- paste(DT3$Lvalue[i], DT3$Rvalue[i], 'm','d', sep='')
  } else if (DT3$operation[i] == 'MUL' & DT3$format[i] == 'WORD') {
    DT3$item[i] <- paste(DT3$Lvalue[i], DT3$Rvalue[i], 'm','w', sep='')
  } else if (DT3$operation[i] == 'MUL' & DT3$format[i] == 'PSEUDO') {
    DT3$item[i] <- paste(DT3$Lvalue[i], DT3$Rvalue[i], 'm','p', sep='')
  }
}

# convert sub and item to factor if they are not already factors
if (!is.factor(DT3$sub)) {
  DT3$sub <- factor(DT3$sub)
}
if (!is.factor(DT3$item)) {
  DT3$item <- factor(DT3$item)
}

# missing data
DT3$missingness <- factor(with(DT3, ifelse(DT3$correct == 1, 0, 1)), levels=c(0, 1))
# only need to run this once
# missing_count <- as.data.frame(count(subset(DT3, standard.set==1), vars=c('missingness', 'operation', 'psizex', 'typexy')))
# ggplot(missing_count, aes(x = psizex, y = freq)) + geom_line(aes(color = missingness)) + facet_grid(typexy ~ operation)
# accuracy decreases on large mul nontie, as high as about 50%
# also looked at missingness across formats, did not spot obvious difference, so collapsed over fmt

# add alpha to avoid overplotting
DT3 = ddply(DT3, .(sub), function(x){x$alpha = ifelse(runif(n = 1) > 0.8, 0.9, 0.1)
x})

# trimming RT >= 250
DT3_RT250 <- subset(DT3, RT >= 250 & correct == 1 & standard.set == 1)
DT3_RT250$RTlog <- log(DT3_RT250$RT)

### deleted multiple lines here

# Diagnostic -------------------------
### scatter plot
ggplot(DT3_RT250, aes(x = psizex, y = RT)) + ylim(0, 4000) + aes(group=factor(sub)) + facet_grid(format~typexy+operation) + geom_line(stat="smooth",method = "loess", colour = 'grey45', aes(alpha=alpha)) + stat_summary(aes(group = 1), geom = "line", fun.y = mean, linetype = "dotted", colour='red') + guides(alpha=F,colour=F)
# some outliers, use ylim to limit.  Nontie problems look bumpy, no obvious turning point, but does not look strictly linear either

### Diagnositc for m410
# maximal str, all trials, nontie coded as the ref group
m410DT3_ML <- lmer(RTlog~c_psizex*typexy*fmtW*operation + (1+c_psizex+typexy+fmtW|sub) + (1|item), REML=F, data=DT3_RT250)

### normality of outcome
hist(DT3_RT250$RT)
# negatively skewed
hist(DT3_RT250$RTlog)
# improved, still negatively skewed
# Level 1 diagnostic
# Homo
# residual ~ sub
plot(m410DT3_ML,sub~resid(., scaled=TRUE))
# outliers on both sides, more residuals on the positive side

tmp <- resid(m410DT3_ML)
DT3_RT250$resid<-tmp

tmp_a <- sort(tmp) # ascending
tmp_d <- sort(tmp, decreasing= T)# descending

# delete up to 50 cases on both head and tail and do sensitivity analysis
# DT3_25x2_RT250 <- subset(DT3_RT250, resid > tmp_a[25] & resid < tmp_d[25])
# delete up to 50 cases on tail and do sensitivity analysis
DT3_50_RT250 <- subset(DT3_RT250, resid < tmp_d[50])

### sensitivity analysis begin
# m410DT3_25x2 <- lmer(RTlog~c_psizex*typexy*fmtW*operation + (1+c_psizex+typexy+fmtW|sub) + (1|item), REML=F, data=DT3_25x2_RT250)
m410DT3_50 <- lmer(RTlog~c_psizex*typexy*fmtW*operation + (1+c_psizex+typexy+fmtW|sub) + (1|item), REML=F, data=DT3_50_RT250)

# plot(m410DT3_25x2,sub~resid(., scaled=TRUE))
# when I deleted 25 residuals on both head and tail, the plot looks a bit more balanced, on the right side, still some residuals as high as 4SD away, so decided to deleted the biggest 50 residuals and plot again
plot(m410DT3_50,sub~resid(., scaled=TRUE))
# looks even more balanced

# compare results
# random effects
randsens_m410DT3_50 <- cbind(as.data.frame(VarCorr(m410DT3_ML))[,c(1:3,5)], as.data.frame(VarCorr(m410DT3_50))[,5])
colnames(randsens_m410DT3_50) <- c('grp','var1', 'var2', 'Full','50OutliersDel')

#fixed effect
summ_m410DT3 <- summary(m410DT3_ML)
summ_m410DT3_50 <- summary(m410DT3_50)
fixedsens_m410DT3_50 <-cbind(as.data.frame(summ_m410DT3$coefficients)[,c(1,3)], as.data.frame(summ_m410DT3_50$coefficients)[,c(1,3)])
colnames(fixedsens_m410DT3_50) <- c('full.est','full.t','50OutliersDel.est','50OutliersDel.t')

write.csv(randsens_m410DT3_50, file='~/Documents/R/Pseudohomophone/randsens_m410DT3_50.csv')
write.csv(fixedsens_m410DT3_50, file='~/Documents/R/Pseudohomophone/fixedsens_m410DT3_50.csv')
# size x op changed from sig. to boarderline in reduced sample, otherwise, everything looks similar enough
########! sensitivity analysis end

# residual ~ pred
plot(m410DT3_ML)
# no discernable relationship between residual and fitted value, suggest appropriate model specification

# residual ~ psize
ggplot(m410DT3_ML, aes(DT3_RT250$psizex,residuals(m410DT3_ML,scaled=T))) + geom_point() + geom_smooth(method = "loess", size=1.5)
# looks linear

plot(m410DT3_ML, operation ~ resid(., scaled=TRUE))

plot(m410DT3_ML, typexy ~ resid(., scaled=TRUE))

plot(m410DT3_ML, fmtW ~ resid(., scaled=TRUE))
# visible residuals in all of them

#residual & predicted
xyplot(resid(m410DT3_ML)~fitted(m410DT3_ML))
# No discernable outliers

#normality of residuals
hist(residuals(m410DT3_ML,scaled=T))
qqnorm(residuals(m410DT3_ML,scaled=T))
qqnorm(residuals(m410DT3_50,scaled=T))
# not perfectly normal, but close enough, in the reduced sample, qq plot looks slightly improved

# Level 2
# extract random effects. ranef is structured as a list; at [1], it is the random int of item; at [2], random int & slps of sub.  See 2015-11-09 12:20 PM in OneNote
Level2_1_m410DT3 <- as.data.frame(ranef(m410DT3_ML)[1])
Level2_2_m410DT3 <- as.data.frame(ranef(m410DT3_ML)[2])
names(Level2_1_m410DT3) <- c("uIntItem")
# random int of items are grouped on items
names(Level2_2_m410DT3) <- c('uIntSub','uSize','utypexy','ufmtDW','ufmtPW')
#

# linearity & normality
hist(Level2_1_m410DT3$uIntItem)
# looks normal
ggpairs(data=Level2_2_m410DT3, title = "Relation between random effects")
# signs of correlation between int, psize, tie, and two dummy for format, no clear sign of outliers
# In DT1, I also merged level-2 random effects back into the dataset and plotted them against level-2 predictor
# Merge level 2 back to data
DT3_L1L2_RT250 <- merge(DT3_RT250, Level2_2_m410DT3, by.x="sub", by.y="row.names")
#only against level-2 predictors
boxplot(uIntSub~operation, data=DT3_L1L2_RT250)
# one potential outlier on mul, lower end
boxplot(uSize~operation, data=DT3_L1L2_RT250)
# one boarderline outlier on mul, higher end
boxplot(utypexy~operation, data=DT3_L1L2_RT250)
# no visible outlier
boxplot(ufmtDW~operation, data=DT3_L1L2_RT250)
# two potential outlier on mul, lower end
boxplot(ufmtPW~operation, data=DT3_L1L2_RT250)
# no visibal outlier

# A couple potential outliers, flag maxima/minima on random effects
uIntSub.minima <- unique(DT3_L1L2_RT250[order(DT3_L1L2_RT250$uIntSub),]$sub)[1:1]
# [1] 3
uSize.maxima <- unique(DT3_L1L2_RT250[order(DT3_L1L2_RT250$uSize,decreasing=T),]$sub)[1:1]
# [1] 47
ufmtDW.minima <- unique(DT3_L1L2_RT250[order(DT3_L1L2_RT250$ufmtDW),]$sub)[1:2]
# [1] 9 31

#sensitivity analysis
m410DT3_n3 <-lmer(RTlog~c_psizex*typexy*fmtW*operation + (1+c_psizex+typexy+fmtW|sub) + (1|item), REML=F, data=subset(DT3_RT250,sub!='3'))
m410DT3_n9 <-lmer(RTlog~c_psizex*typexy*fmtW*operation + (1+c_psizex+typexy+fmtW|sub) + (1|item), REML=F, data=subset(DT3_RT250,sub!='9'))
m410DT3_n31 <-lmer(RTlog~c_psizex*typexy*fmtW*operation + (1+c_psizex+typexy+fmtW|sub) + (1|item), REML=F, data=subset(DT3_RT250,sub!='31'))
m410DT3_n47 <-lmer(RTlog~c_psizex*typexy*fmtW*operation + (1+c_psizex+typexy+fmtW|sub) + (1|item), REML=F, data=subset(DT3_RT250,sub!='47'))

# compare results
# random effects omitted
#fixed effects
summ_m410DT3_n3 <- summary(m410DT3_n3)
summ_m410DT3_n9 <- summary(m410DT3_n9)
summ_m410DT3_n31 <- summary(m410DT3_n31)
summ_m410DT3_n47 <- summary(m410DT3_n47)
fixed_sensitivity_L2<-cbind(
  as.data.frame(summ_m410DT3$coefficients)[,c(1,3)],
  as.data.frame(summ_m410DT3_n3$coefficients)[,c(1,3)],
  as.data.frame(summ_m410DT3_n9$coefficients)[,c(1,3)],
  as.data.frame(summ_m410DT3_n31$coefficients)[,c(1,3)],
  as.data.frame(summ_m410DT3_n47$coefficients)[,c(1,3)])
colnames(fixed_sensitivity_L2) <- c('full.est','full.t','full_3.est','full_3.t', 'full_9.est','full_9.t', 'full_31.est','full_31.t', 'full_47.est','full_47.t')

write.csv(fixed_sensitivity_L2, file='~/Documents/R/Pseudohomophone/fixed_sensitivity_L2.csv')
# n47 seems to have impact, size x operation changes from boarderline to insig when this case is deleted
###! diagnostic for m410

### take a closer look at these outliers identifed so far
ggplot(subset(DT3_L1L2_RT250, resid>=tmp_d[50]), aes(x = psizex, y = RT)) + facet_grid(format~typexy+operation) + geom_point()
# some RT look incredibly long, especially on word nontie multiplication
ggplot(subset(DT3_L1L2_RT250, sub == '47'), aes(x = psizex, y = RT)) + facet_grid(format~typexy+operation) + geom_line()
# this participant did multiplication problems, seems to have extreme RT data, consider deleting. Once deleted, the experiment design would be balanced between two operations (38 each)
# count(count(subset(DT3, sub!='47'), vars=c('sub', 'operation')), 'operation')
# operation freq
# 1       add 9180
# 2       mul 9180
###! diagnostic

### Anaylsis proper ------------
# trimming RT >= 250, deleting sub47, see diagnostic
DT3_trim <- subset(DT3, RT >= 250 & correct == 1 & standard.set == 1 & sub != '47')
DT3_trim$RTlog <- log(DT3_trim$RT)

# models in stage one (replication) ----
### m21*, only nontie, digit v word, separate operation
### addition, raw RT ----
m210dwDT3atr_REML <- lmer(RT ~ problemsize*fmtW + (1|sub), data=subset(DT3_trim, typexy=='NTIE' & format!= "PSEUDO" & operation == "ADD"))

m211dwDT3atr_REML <- lmer(RT ~ problemsize*fmtW + (1 + problemsize + fmtW|sub) + (1|item), data=subset(DT3_trim, typexy=='NTIE' & format!= "PSEUDO" & operation == "ADD"))

# compared w/o random intercept on item, random intercept on item is necessary.  See untidy version for code

m220dwDT3atr_REML <- lmer(RT~c_psizex*fmtW + (1|sub), data=subset(DT3_trim, typexy=='NTIE' & format!= "PSEUDO" & operation == "ADD"))

m221dwDT3atr_REML <- lmer(RT~c_psizex*fmtW + (1 + c_psizex + fmtW|sub) + (1|item), data=subset(DT3_trim, typexy=='NTIE' & format!= "PSEUDO" & operation == "ADD"))

fixef_m2xxdwDT3atr <- cbind(
  as.data.frame(summary(m210dwDT3atr_REML)$coefficients)[,c(1,3)],
  as.data.frame(summary(m211dwDT3atr_REML)$coefficients)[,c(1,3)],
  as.data.frame(summary(m220dwDT3atr_REML)$coefficients)[,c(1,3)],
  as.data.frame(summary(m221dwDT3atr_REML)$coefficients)[,c(1,3)]
)
colnames(fixef_m2xxdwDT3atr) <- c('m210_est','m210_t','m211_est','m211_t', 'm220_est','m220_t', 'm221_est','m221_t')

write.csv(fixef_m2xxdwDT3atr, file='~/Documents/R/Pseudohomophone/fixef_m2xxdwDT3atr.csv')

### addition, log RT ------
m210dwDT3at_REML <- lmer(RTlog ~ problemsize*fmtW + (1|sub), data=subset(DT3_trim, typexy=='NTIE' & format!= "PSEUDO" & operation == "ADD"))

m211dwDT3at_REML <- lmer(RTlog ~ problemsize*fmtW + (1 + problemsize + fmtW|sub) + (1|item), data=subset(DT3_trim, typexy=='NTIE' & format!= "PSEUDO" & operation == "ADD"))

m220dwDT3at_REML <- lmer(RTlog ~ c_psizex*fmtW + (1|sub), data=subset(DT3_trim, typexy=='NTIE' & format!= "PSEUDO" & operation == "ADD"))

m221dwDT3at_REML <- lmer(RTlog ~ c_psizex*fmtW + (1 + c_psizex + fmtW|sub) + (1|item), data=subset(DT3_trim, typexy=='NTIE' & format!= "PSEUDO" & operation == "ADD"))

fixef_m2xxdwDT3at <- cbind(
  as.data.frame(summary(m210DT3at_REML)$coefficients)[,c(1,3)],
  as.data.frame(summary(m211DT3at_REML)$coefficients)[,c(1,3)],
  as.data.frame(summary(m220DT3at_REML)$coefficients)[,c(1,3)],
  as.data.frame(summary(m221DT3at_REML)$coefficients)[,c(1,3)]
)
colnames(fixef_m2xxdwDT3at) <- c('m210_est','m210_t','m211_est','m211_t', 'm220_est','m220_t', 'm221_est','m221_t')

write.csv(fixef_m2xxdwDT3at, file='~/Documents/R/Pseudohomophone/fixef_m2xxdwDT3at.csv')


### multiplication, raw RT
m210DT3mtr_REML <- lmer(RT ~ problemsize*fmtW + (1|sub), data=subset(DT3_trim, typexy=='NTIE' & format!= "PSEUDO" & operation == "MUL"))

m211DT3mtr_REML <- lmer(RT ~ problemsize*fmtW + (1 + problemsize + fmtW|sub) + (1|item), data=subset(DT3_trim, typexy=='NTIE' & format!= "PSEUDO" & operation == "MUL"))

m220DT3mtr_REML <- lmer(RT~c_psizex*fmtW + (1|sub), data=subset(DT3_trim, typexy=='NTIE' & format!= "PSEUDO" & operation == "MUL"))

m221DT3mtr_REML <- lmer(RT~c_psizex*fmtW + (1 + c_psizex + fmtW|sub) + (1|item), data=subset(DT3_trim, typexy=='NTIE' & format!= "PSEUDO" & operation == "MUL"))

fixef_m2xxDT3mtr <- cbind(
  as.data.frame(summary(m210DT3mtr_REML)$coefficients)[,c(1,3)],
  as.data.frame(summary(m211DT3mtr_REML)$coefficients)[,c(1,3)],
  as.data.frame(summary(m220DT3mtr_REML)$coefficients)[,c(1,3)],
  as.data.frame(summary(m221DT3mtr_REML)$coefficients)[,c(1,3)]
)
colnames(fixef_m2xxDT3mtr) <- c('m210_est','m210_t','m211_est','m211_t', 'm220_est','m220_t', 'm221_est','m221_t')

write.csv(fixef_m2xxDT3mtr, file='~/Documents/R/Pseudohomophone/fixef_m2xxDT3mtr.csv')

### multiplication, log RT -----
m210DT3mt_REML <- lmer(RTlog ~ problemsize*fmtW + (1|sub), data=subset(DT3_trim, typexy=='NTIE' & format!= "PSEUDO" & operation == "MUL"))

m211DT3mt_REML <- lmer(RTlog ~ problemsize*fmtW + (1 + problemsize + fmtW|sub) + (1|item), data=subset(DT3_trim, typexy=='NTIE' & format!= "PSEUDO" & operation == "MUL"))

m220DT3mt_REML <- lmer(RTlog ~ c_psizex*fmtW + (1|sub), data=subset(DT3_trim, typexy=='NTIE' & format!= "PSEUDO" & operation == "MUL"))

m221DT3mt_REML <- lmer(RTlog ~ c_psizex*fmtW + (1 + c_psizex + fmtW|sub) + (1|item), data=subset(DT3_trim, typexy=='NTIE' & format!= "PSEUDO" & operation == "MUL"))

fixef_m2xxDT3mt <- cbind(
  as.data.frame(summary(m210DT3mt_REML)$coefficients)[,c(1,3)],
  as.data.frame(summary(m211DT3mt_REML)$coefficients)[,c(1,3)],
  as.data.frame(summary(m220DT3mt_REML)$coefficients)[,c(1,3)],
  as.data.frame(summary(m221DT3mt_REML)$coefficients)[,c(1,3)]
)
colnames(fixef_m2xxDT3mt) <- c('m210_est','m210_t','m211_est','m211_t', 'm220_est','m220_t', 'm221_est','m221_t')

write.csv(fixef_m2xxDT3mt, file='~/Documents/R/Pseudohomophone/fixef_m2xxDT3mt.csv')

### add retrieval to the picture ----
# log RT -----
m321dwDT3at_REML <- lmer(RTlog ~ c_psizex*fmtW*retr + (1 + c_psizex + fmtW + retr|sub) + (1|item), data=subset(DT3_trim, typexy=='NTIE' & format!= "PSEUDO" & operation == "ADD"))

m321dwDT3atr_REML <- lmer(RT ~ c_psizex*fmtW*retr + (1 + c_psizex + fmtW + retr|sub) + (1|item), data=subset(DT3_trim, typexy=='NTIE' & format!= "PSEUDO" & operation == "ADD"))

m321dwDT3mt_REML <- lmer(RTlog ~ c_psizex*fmtW*retr + (1 + c_psizex + fmtW + retr|sub) + (1|item), data=subset(DT3_trim, typexy=='NTIE' & format!= "PSEUDO" & operation == "MUL"))

m321dwDT3mtr_REML <- lmer(RT ~ c_psizex*fmtW*retr + (1 + c_psizex + fmtW + retr|sub) + (1|item), data=subset(DT3_trim, typexy=='NTIE' & format!= "PSEUDO" & operation == "MUL"))

# try random intx, but model failed to converge

fixef_m321dwDT3t <- cbind(
  as.data.frame(summary(m321dwDT3at_REML)$coefficients)[,c(1,3)],
  as.data.frame(summary(m321dwDT3atr_REML)$coefficients)[,c(1,3)],
  as.data.frame(summary(m321dwDT3mt_REML)$coefficients)[,c(1,3)],
  as.data.frame(summary(m321dwDT3mtr_REML)$coefficients)[,c(1,3)]
)
colnames(fixef_m321dwDT3t) <- c('m321dwDT3at_est','m321dwDT3at_t','m321dwDT3atr_est','m321dwDT3atr_t', 'm321dwDT3mt_est','m321dwDT3mt_t', 'm321dwDT3mtr_est','m321dwDT3mtr_t')

write.csv(fixef_m321dwDT3t, file='~/Documents/R/Pseudohomophone/fixef_m321dwDT3t.csv')

# models in stage two (Digit vs Pseudo) -----
m321dpDT3at_REML <- lmer(RTlog ~ c_psizex*fmtP*retr + (1 + c_psizex + fmtP + retr|sub) + (1|item), data=subset(DT3_trim, typexy=='NTIE' & format!= "WORD" & operation == "ADD"))

m321dpDT3atr_REML <- lmer(RT ~ c_psizex*fmtP*retr + (1 + c_psizex + fmtP + retr|sub) + (1|item), data=subset(DT3_trim, typexy=='NTIE' & format!= "WORD" & operation == "ADD"))

m321dpDT3mt_REML <- lmer(RTlog ~ c_psizex*fmtP*retr + (1 + c_psizex + fmtP + retr|sub) + (1|item), data=subset(DT3_trim, typexy=='NTIE' & format!= "WORD" & operation == "MUL"))

m321dpDT3mtr_REML <- lmer(RT ~ c_psizex*fmtP*retr + (1 + c_psizex + fmtP + retr|sub) + (1|item), data=subset(DT3_trim, typexy=='NTIE' & format!= "WORD" & operation == "MUL"))


fixef_m321dpDT3t <- cbind(
  as.data.frame(summary(m321dpDT3at_REML)$coefficients)[,c(1,3)],
  as.data.frame(summary(m321dpDT3atr_REML)$coefficients)[,c(1,3)],
  as.data.frame(summary(m321dpDT3mt_REML)$coefficients)[,c(1,3)],
  as.data.frame(summary(m321dpDT3mtr_REML)$coefficients)[,c(1,3)]
)
colnames(fixef_m321dpDT3t) <- c('m321dpDT3at_est','m321dpDT3at_t','m321dpDT3atr_est','m321dpDT3atr_t', 'm321dpDT3mt_est','m321dpDT3mt_t', 'm321dpDT3mtr_est','m321dpDT3mtr_t')

write.csv(fixef_m321dpDT3t, file='~/Documents/R/Pseudohomophone/fixef_m321dpDT3t.csv')

### m2xx
m221dpDT3at_REML <- lmer(RTlog ~ c_psizex*fmtP + (1 + c_psizex + fmtP|sub) + (1|item), data=subset(DT3_trim, typexy=='NTIE' & format!= "WORD" & operation == "ADD"))

m221dpDT3atr_REML <- lmer(RT ~ c_psizex*fmtP + (1 + c_psizex + fmtP|sub) + (1|item), data=subset(DT3_trim, typexy=='NTIE' & format!= "WORD" & operation == "ADD"))

m221dpDT3mt_REML <- lmer(RTlog ~ c_psizex*fmtP + (1 + c_psizex + fmtP|sub) + (1|item), data=subset(DT3_trim, typexy=='NTIE' & format!= "WORD" & operation == "MUL"))

m221dpDT3mtr_REML <- lmer(RT ~ c_psizex*fmtP + (1 + c_psizex + fmtP|sub) + (1|item), data=subset(DT3_trim, typexy=='NTIE' & format!= "WORD" & operation == "MUL"))


fixef_m221dpDT3t <- cbind(
  as.data.frame(summary(m221dpDT3atr_REML)$coefficients)[,c(1,3)],
  as.data.frame(summary(m221dpDT3at_REML)$coefficients)[,c(1,3)],
  as.data.frame(summary(m221dpDT3mtr_REML)$coefficients)[,c(1,3)],
  as.data.frame(summary(m221dpDT3mt_REML)$coefficients)[,c(1,3)]
)
colnames(fixef_m221dpDT3t) <- c('m221dpDT3atr_est', 'm221dpDT3atr_t', 'm221dpDT3at_est', 'm221dpDT3at_t', 'm221dpDT3mtr_est', 'm221dpDT3mtr_t', 'm221dpDT3mt_est', 'm221dpDT3mt_t')

write.csv(fixef_m221dpDT3t, file='~/Documents/R/Pseudohomophone/fixef_m221dpDT3t.csv')



### prepare data for spss ------
DT3w_trim <- dcast(DT3_trim, sub + operation ~ format + typexy + problemsize, mean, value.var = "RT")
DT3w_trimlog <- dcast(DT3_trim, sub + operation ~ format + typexy + problemsize, mean, value.var = "RTlog")

# export DT3w and run anova on them using SPSS
write.csv(DT3w_trim, "~/Desktop/DT3w_trim.csv")
write.csv(DT3w_trimlog, "~/Desktop/DT3w_trimlog.csv")

## I could also aggregate data for 4-way ANOVA, with retr in it, but there are empty cells, cannot run ANOVA
# DT3w4_trim <- dcast(subset(DT3_trim, typexy = 'NTIE'), sub + operation ~ problemsize + format + retr, mean, value.var = "RT")
# DT3w4_trimlog <- dcast(subset(DT3_trim, typexy = 'NTIE'), sub + operation ~ problemsize + format + retr, mean, value.var = "RTlog")
# 
# # export DT3w4
# write.csv(DT3w4_trim, "~/Desktop/DT3w4_trim.csv")
# write.csv(DT3w4_trimlog, "~/Desktop/DT3w4_trimlog.csv")


### compute point estimates ------
# DV == logRT ---
fixef_m221dwDT3at <- as.data.frame(summary(m221dwDT3at_REML)$coefficients)
fixef_m221dpDT3at <- as.data.frame(summary(m221dpDT3at_REML)$coefficients)
fixef_m221dwDT3mt <- as.data.frame(summary(m221dwDT3mt_REML)$coefficients)
fixef_m221dpDT3mt <- as.data.frame(summary(m221dpDT3mt_REML)$coefficients)


ADD_NTIE_DIGIT <- fixef_m221dwDT3at['(Intercept)', 1] + fixef_m221dwDT3at["c_psizex", 1]*c_psizex

ADD_NTIE_WORD <- fixef_m221dwDT3at['(Intercept)', 1] + fixef_m221dwDT3at["c_psizex", 1]*c_psizex + fixef_m221dwDT3at["fmtWWORD", 1]*1 + fixef_m221dwDT3at["c_psizex:fmtWWORD",1]*c_psizex

ADD_NTIE_PSEUDO <- fixef_m221dpDT3at['(Intercept)', 1] + fixef_m221dpDT3at["c_psizex", 1]*c_psizex + fixef_m221dpDT3at["fmtPPSEUDO", 1]*1 + fixef_m221dpDT3at["c_psizex:fmtPPSEUDO",1]*c_psizex

MUL_NTIE_DIGIT <- fixef_m221dwDT3mt['(Intercept)', 1] + fixef_m221dwDT3mt["c_psizex", 1]*c_psizex

MUL_NTIE_WORD <- fixef_m221dwDT3mt['(Intercept)', 1] + fixef_m221dwDT3mt["c_psizex", 1]*c_psizex + fixef_m221dwDT3mt["fmtWWORD", 1]*1 + fixef_m221dwDT3mt["c_psizex:fmtWWORD",1]*c_psizex

MUL_NTIE_PSEUDO <- fixef_m221dpDT3mt['(Intercept)', 1] + fixef_m221dpDT3mt["c_psizex", 1]*c_psizex + fixef_m221dpDT3mt["fmtPPSEUDO", 1]*1 + fixef_m221dpDT3mt["c_psizex:fmtPPSEUDO",1]*c_psizex

pointest_m221DT3t <- data.frame(psizex,c_psizex,ADD_NTIE_DIGIT, ADD_NTIE_WORD, ADD_NTIE_PSEUDO, MUL_NTIE_DIGIT, MUL_NTIE_WORD, MUL_NTIE_PSEUDO)


pointestlog_m221DT3t <- pointest_m221DT3t
pointestlog_m221DT3t[3:8] <- lapply(pointest_m221DT3t[3:8], function(x) exp(x))

pointestlogL_m221DT3t <- melt(pointestlog_m221DT3t, id.vars=c('psizex','c_psizex'),variable.name='pfeature',value.name='RT_fixpred')

pointestlogL_m221DT3t <- separate(pointestlogL_m221DT3t, pfeature, into = c("operation", 'typexy', "format"), sep = '_')
pointestlogL_m221DT3t$operation <- factor(pointestlogL_m221DT3t$operation)

pointestlogL_m221DT3t$format <- factor(pointestlogL_m221DT3t$format)


# I looked at point estimate for digit based on m221dw and m221dp, they are very close (+/- 1 ms)

## point estimate with retrieval included
fixef_m321dwDT3at <- as.data.frame(summary(m321dwDT3at_REML)$coefficients)
fixef_m321dpDT3at <- as.data.frame(summary(m321dpDT3at_REML)$coefficients)
fixef_m321dwDT3mt <- as.data.frame(summary(m321dwDT3mt_REML)$coefficients)
fixef_m321dpDT3mt <- as.data.frame(summary(m321dpDT3mt_REML)$coefficients)

RETR_DIGIT_ADD_NTIE <- fixef_m321dwDT3at['(Intercept)', 1] + fixef_m321dwDT3at["c_psizex", 1]*c_psizex 

RETR_WORD_ADD_NTIE <- fixef_m321dwDT3at['(Intercept)', 1] + fixef_m321dwDT3at["c_psizex", 1]*c_psizex + fixef_m321dwDT3at["fmtWWORD", 1]*1 + fixef_m321dwDT3at["c_psizex:fmtWWORD",1]*c_psizex

RETR_PSEUDO_ADD_NTIE <- fixef_m321dpDT3at['(Intercept)', 1] + fixef_m321dpDT3at["c_psizex", 1]*c_psizex + fixef_m321dpDT3at["fmtPPSEUDO", 1]*1 + fixef_m321dpDT3at["c_psizex:fmtPPSEUDO",1]*c_psizex

NRETR_DIGIT_ADD_NTIE <- fixef_m321dwDT3at['(Intercept)', 1] + fixef_m321dwDT3at["c_psizex", 1]*c_psizex  + fixef_m321dwDT3at["retr2",1]*1 +  fixef_m321dwDT3at["c_psizex:retr2",1]*c_psizex 

NRETR_WORD_ADD_NTIE <- fixef_m321dwDT3at['(Intercept)', 1] + fixef_m321dwDT3at["c_psizex", 1]*c_psizex + 
  fixef_m321dwDT3at["fmtWWORD", 1]*1 + fixef_m321dwDT3at["c_psizex:fmtWWORD",1]*c_psizex + fixef_m321dwDT3at["retr2",1]*1 + 
  fixef_m321dwDT3at["c_psizex:retr2",1]*c_psizex + fixef_m321dwDT3at["fmtWWORD:retr2",1]*1 + fixef_m321dwDT3at["c_psizex:fmtWWORD:retr2",1]*c_psizex 

NRETR_PSEUDO_ADD_NTIE <- fixef_m321dpDT3at['(Intercept)', 1] + fixef_m321dpDT3at["c_psizex", 1]*c_psizex + 
  fixef_m321dpDT3at["fmtPPSEUDO", 1]*1 + 
  fixef_m321dpDT3at["retr2",1]*1 + 
  fixef_m321dpDT3at["c_psizex:fmtPPSEUDO",1]*c_psizex + fixef_m321dpDT3at["c_psizex:retr2",1]*c_psizex + fixef_m321dpDT3at["fmtPPSEUDO:retr2",1]*1 + fixef_m321dpDT3at["c_psizex:fmtPPSEUDO:retr2",1]*c_psizex


RETR_DIGIT_MUL_NTIE <- fixef_m321dwDT3mt['(Intercept)', 1] + fixef_m321dwDT3mt["c_psizex", 1]*c_psizex 

RETR_WORD_MUL_NTIE <- fixef_m321dwDT3mt['(Intercept)', 1] + fixef_m321dwDT3mt["c_psizex", 1]*c_psizex + fixef_m321dwDT3mt["fmtWWORD", 1]*1 + fixef_m321dwDT3mt["c_psizex:fmtWWORD",1]*c_psizex

RETR_PSEUDO_MUL_NTIE <- fixef_m321dpDT3mt['(Intercept)', 1] + fixef_m321dpDT3mt["c_psizex", 1]*c_psizex + fixef_m321dpDT3mt["fmtPPSEUDO", 1]*1 + fixef_m321dpDT3mt["c_psizex:fmtPPSEUDO",1]*c_psizex

NRETR_DIGIT_MUL_NTIE <- fixef_m321dwDT3mt['(Intercept)', 1] + fixef_m321dwDT3mt["c_psizex", 1]*c_psizex  + fixef_m321dwDT3mt["retr2",1]*1 +  fixef_m321dwDT3mt["c_psizex:retr2",1]*c_psizex 

NRETR_WORD_MUL_NTIE <- fixef_m321dwDT3mt['(Intercept)', 1] + fixef_m321dwDT3mt["c_psizex", 1]*c_psizex + 
  fixef_m321dwDT3mt["fmtWWORD", 1]*1 + fixef_m321dwDT3mt["c_psizex:fmtWWORD",1]*c_psizex + fixef_m321dwDT3mt["retr2",1]*1 + 
  fixef_m321dwDT3mt["c_psizex:retr2",1]*c_psizex + fixef_m321dwDT3mt["fmtWWORD:retr2",1]*1 + fixef_m321dwDT3mt["c_psizex:fmtWWORD:retr2",1]*c_psizex 

NRETR_PSEUDO_MUL_NTIE <- fixef_m321dpDT3mt['(Intercept)', 1] + fixef_m321dpDT3mt["c_psizex", 1]*c_psizex + 
  fixef_m321dpDT3mt["fmtPPSEUDO", 1]*1 + 
  fixef_m321dpDT3mt["retr2",1]*1 + 
  fixef_m321dpDT3mt["c_psizex:fmtPPSEUDO",1]*c_psizex + fixef_m321dpDT3mt["c_psizex:retr2",1]*c_psizex + fixef_m321dpDT3mt["fmtPPSEUDO:retr2",1]*1 + fixef_m321dpDT3mt["c_psizex:fmtPPSEUDO:retr2",1]*c_psizex


pointest_m321DT3t <- data.frame(psizex,c_psizex,RETR_DIGIT_ADD_NTIE, RETR_WORD_ADD_NTIE, RETR_PSEUDO_ADD_NTIE, RETR_DIGIT_MUL_NTIE, RETR_WORD_MUL_NTIE, RETR_PSEUDO_MUL_NTIE, NRETR_DIGIT_ADD_NTIE, NRETR_WORD_ADD_NTIE, NRETR_PSEUDO_ADD_NTIE, NRETR_DIGIT_MUL_NTIE, NRETR_WORD_MUL_NTIE, NRETR_PSEUDO_MUL_NTIE)

pointestlog_m321DT3t <- pointest_m321DT3t
pointestlog_m321DT3t[3:14] <- lapply(pointest_m321DT3t[3:14], function(x) exp(x))

pointestlogL_m321DT3t <- melt(pointestlog_m321DT3t, id.vars=c('psizex','c_psizex'),variable.name='pfeature',value.name='RT_fixpred')

pointestlogL_m321DT3t <- separate(pointestlogL_m321DT3t, pfeature, into = c("retr", "format", "operation", "typexy"), sep = '_')

# 2017-05-26   If I don't specify the lables, the figure would put NRETR above RETR, as in alphabetical
pointestlogL_m321DT3t$retr <- factor(with(pointestlogL_m321DT3t, ifelse(pointestlogL_m321DT3t$retr == 'RETR', 0 ,1)), levels = c(0 ,1), labels = c("RETR","NRETR"))

pointestlogL_m321DT3t$operation <- factor(pointestlogL_m321DT3t$operation)

pointestlogL_m321DT3t$format <- factor(pointestlogL_m321DT3t$format)


### ggplot -----
# In phase2, I merged point estimate with DT3_trim before plotting, but in this version, I no longer need to merge.
# DV == log RT -----
# two-way
pSizeFormat_m221dwDT3at <- ggplot(subset(pointestlogL_m221DT3t, format != 'PSEUDO' & operation == 'ADD'), aes(x = psizex, y=RT_fixpred, linetype = format, color = format)) + geom_line() + scale_color_manual(values=c('black','black', 'black')) + labs(x='Problem Size', y='Mean Response Time (ms)') + ylim(800,4100) + apatheme

pSizeFormat_m221dpDT3at <- ggplot(subset(pointestlogL_m221DT3t, format != 'WORD' & operation == 'ADD'), aes(x = psizex, y=RT_fixpred, linetype = format, color = format)) + geom_line() + scale_color_manual(values=c('black','black', 'black')) + labs(x='Problem Size', y='Mean Response Time (ms)') + ylim(800,4100) + apatheme

pSizeFormat_m221dwDT3mt <- ggplot(subset(pointestlogL_m221DT3t, format != 'PSEUDO' & operation == 'MUL'), aes(x = psizex, y=RT_fixpred, linetype = format, color = format)) + geom_line() + scale_color_manual(values=c('black','black', 'black')) + labs(x='Problem Size', y='Mean Response Time (ms)') + ylim(800,4100) + apatheme

pSizeFormat_m221dpDT3mt <- ggplot(subset(pointestlogL_m221DT3t, format != 'WORD' & operation == 'MUL'), aes(x = psizex, y=RT_fixpred, linetype = format, color = format)) + geom_line() + scale_color_manual(values=c('black','black', 'black')) + labs(x='Problem Size', y='Mean Response Time (ms)') + ylim(800,4100) + apatheme

# three-way
pSizeFormat_m321dwDT3at <- ggplot(subset(pointestlogL_m321DT3t, format != 'PSEUDO' & operation == 'ADD'), aes(x = psizex, y=RT_fixpred, linetype = format, color = format)) + geom_line() + scale_color_manual(values=c('black','black')) + labs(x='Problem Size', y='Mean Response Time (ms)') + ylim(800,4100) + apatheme + facet_grid(retr~.)

pSizeFormat_m321dpDT3at <- ggplot(subset(pointestlogL_m321DT3t, format != 'WORD' & operation == 'ADD'), aes(x = psizex, y=RT_fixpred, linetype = format, color = format)) + geom_line() + scale_color_manual(values=c('black','black')) + labs(x='Problem Size', y='Mean Response Time (ms)') + ylim(800,4100) + apatheme + facet_grid(retr~.)

pSizeFormat_m321dwDT3mt <- ggplot(subset(pointestlogL_m321DT3t, format != 'PSEUDO' & operation == 'MUL'), aes(x = psizex, y=RT_fixpred, linetype = format, color = format)) + geom_line() + scale_color_manual(values=c('black','black')) + labs(x='Problem Size', y='Mean Response Time (ms)') + ylim(800,4100) + apatheme + facet_grid(retr~.)

pSizeFormat_m321dpDT3mt <- ggplot(subset(pointestlogL_m321DT3t, format != 'WORD' & operation == 'MUL'), aes(x = psizex, y=RT_fixpred, linetype = format, color = format)) + geom_line() + scale_color_manual(values=c('black','black')) + labs(x='Problem Size', y='Mean Response Time (ms)') + ylim(800,4100) + apatheme + facet_grid(retr~.)

# DV == raw RT -----
DT3_plotm3xx <- aggregate(RT ~ format + psizex + c_psizex + operation + retr, subset(DT3_trim, typexy == 'NTIE'), mean)

DT3_plotm2xx <- aggregate(RT ~ format + psizex + c_psizex + operation, subset(DT3_trim, typexy == 'NTIE'), mean)

pSizeFormat_m221dwDT3atr <- ggplot(subset(DT3_plotm2xx, format != "PSEUDO" & operation == "ADD"), aes(x = psizex, y=RT, linetype = format, color = format)) + geom_smooth(method = lm, se = F) + scale_color_manual(values=c('black','black')) + labs(x='Problem Size', y='Mean Response Time (ms)') + ylim(800,4100) + apatheme

pSizeFormat_m221dpDT3atr <- ggplot(subset(DT3_plotm2xx, format != "WORD" & operation == "ADD"), aes(x = psizex, y=RT, linetype = format, color = format)) + geom_smooth(method = lm, se = F) + scale_color_manual(values=c('black','black')) + labs(x='Problem Size', y='Mean Response Time (ms)') + ylim(800,4100) + apatheme

pSizeFormat_m221dwDT3mtr <- ggplot(subset(DT3_plotm2xx, format != "PSEUDO" & operation == "MUL"), aes(x = psizex, y=RT, linetype = format, color = format)) + geom_smooth(method = lm, se = F) + scale_color_manual(values=c('black','black')) + labs(x='Problem Size', y='Mean Response Time (ms)') + ylim(800,4100) + apatheme

pSizeFormat_m221dpDT3mtr <- ggplot(subset(DT3_plotm2xx, format != "WORD" & operation == "MUL"), aes(x = psizex, y=RT, linetype = format, color = format)) + geom_smooth(method = lm, se = F) + scale_color_manual(values=c('black','black')) + labs(x='Problem Size', y='Mean Response Time (ms)') + ylim(800,4100) + apatheme

pSizeFormat_m321dwDT3atr <- ggplot(subset(DT3_plotm3xx, format != "PSEUDO" & operation == "ADD"), aes(x = psizex, y=RT, linetype = format, color = format)) + geom_smooth(method = lm, se = F) + scale_color_manual(values=c('black','black')) + labs(x='Problem Size', y='Mean Response Time (ms)') + ylim(800,4100) + apatheme + facet_grid(retr~.)

pSizeFormat_m321dpDT3atr <- ggplot(subset(DT3_plotm3xx, format != "WORD" & operation == "ADD"), aes(x = psizex, y=RT, linetype = format, color = format)) + geom_smooth(method = lm, se = F) + scale_color_manual(values=c('black','black')) + labs(x='Problem Size', y='Mean Response Time (ms)') + ylim(800,4100) + apatheme + facet_grid(retr~.)

pSizeFormat_m321dwDT3mtr <- ggplot(subset(DT3_plotm3xx, format != "PSEUDO" & operation == "MUL"), aes(x = psizex, y=RT, linetype = format, color = format)) + geom_smooth(method = lm, se = F) + scale_color_manual(values=c('black','black')) + labs(x='Problem Size', y='Mean Response Time (ms)') + ylim(800,4100) + apatheme + facet_grid(retr~.)

pSizeFormat_m321dpDT3mtr <- ggplot(subset(DT3_plotm3xx, format != "WORD" & operation == "MUL"), aes(x = psizex, y=RT, linetype = format, color = format)) + geom_smooth(method = lm, se = F) + scale_color_manual(values=c('black','black')) + labs(x='Problem Size', y='Mean Response Time (ms)') + ylim(800,4100) + apatheme + facet_grid(retr~.)


grid.arrange(pSizeFormat_m221dwDT3atr, pSizeFormat_m221dpDT3atr, pSizeFormat_m221dwDT3at, pSizeFormat_m221dpDT3at, pSizeFormat_m321dwDT3atr, pSizeFormat_m321dpDT3atr, pSizeFormat_m321dwDT3at, pSizeFormat_m321dpDT3at, ncol=4, nrow =2)

grid.arrange(pSizeFormat_m221dwDT3mtr, pSizeFormat_m221dpDT3mtr, pSizeFormat_m221dwDT3mt, pSizeFormat_m221dpDT3mt, pSizeFormat_m321dwDT3mtr, pSizeFormat_m321dpDT3mtr, pSizeFormat_m321dwDT3mt, pSizeFormat_m321dpDT3mt, ncol=4, nrow =2)
# addjust lengend, y-axis [ ]





# ANOVA plots ----
# raw RT
DT3_anovaplot <- aggregate(RT~format + operation + problemsize, data=subset(DT3_trim, typexy == 'NTIE' & format != 'PSEUDO'), mean)
tcrit <- qt(.975, 33)
# by running two separate anova, this tcrit is also smaller. If opration is included (i.e., a 3-way ANOVA), this df would be twice as large. CI for the same effect in a 3-way ANOVA is: qt(.975, 66)*sqrt(26365/68) = 39.31359, this number is a bit larger than CIa, but a lot smaller than CIm
CIa <- tcrit*sqrt(4913/34) # CI for addition
CIm <- tcrit*sqrt(47818/34) # CI for multiplication

DT3_anovaplot$CI <- ifelse(DT3_anovaplot$operation == 'ADD', CIa, CIm)

pd<-position_dodge(0.1)

pDigitWord_DT3anova <- ggplot(DT3_anovaplot, aes(x=problemsize, y=RT, color=format, linetype=format, group = format)) + geom_point(position=pd, aes(shape=format)) + geom_errorbar(aes(ymin=RT-CI, ymax=RT+CI), colour="black", width=.1, linetype = 1, position=pd) + geom_line(position=pd) + ylim(900,2500) + scale_color_manual(values=c('black','black')) + labs(x='Problem Size', y='Mean Response Time (ms)') + facet_grid(.~operation) + apatheme
ggsave(file='~/Documents/R/Pseudohomophone/pDigitWord_DT3anova.png')


# log RT
DT3_anovaplot_log <- aggregate(RTlog~format + operation + problemsize, data=subset(DT3_trim, typexy == 'NTIE' & format != 'PSEUDO'), mean)
tcrit <- qt(.975, 33)

CIa_log <- tcrit*sqrt(.002/34) # CI for addition
CIm_log <- tcrit*sqrt(.006/34) # CI for multiplication

DT3_anovaplot_log$CI <- ifelse(DT3_anovaplot_log$operation == 'ADD', CIa_log, CIm_log)
# > levels(DT3_anovaplot_log$problemsize)
# [1] "1" "2"
# Not sure why this happened, need to relabel psize
DT3_anovaplot_log$problemsize<-factor(DT3_anovaplot_log$problemsize, labels=c("S","L"))

pd<-position_dodge(0.1)

pDigitWord_DT3anova_log <- ggplot(DT3_anovaplot_log, aes(x=problemsize, y=RTlog, color=format, linetype=format, group = format)) + geom_point(position=pd, aes(shape=format)) + geom_errorbar(aes(ymin=RTlog-CI, ymax=RTlog+CI), colour="black", width=.1, linetype = 1, position=pd) + geom_line(position=pd) + ylim(log(900),log(2500)) + scale_color_manual(values=c('black','black')) + labs(x='Problem Size', y='Mean Response Time (log-transformed)') + facet_grid(.~operation) + apatheme
ggsave(file='~/Documents/R/Pseudohomophone/pDigitWord_DT3anova_log.png')


# extra models -----
# test size x fmt on d-p for nonretrieval
m322dpDT3mt_REML <- lmer(RTlog ~ c_psizex*fmtP*nretr + (1 + c_psizex + fmtP + nretr|sub) + (1|item), data=subset(DT3_trim, typexy=='NTIE' & format!= "WORD" & operation == "MUL"))
# no size x fmt on nonretrieval trials

# test size x fmt on d-w for nonretrieval
m322dwDT3atr_REML <- lmer(RT ~ c_psizex*fmtW*nretr + (1 + c_psizex + fmtW + nretr|sub) + (1|item), data=subset(DT3_trim, typexy=='NTIE' & format!= "PSEUDO" & operation == "ADD"))
# ns size x fmt on nonretrieval trials (boarderline on retrieval trials)

### plot retr% and nretr% -----
# get rid of strategy = NA trials.  Only a handful.
retr_count <- as.data.frame(count(subset(DT3_trim, typexy == 'NTIE' & retr != 'NA'), vars=c('operation', 'format', 'psizex', 'retr')))

# rename columns
retr_count$count <- retr_count$freq
retr_count$freq <- NULL

# calculate number of trials irrespective of strategy
retr_freq <- as.data.frame(count(subset(DT3_trim, typexy == 'NTIE' & retr != 'NA'), vars=c('operation', 'format', 'psizex')))

# rename columns
retr_freq$total <- retr_freq$freq
retr_freq$freq <- NULL

retr_freq <- merge(retr_count, retr_freq, by = c('operation', 'format', 'psizex'))

retr_freq$perc <- retr_freq$count/retr_freq$total

# ggplot(retr_freq, aes(x = psizex, y = perc)) + geom_line(aes(color = retr)) + facet_grid(operation ~ format)

(pRETR_freq <- ggplot(retr_freq, aes(x = psizex, y = perc, color = retr, linetype = retr, group = retr)) + geom_line() + scale_color_manual(values=c('black','black')) + labs(x='Problem Size', y='Proportions') + scale_y_continuous(breaks=c(.5, 1)) + facet_grid(operation ~ format) + apatheme)
ggsave(file='~/Documents/R/Pseudohomophone/pRETR_freq.png')


# count how many participants would be deleted if apply ANOVA
emptycell_add_count <- as.data.frame(count(subset(DT3_trim, operation == 'ADD' & typexy == 'NTIE' & retr != 'NA'), vars=c('sub', 'format', 'problemsize', 'retr')))

emptycell_mul_count <- as.data.frame(count(subset(DT3_trim, operation == 'MUL' & typexy == 'NTIE' & retr != 'NA'), vars=c('sub', 'format', 'problemsize', 'retr')))
# no empty cell, but some cells with very few observations

# effect sizes for size x fmt -----
cellmeans <- with(subset(DT3_trim, typexy == 'NTIE'), aggregate(x = list(Y = RT), by = list(A = format, B = retr, C = problemsize, D = operation), FUN = mean))
