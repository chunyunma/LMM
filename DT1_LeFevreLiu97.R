##### Appendix   R scripts for LeFevre & Liu 97 #####

data <- read.table("~/Documents/R/LiuDataFull.txt",header=T)
data.work <- unique(subset(data, correct==1 & type!='Zeros' & type!='TimeOne')) #correct trials only, no duplicates
# create dummy variables
data.work$typeD1 <- factor(with(data.work,ifelse(data.work$type=="Tie",1,0)))
data.work$typeD1<-factor(data.work$typeD1,levels=c(0,1), labels=c("Reg","Tie"))
data.work$typeD2 <- factor(with(data.work,ifelse(data.work$type=="Fives",1,0)))
data.work$typeD2<-factor(data.work$typeD2,levels=c(0,1), labels=c("Reg","Fives"))
# trimming RT >= 500
data.work.RT500 <- subset(data.work,RT >= 500)
# add alpha to avoid overplotting
data.work.RT500 = ddply(data.work.RT500, .(sub), function(x){x$alpha = ifelse(runif(n = 1) > 0.9, 1, 0.1)
+ x
+ })

#Data exploration, without running any model, using RT, not RTlog
#spaghetti template, tspag
tspag<-ggplot(data.work.RT500,aes(x=product,y=RT))+ geom_point(size=0)
# plot size by probType by culture
# add smoother to each individual line
bwspag <- tspag + aes(group=factor(sub)) + guides(colour=F)
bwspag + facet_grid(culture~type) + geom_line(stat="smooth",method = "loess", aes(alpha=alpha)) + stat_summary(aes(group = 1), geom = "point", fun.y = mean, shape = 17, size = 2, colour='red') + guides(alpha=F)
# Will polish this plot later, but it has the key elements I want: individual spaghetti and group mean
# Crude way of telling how many random effects I need (three)


# test models
###### model.RT500.D1D2 ######
model.RT500.D1D2 <-lmer(RT~product + typeD1 + typeD2 + product:typeD1 + product:typeD2 + culture + product:culture + typeD1:culture + typeD2:culture + product:typeD1:culture + product:typeD2:culture + (1 + typeD1 + typeD2|sub), data=data.work.RT500)
# also run model.work.D1D2, model.work.D1, model.work.D2 etc
# What random effects to include? Plot spaghetti of predicted values for different models
p.RT500.D1D2<-ggplot(data.work.RT500, aes(x=data.work.RT500$product,y=fitted(model.RT500.D1D2),group=sub))
p.RT500.D1D2 + aes(alpha=alpha) + guides(alpha=F) + geom_line()


#transform RT in data.work.RT500, plot histogram
data.work.RT500$RTlog <- log(data.work.RT500$RT)
par(mfrow=c(2,2))
hist(data.work.RT500$RT)
hist(data.work.RT500$RTlog)
qqnorm(data.work.RT500$RT)
qqnorm(data.work.RT500$RTlog)

##### model.RT500.log #######
model.RT500.log <-lmer(RTlog~product + typeD1 + typeD2 + product:typeD1 + product:typeD2 + culture + product:culture + typeD1:culture + typeD2:culture + product:typeD1:culture + product:typeD2:culture + (1 + product + typeD1 + typeD2|sub), data=data.work.RT500)

# Should I include three random effects ? spaghetti of predicted values vs. each random effect
#size
p.PredSize.RT500.log<-ggplot(data.work.RT500, aes(x=data.work.RT500$product,y=fitted(model.RT500.log),group=sub))
p.PredSize.RT500.log + aes(alpha=alpha) + guides(alpha=FALSE) + geom_line()
ggsave(file="~/Documents/R/Diagnostic/RT500log/pred|size.png")
#D1
p.PredD1.RT500.log<-ggplot(data.work.RT500, aes(x=data.work.RT500$typeD1,y=fitted(model.RT500.log),group=sub))
p.PredD1.RT500.log + aes(alpha=alpha) + guides(alpha=FALSE) + geom_line()
ggsave(file="~/Documents/R/Diagnostic/RT500log/pred|D1.RT500log.png")
#D2
p.PredD2.RT500.log<-ggplot(data.work.RT500, aes(x=data.work.RT500$typeD2,y=fitted(model.RT500.log),group=sub))
p.PredD2.RT500.log + aes(alpha=alpha) + guides(alpha=FALSE) + geom_line()
ggsave(file="~/Documents/R/Diagnostic/RT500log/pred|D2.RT500log.png")

#normality of residuals
# I can use par(mfrow=c(2,2)) to put normality plots for RTlog and RT in one graph
hist(residuals(model.RT500.log,scaled=T))
qqnorm(residuals(model.RT500.log,scaled=T))

#Level 1 diagnostic
# Homo
plot(model.RT500.log,sub ~ resid(., scaled=TRUE))
ggplot(data.work.RT500, aes(data.work.RT500$product,residuals(model.RT500.log,scaled=T))) + geom_point() + geom_smooth(method = "loess", size=1.5)
ggsave(file="HomoSizeResid.RT500log.png")
#residual & predicted
xyplot(resid(model.RT500.log)~fitted(model.RT500.log))
# No discernable outliers

# Level 2
Level2.RT500.log <- ranef(model.RT500.log)
p_load(GGally)
Level2.RT500.log <- as.data.frame(ranef(model.RT500.log)[[1]])
names(Level2.RT500.log) <- c("u0","uSize","uD1","uD2")
# linearity & normality
ggpairs(data=Level2.RT500.log, title = "Relation between random effects")
# Merge level 2 back to data
data.work.RT500.L2log <- merge(data.work.RT500, Level2.RT500.log, by.x="sub", by.y="row.names")
#only against level-2 predictors
boxplot(u0~culture, data=data.work.RT500.L2log)
boxplot(uSize~culture, data=data.work.RT500.L2log)
boxplot(uD1~culture, data=data.work.RT500.L2log)
boxplot(uD2~culture, data=data.work.RT500.L2log)
# A couple potential outliers, flag maxima/minima on random effects
u0.minima <- unique(data.work.RT500.L2log[order(data.work.RT500.L2log$u0),]$sub)[1:1]
# [1] 22
uSize.maxima <- unique(data.work.RT500.L2log[order(data.work.RT500.L2log$uSize,decreasing=T),]$sub)[1:1]
# [1] 17
uD1.maxima <- unique(data.work.RT500.L2log[order(data.work.RT500.L2log$uD1,decreasing=T),]$sub)[1:1]
uD1.minima <- unique(data.work.RT500.L2log[order(data.work.RT500.L2log$uD1),]$sub)[1:2]
uD1.minima
# [1] 17 38
uD1.maxima
# [1] 22
uD2.minima <- unique(data.work.RT500.L2log[order(data.work.RT500.L2log$uD2),]$sub)[1:2]
uD2.maxima <- unique(data.work.RT500.L2log[order(data.work.RT500.L2log$uD2,decreasing=T),]$sub)[1:1]
uD2.maxima
# [1] 22
uD2.minima
# [1] 38 17

#sensitivity analysis
model.RT500.log_17 <-lmer(RTlog~product + typeD1 + typeD2 + product:typeD1 + product:typeD2 + culture + product:culture + typeD1:culture + typeD2:culture + product:typeD1:culture + product:typeD2:culture + (1 + product + typeD1 + typeD2|sub), data=subset(data.work.RT500,sub!='17'))
model.RT500.log_22 <-lmer(RTlog~product + typeD1 + typeD2 + product:typeD1 + product:typeD2 + culture + product:culture + typeD1:culture + typeD2:culture + product:typeD1:culture + product:typeD2:culture + (1 + product + typeD1 + typeD2|sub), data=subset(data.work.RT500,sub!='22'))
model.RT500.log_38 <-lmer(RTlog~product + typeD1 + typeD2 + product:typeD1 + product:typeD2 + culture + product:culture + typeD1:culture + typeD2:culture + product:typeD1:culture + product:typeD2:culture + (1 + product + typeD1 + typeD2|sub), data=subset(data.work.RT500,sub!='38'))

# compare results
# random effects
random.sensitivity <- cbind(
  + as.data.frame(VarCorr(model.RT500.log))[,c(1:3,5)],
  + as.data.frame(VarCorr(model.RT500.log_17))[,5],
  + as.data.frame(VarCorr(model.RT500.log_22))[,5],
  + as.data.frame(VarCorr(model.RT500.log_38))[,5])
colnames(random.sensitivity) <- c('grp','var1', 'var2', 'full','full_17','full_22','full_38')

#fixed effects
summ.model.RT500.log <- summary(model.RT500.log)
summ.model.RT500.log_17 <- summary(model.RT500.log_17)
summ.model.RT500.log_22 <- summary(model.RT500.log_22)
summ.model.RT500.log_38 <- summary(model.RT500.log_38)
fixed.sensitivity<-cbind(
  + as.data.frame(summ.model.RT500.log$coefficients)[,c(1,3)],
  + as.data.frame(summ.model.RT500.log_17$coefficients)[,c(1,3)],
  + as.data.frame(summ.model.RT500.log_22$coefficients)[,c(1,3)],
  + as.data.frame(summ.model.RT500.log_38$coefficients)[,c(1,3)])
colnames(fixed.sensitivity) <- c('full.est','full.t','full_17.est','full_17.t', 'full_22.est','full_22.t', 'full_38.est','full_38.t')

write.table(random.sensitivity, file='~/Documents/R/Diagnostic/RT500log/rand.sensitivity.txt', quote=F, sep="\t")
write.table(fixed.sensitivity, file='~/Documents/R/Diagnostic/RT500log/fixef.sensitivity.txt', quote=F, sep="\t")

#Analysis proper
# export results of model.RT500.log
write.table(as.data.frame(summ.model.RT500.log$coefficients),"~/Documents/R/Diagnostic/RT500log/fixef.RT500.log.txt", quote=F, sep="\t")
write.table(as.data.frame(summ.model.RT500.log$varcor),"~/Documents/R/Diagnostic/RT500log/ranef.RT500.log.txt", quote=F, sep="\t")
#Tabulate fixed effects of model.RT500.D1D2, model.RT500.D1. I could combine previous export with this step
#Double check which raw-RT models converged
summary(model.RT500.D1D2)$optinfo$conv$lme4$messages
NULL
summary(model.RT500.D1)$optinfo$conv$lme4$messages
[[1]]
[1] "Model is nearly unidentifiable: very large eigenvalue\n - Rescale variables?"
summary(model.RT500.D2)$optinfo$conv$lme4$messages
[1] "Model failed to converge with max|grad| = 0.00201718 (tol = 0.002, component 1)"
summ.model.RT500.D1D2 <- summary(model.RT500.D1D2)
summ.model.RT500.D1 <- summary(model.RT500.D1)
fixed.rawRT <- cbind(
  + as.data.frame(summ.model.RT500.D1D2$coefficients)[,1:3],
  + as.data.frame(summ.model.RT500.D1$coefficients)[,1:3])
colnames(fixed.rawRT) <- c('rawRT.D1D2.est','rawRT.D1D2.se','rawRT.D1D2.t', 'rawRT.D1.est','rawRT.D1.se','rawRT.D1.t')
write.table(fixed.rawRT,"~/Documents/R/Diagnostic/RT500/fixef.RT500.txt", quote=F, sep="\t")
ranef.rawRT <- cbind(
  + as.data.frame(VarCorr(model.RT500.D1D2))[,c(1:3,5)],
  + as.data.frame(VarCorr(model.RT500.D1))[,c(1:3,5)])
#two models have different randome effects, be carefull when combining
colnames(ran.rawRT) <- c('grp','var1', 'var2', 'RT500.D1D2','grp','var1', 'var2','RT500.D1')
write.table(ranef.rawRT,'~/Documents/R/Diagnostic/RT500/ranef.RT500.txt', quote=F, sep="\t")

#plot three-way of model.RT500.D1
p_load(effects)
ef.SizeTieCulture.RT500.D1 <- effect("product:typeD1:culture", model.RT500.D1)
plot(ef.SizeTieCulture.RT500.D1, multiline=T)

#plot two-way of model.RT500.D1
ef.FiveCulture.RT500.D1 <- effect("typeD2:culture", model.RT500.D1)
x <- as.data.frame(ef.FiveCulture.RT500.D1)
(p.2way.FiveCulture.RT500.D1 <- ggplot(x,aes(typeD2, fit)) + geom_bar(stat='identity') + facet_wrap(~ culture))
ggsave(file='~/Documents/R/Diagnostic/RT500/Twoway.FiveCulture.RT500.D1.png')

# visualize random slope of size
p.PredSize.RT500.D1<-ggplot(data.work.RT500, aes(x=data.work.RT500$product,y=fitted(model.RT500.D1),group=sub))
p.PredSize.RT500.D1 + aes(alpha=alpha) + guides(alpha=FALSE) + geom_line()
ggsave(file='~/Documents/R/Diagnostic/RT500/pred.RT500.D1|size.png')