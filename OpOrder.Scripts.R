# deleted multiple lines here, data exploration including histogram and spagh plot

# deleted 100+ lines here, models with RT trimmed at 500ms

#------------------------------------------------------------------------------------#
# multiplication & addition 2016-11-17
# maximal str
model400DT2.RT250 <- lmer(RTlog~c_psize*odct*op + c_psize*tpct*op + (1+c_psize+odct+tpct+op|sub) + (1|item), REML=F, data=DT2.RT250)
# failed to converge

# max str, reduced sample
model401DT2.RT250 <- lmer(RTlog~c_psize*odct*op + c_psize*tpct*op + (1+c_psize+odct+tpct+op|sub) + (1|item), REML=F, data=DT2_50.RT250)
# failed to converge

# maximal str, zero corelation
model410DT2.RT250 <- lmer(RTlog~c_psize*odct*op + c_psize*tpct*op + (1|sub) + (0+c_psize|sub) + (0+odct+tpct|sub) + (0+op|sub) + (1|item), REML=F, data=DT2.RT250)
# a bizzar thing about this model is it generates two random effects for op: opadd and opmul, only when I impose zero correlation on int and op.  Only opmul would be generated if I don't impose this restriction.  Is this expected whenever an IV is a factor?
# converged
display(model410DT2.RT250)
# AIC = 1341.6, DIC = 1297.6, deviance = 1297.6


# partial correlation, between int and psize
model420DT2.RT250 <- lmer(RTlog~c_psize*odct*op + c_psize*tpct*op + (1+c_psize|sub) + (0+odct+tpct|sub) + (0+op|sub) + (1|item), REML=F, data=DT2.work.RT250)
# failed to converge

# partial correlation, between int, psize and two contrasts for type and order
model430DT2.RT250 <- lmer(RTlog~c_psize*odct*op + c_psize*tpct*op + (1+c_psize+odct+tpct|sub) + (0+op|sub) + (1|item), REML=F, data=DT2.work.RT250)
# failed to converge

# maximal str, separate operation
model401DT2M.RT250 <- lmer(RTlog~c_psize*odct + c_psize*tpct + (1+c_psize+odct+tpct|sub) + (1|item), REML=F, data=subset(DT2.work.RT250,op='add'))
model401DT2A.RT250 <- lmer(RTlog~c_psize*odct + c_psize*tpct + (1+c_psize+odct+tpct|sub) + (1|item), REML=F, data=subset(DT2.work.RT250,op='mul'))
# both converged

# diagnostic of model410DT2.RT250 begins
# fitted ~ random effects
# (p.pred.psize.m410DT2 <-ggplot(DT2.work.RT250, aes(x=DT2.work.RT250$c_psize,y=fitted(model410DT2.RT250),group=sub)) + aes(alpha=alpha) + guides(alpha=F) + geom_line() + facet_grid(op~type))
# (p.pred.type.m410DT2 <-ggplot(DT2.work.RT250, aes(x=DT2.work.RT250$type,y=fitted(model410DT2.RT250),group=sub)) + aes(alpha=alpha) + guides(alpha=F) + geom_line() + facet_grid(op~order))
# (p.pred.op.m410DT2 <-ggplot(DT2.work.RT250, aes(x=DT2.work.RT250$op,y=fitted(model410DT2.RT250),group=sub)) + aes(alpha=alpha) + guides(alpha=F) + geom_point() + facet_grid(type~order))
# # the last plot for op looks weird, no fitted value for mul&LS or Tie&add, not sure why.  But It seems intuitively justifiable to set op as a random factor

# Level 1 diagnostic
# Homo
# residual ~ sub
plot(model410DT2.RT250,sub~resid(., scaled=TRUE))
# outliers on both sides

tmp <- resid(model410DT2.RT250)
DT2.RT250$resid<-tmp

tmp_a <- sort(tmp) # ascending
tmp_d <- sort(tmp, decreasing= T)# descending

# delete up to 50 cases on both head and tail and do sensitivity analysis
DT2_50.RT250 <- subset(DT2.RT250, resid > tmp_a[25] & resid < tmp_d[25])

######## sensitivity analysis ###########
# maximal str, zero corelation, with reduced sample
model411DT2.RT250 <- lmer(RTlog~c_psize*odct*op + c_psize*tpct*op + (1|sub) + (0+c_psize|sub) + (0+odct+tpct|sub) + (0+op|sub) + (1|item), REML=F, data=DT2_50.RT250)

plot(model411DT2.RT250,sub~resid(., scaled=TRUE))

# compare results
# random effects
randsens.m410DT2 <- cbind(as.data.frame(VarCorr(model410DT2.RT250))[,c(1:3,5)], as.data.frame(VarCorr(model411DT2.RT250))[,5])
colnames(randsens.m410DT2) <- c('grp','var1', 'var2', 'Full','50OutliersDel')

#fixed effect
summ.model410DT2 <- summary(model410DT2.RT250)
summ.model411DT2 <- summary(model411DT2.RT250)
fixedsens.m410DT2 <-cbind(as.data.frame(summ.model410DT2$coefficients)[,c(1,3)], as.data.frame(summ.model411DT2$coefficients)[,c(1,3)])
colnames(fixedsens.m410DT2) <- c('full.est','full.t','50OutliersDel.est','50OutliersDel.t')

write.table(randsens.m410DT2, file='~/Documents/R/OpOrder/randsens.m410DT2.txt', quote=F, sep="\t")
write.table(fixedsens.m410DT2, file='~/Documents/R/OpOrder/fixedsens.m410DT2.txt', quote=F, sep="\t")
##### sensitivity analysis end #######

# residual ~ pred
plot(model410DT2.RT250)
# no discernable relationship between residual and fitted value, suggest appropriate model specification

# residual ~ psize
ggplot(model410DT2.RT250, aes(DT2.work.RT250$psize,residuals(model410DT2.RT250,scaled=T))) + geom_point() + geom_smooth(method = "loess", size=1.5)
ggsave('~/Documents/R/OpOrder/add.HomoSizeResid.png')

plot(model410DT2.RT250,type ~ resid(., scaled=TRUE))

plot(model410DT2.RT250,order ~ resid(., scaled=TRUE))

plot(model410DT2.RT250,op ~ resid(., scaled=TRUE))

#residual & predicted
xyplot(resid(model410DT2.RT250)~fitted(model410DT2.RT250))
# No discernable outliers

#normality of residuals
hist(residuals(model410DT2.RT250,scaled=T))
qqnorm(residuals(model410DT2.RT250,scaled=T))
# not perfectly normal

# Level 2
# extract random effects. ranef is structured as a list; at [1], it is the random int of item; at [2], random int & slps of sub.  See 2015-11-09 12:20 PM in OneNote
Level2_1.m410DT2 <- as.data.frame(ranef(model410DT2.RT250)[1])
Level2_2.m410DT2 <- as.data.frame(ranef(model410DT2.RT250)[2])
names(Level2_1.m410DT2) <- c("uIntItem")
# random int of items are grouped on items
names(Level2_2.m410DT2) <- c('uAdd','uMul', 'uodct', 'utpct', 'uSize','uIntSub')
# [ ] I don't quite understand why two random effects for two levels of op. Only one random effect would be generated if I allow correlation between op and int

# linearity & normality
hist(Level2_1.m410DT2$uIntItem)
ggpairs(data=Level2_2.m410DT2, title = "Relation between random effects")
# signs of correlation between int, psize, and two contrasts (odct, tpct)

# diagnostic of model410DT2.RT250 ends

# Analysis proper, turn off REML=F
model410DT2.ML <- lmer(RTlog~c_psize*odct*op + c_psize*tpct*op + (1|sub) + (0+c_psize|sub) + (0+odct+tpct|sub) + (0+op|sub) + (1|item), data=DT2.work.RT250)
# failed to converge, try the same analysis with reduced sample (50 outliers deleted).
model411DT2.ML <- lmer(RTlog~c_psize*odct*op + c_psize*tpct*op + (1|sub) + (0+c_psize|sub) + (0+odct+tpct|sub) + (0+op|sub) + (1|item), data=DT2_50.work.RT250)
# converged successfully. Even though the reduced sample seem to have no impact on estimation, these outliers seem to affect the algorithm.  The max model with reduced sample would still not converge (model401DT2.RT250)
display(model411DT2.ML)

# reverse coding of op, to test order by operation for addition 2016-11-25
# recode op, using mul as the ref group
DT2_50.work.RT250$op_mul <- factor(with(DT2_50.work.RT250,ifelse(DT2_50.work.RT250$op == 'mul', 0, 1)), levels=c(0,1), labels=c('mul','add'))
model412DT2.ML <- lmer(RTlog~c_psize*odct*op_mul + c_psize*tpct*op_mul + (1|sub) + (0+c_psize|sub) + (0+odct+tpct|sub) + (0+op_mul|sub) + (1|item), REML = F, data=DT2_50.work.RT250)
# failed to converge!  But the only difference between model411 and model412 is dummy coding of operation! Another sign of how finicky the model fitting process is.
model412DT2 <- lmer(RTlog~c_psize*odct*op_mul + c_psize*tpct*op_mul + (1|sub) + (0+c_psize|sub) + (0+odct+tpct|sub) + (0+op_mul|sub) + (1|item), REML = F, data=DT2_50.work.RT250)
# with REML turned off, however, the model converged, but gave "rescale var" warning, I don't understand ... regardless, the coefficient of operandXoperation is significant
display(model412DT2)

# export results of model411DT2.ML
write.table(as.data.frame(summary(model411DT2.ML)$coefficients),"~/Documents/R/OpOrder/fixef.m411DT2ML.txt", quote=F, sep="\t")

# I need to plug in c_psize into the equation when calculating point estimate, then plot againg psize.  To do that, I need to extract each pair of c_psize & psize
psizeAdd <- sort(unique(DT2add$psize))
c_psizeAdd <- sort(unique(DT2add$c_psize))

psizeMul <- sort(unique(DT2mul$psize))
c_psizeMul <- sort(unique(DT2mul$c_psize))
psize <- as.data.frame(cbind(psizeAdd, c_psizeAdd, psizeMul, c_psizeMul), col.names=c('psizeAdd', 'c_psizeAdd', 'psizeMul', 'c_psizeMul'))
write.table(psize,"~/Documents/R/OpOrder/psize.txt", quote=F, sep="\t")

# write.table(as.data.frame(summary(model411DT2.ML)$varcor),".txt", quote=F, sep="\t")

#------------------------------------------------------------------------------------#
# deleted multiple lines on add, mul modeled separately

###### use product as psize for addition, re-plot 2017-01-04 ##########
#scale psizex after concatecate
DT2.RT250$c_psizex <- scale(DT2.RT250$psizex, center=T, scale=T)
# This indeed improved model convergence!!!

# maximal str, zero corelation
model410DT2.psizex.RT250 <- lmer(RTlog~c_psizex*odct*op + c_psizex*tpct*op + (1|sub) + (0+c_psizex|sub) + (0+odct+tpct|sub) + (0+op|sub) + (1|item), REML=F, data=DT2.RT250)
# converged
display(model410DT2.RT250)

# Analysis proper, turn off REML=F
model410DT2.psizex.ML <- lmer(RTlog~c_psizex*odct*op + c_psizex*tpct*op + (1|sub) + (0+c_psizex|sub) + (0+odct+tpct|sub) + (0+op|sub) + (1|item), data=DT2.RT250)
# failed to converge, try the same analysis with reduced sample (50 outliers deleted).
model411DT2.psizex.ML <- lmer(RTlog~c_psizex*odct*op + c_psizex*tpct*op + (1|sub) + (0+c_psizex|sub) + (0+odct+tpct|sub) + (0+op|sub) + (1|item), data=DT2_50.RT250)

# export results and plot
write.table(as.data.frame(summary(model411DT2.psizex.ML)$coefficients),"~/Documents/R/OpOrder/fixef.m411DT2psizexML.txt", quote=F, sep="\t")

# I need to plug in c_psize into the equation when calculating point estimate, then plot againg psize.  To do that, I need to extract each pair of c_psize & psize
psizex <- sort(unique(DT2.RT250$psizex))
c_psizex <- sort(unique(DT2.RT250$c_psizex))

psizex <- as.data.frame(cbind(psizex, c_psizex), col.names=c('psizex', 'c_psizex'))
write.table(psizex,"~/Documents/R/OpOrder/psizex.txt", quote=F, sep="\t")


###### [obsolete] START: use product as psize for addition, re-plot 2016-12-05  ##########
# maximal str
model400DT2.psizex.RT250<- lmer(RTlog~c_psizex*odct*op + c_psizex*tpct*op + (1+c_psizex+odct+tpct+op|sub) + (1|item), REML=F, data=DT2.RT250)
# failed to converge

# max str, reduced sample
model401DT2.psizex.RT250 <- lmer(RTlog~c_psizex*odct*op + c_psizex*tpct*op + (1+c_psizex+odct+tpct+op|sub) + (1|item), REML=F, data=DT2_50.RT250)
# failed to converge

# maximal str, zero corelation
model410DT2.psizex.RT250 <- lmer(RTlog~c_psizex*odct*op + c_psizex*tpct*op + (1|sub) + (0+c_psizex|sub) + (0+odct+tpct|sub) + (0+op|sub) + (1|item), REML=F, data=DT2.RT250)
# failed to converge!
display(model410DT2.RT250)

# maximal str, zero corelation, reduced sample
# note that outliers were determined using information from model410DT2.RT250
model411DT2.psizex.RT250 <- lmer(RTlog~c_psizex*odct*op + c_psizex*tpct*op + (1|sub) + (0+c_psizex|sub) + (0+odct+tpct|sub) + (0+op|sub) + (1|item), REML=F, data=DT2_50.RT250)
# failed to converge

# drop item as a random effect, zero correlation, full sample
model310DT2.psizex.RT250 <- lmer(RTlog~c_psizex*odct*op + c_psizex*tpct*op + (1|sub) + (0+c_psizex|sub) + (0+odct+tpct|sub) + (0+op|sub), REML=F, data=DT2.RT250)
# converged

# drop item as a random effect, full sample
model300DT2.psizex.RT250 <- lmer(RTlog~c_psizex*odct*op + c_psizex*tpct*op + (1+c_psizex+odct+tpct+op|sub), REML=F, data=DT2.RT250)
# failed to converge

# drop item as a random effect, reduced sample
model301DT2.psizex.RT250 <- lmer(RTlog~c_psizex*odct*op + c_psizex*tpct*op + (1+c_psizex+odct+tpct+op|sub), REML=F, data=DT2_50.RT250)
# failed to converge

# I need to plot the results of model310DT2.psizex.RT250 to see if the pattern remain the same as previously when addition psize was coded as sum

# REML=F turned off, zero correlation, reduced sample, model310DT2.psizex.RT250.ML on full sample failed to converge
model311DT2.psizex.RT250.ML <- lmer(RTlog~c_psizex*odct*op + c_psizex*tpct*op + (1|sub) + (0+c_psizex|sub) + (0+odct+tpct|sub) + (0+op|sub), data=DT2_50.RT250)

###### [obsolete] END: use product as psize for addition, re-plot 2016-12-05 ##########