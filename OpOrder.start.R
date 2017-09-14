library(pacman)
p_load(tcltk2,lme4,arm,plyr,ggplot2,reshape2, lattice, GGally)

tclTaskSchedule(1000, {options(prompt=paste(Sys.time(),"> "))}, id = "ticktock", redo = TRUE)

###### multiplication #########
DT2mul <- read.table("~/Documents/R/OpOrder/OpOrderMul.csv",header=T, sep=',')

##### data prep ##### 
# rename columns, make them the same as Liu97 
names(DT2mul)[1]<-'sub'
names(DT2mul)[2]<-"trial"
names(DT2mul)[3]<-'Lvalue'
names(DT2mul)[4]<-'Rvalue'
names(DT2mul)[5]<-"psize"
names(DT2mul)[6]<-"RT"
names(DT2mul)[7]<-"response"
# getting rid of S in front of responses
DT2mul$responseN = as.numeric(sapply (DT2mul$response, function (x) {gsub('S','',x)}))

# try rescale psize.  It would have made a difference whether I rescale psize in DT2mul or DT2mul.work.RT500, or rescale psize separately.  But I am not going to interpret raw coefficients, so as long as I keep the correspondence between psize and c_psize, it's fine
DT2mul$c_psize <- scale(DT2mul$psize, center=T, scale=T)

# add columns for accuracy, type, order
# accuracy
DT2mul$correct <- factor(with(DT2mul,ifelse(DT2mul$psize == DT2mul$responseN, 1, 0)))

# deleted script for coding tie, five, and reg with two dummy ... see the backup version

# order
for (i in 1:length(DT2mul$sub)) {
  if (DT2mul$Lvalue[i] < DT2mul$Rvalue[i]) {
    DT2mul$order[i] <- 0
  }
  else if (DT2mul$Lvalue[i] > DT2mul$Rvalue[i]) {
    DT2mul$order[i] <- 1
  }
  else DT2mul$order[i] <- NA
}
DT2mul$order <- factor(DT2mul$order, levels = c(0,1), labels=c('SL', 'LS'))


# type, tie vs. non-tie
for (i in 1:length(DT2mul$sub)) {
  if (DT2mul$Lvalue[i] == DT2mul$Rvalue[i]) {
    DT2mul$type[i] <- 1
  }
  else DT2mul$type[i] <- 2
}
DT2mul$type <- factor(DT2mul$type, levels = c(1,2), labels=c('Tie','NTie'))

# type_order, contrast coding. odct: contrast order; tpct: contrast type
# making sure these columns have been deleted before running the code
DT2mul$odct <- NULL
DT2mul$tpct <- NULL

for (i in 1:length(DT2mul$sub)) {
  if (DT2mul$type[i] == 'NTie') {
    if (DT2mul$Lvalue[i] < DT2mul$Rvalue[i]) {
      DT2mul$odct[i] <- 1
      DT2mul$tpct[i] <- .5
    }
    else {
      DT2mul$odct[i] <- -1
      DT2mul$tpct[i] <- .5
    }
  }
  
  else {
    DT2mul$odct[i] <- 0
    DT2mul$tpct[i] <- -1
  }
}


# add random effects of item, need to create item index first. 
# name items, 2x3 is different from 3x2
DT2mul$item <- NULL
for (i in 1:length(DT2mul$sub)) {
  DT2mul$item[i] <- paste(DT2mul$Lvalue[i], DT2mul$Rvalue[i], 'm', sep='')
}

# convert sub and item to factor if they are not already factors
if (!is.factor(DT2mul$sub)) {
  DT2mul$sub <- factor(DT2mul$sub)
  }
if (!is.factor(DT2mul$item)) {
  DT2mul$item <- factor(DT2mul$item)
}

DT2mul$op = 'mul'

# add alpha to avoid overplotting, turned it on when analyzing mul and add separately
DT2mul = ddply(DT2mul, .(sub), function(x){x$alpha = ifelse(runif(n = 1) > 0.8, 0.9, 0.1)
x})

# trimming RT >= 0 (no trimming)
# DT2mul.work <- subset(DT2mul, RT >= 0 & correct == 1)
# DT2mul.work$RTlog <- log(DT2mul.work$RT)

# trimming RT >= 250
DT2mul.RT250 <- subset(DT2mul, RT >= 250 & correct == 1)
DT2mul.RT250$RTlog <- log(DT2mul.RT250$RT)

################### addition #########################
DT2add <- read.table("~/Documents/R/OpOrder/OpOrderAdd.csv",header=T, sep=',')

##### data prep ##### 
# rename columns, make them the same as Liu97 
names(DT2add)[1]<-'sub'
names(DT2add)[2]<-"trial"
names(DT2add)[3]<-'Lvalue'
names(DT2add)[4]<-'Rvalue'
names(DT2add)[5]<-"psize"
names(DT2add)[6]<-"RT"
names(DT2add)[7]<-"response"
# getting rid of S in front of responses
DT2add$responseN = as.numeric(sapply (DT2add$response, function (x) {gsub('S','',x)}))

# try rescale psize.  It would have made a difference whether I rescale psize in DT2add or DT2add.work.RT500, or rescale psize separately.  But I am not going to interpret raw coefficients, so as long as I keep the correspondence between psize and c_psize, it's fine
DT2add$c_psize <- scale(DT2add$psize, center=T, scale=T)

# add columns for accuracy, type, order
# accuracy
DT2add$correct <- factor(with(DT2add,ifelse(DT2add$psize == DT2add$responseN, 1, 0)))

# type, turn it on when analyzing mul and add separately
# for (i in 1:length(DT2add$sub)) {
#   if (DT2add$Lvalue[i] == DT2add$Rvalue[i]) {
#     DT2add$type[i] <- 1
#   }
#   else DT2add$type[i] <- 2
# }
# DT2add$type <- factor(DT2add$type, levels = c(1,2), labels=c('Tie','Reg'))

# order
for (i in 1:length(DT2add$sub)) {
  if (DT2add$Lvalue[i] < DT2add$Rvalue[i]) {
    DT2add$order[i] <- 0
  }
  else if (DT2add$Lvalue[i] > DT2add$Rvalue[i]) {
    DT2add$order[i] <- 1
  }
  else DT2add$order[i] <- NA
}
DT2add$order <- factor(DT2add$order, levels = c(0,1), labels=c('SL', 'LS'))


# type, tie vs. non-tie
for (i in 1:length(DT2add$sub)) {
  if (DT2add$Lvalue[i] == DT2add$Rvalue[i]) {
    DT2add$type[i] <- 1
  }
  else DT2add$type[i] <- 2
}
DT2add$type <- factor(DT2add$type, levels = c(1,2), labels=c('Tie','NTie'))

# type_order, contrast coding. odct: contrast order; tpct: contrast type
# making sure these columns have been deleted before running the code
DT2add$odct <- NULL
DT2add$tpct <- NULL

for (i in 1:length(DT2add$sub)) {
  if (DT2add$type[i] == 'NTie') {
    if (DT2add$Lvalue[i] < DT2add$Rvalue[i]) {
      DT2add$odct[i] <- 1
      DT2add$tpct[i] <- .5
    }
    else {
      DT2add$odct[i] <- -1
      DT2add$tpct[i] <- .5
    }
  }
  
  else {
    DT2add$odct[i] <- 0
    DT2add$tpct[i] <- -1
  }
}


# add random effects of item, need to create item index first. 
# name items, 2+3 is different from 3+2
DT2add$item <- NULL
for (i in 1:length(DT2add$sub)) {
  DT2add$item[i] <- paste(DT2add$Lvalue[i], DT2add$Rvalue[i], 'a', sep='')
}

# convert sub and item to factor if they are not already factors
if (!is.factor(DT2add$sub)) {
  DT2add$sub <- factor(DT2add$sub)
  }
if (!is.factor(DT2add$item)) {
  DT2add$item <- factor(DT2add$item)
  }

DT2add$op = 'add'

# add alpha to avoid overplotting
DT2add = ddply(DT2add, .(sub), function(x){x$alpha = ifelse(runif(n = 1) > 0.8, 0.9, 0.1)
x})

# trimming RT >= 250
DT2add.RT250 <- subset(DT2add, RT >= 250 & correct==1)
DT2add.RT250$RTlog <- log(DT2add.RT250$RT)

# create a copy of DT2add.work.RT250 and recale the psize as product, then centre it (2016-12-05)
# deleted a few lines here
# 2017-01-02 instead of creating a second copy for coding addtion psize as product, create a second column in the same file. Correspondingly, create a second column in mul data (psizex, c_psizex), so that I can concatenate add and mul files. I need to delete some outliers, but I also need to select these outliers based on a converged model...
DT2add.RT250$psizex <- DT2add.RT250$Lvalue * DT2add.RT250$Rvalue
DT2add.RT250$c_psizex <- scale(DT2add.RT250$psizex, center=T, scale=T)

DT2mul.RT250$psizex <- DT2mul.RT250$psize
DT2mul.RT250$c_psizex <- DT2mul.RT250$c_psize

######### addition & multiplication ########
DT2.RT250 <- rbind(DT2mul.RT250, DT2add.RT250)
# code op with add as the ref group
DT2.RT250$op <- factor(with(DT2.RT250,ifelse(DT2.RT250$op == 'mul', 1, 0)), levels=c(0,1), labels=c("add","mul"))
# recode op to test order by operation for addition, mul as the reference group
# DT2.work.RT250$op_mul <- factor(with(DT2.work.RT250,ifelse(DT2.work.RT250$op == 'mul', 0, 1)), levels=c(0,1), labels=c("mul","add"))