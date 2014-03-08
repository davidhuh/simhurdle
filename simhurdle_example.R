library(devtools)
library(glmmADMB)


## Load function for simulating zero-altered data from fitted hurdle model estimates
source_url("https://raw.github.com/davidhuh/simhurdle/master/simhurdle.R")

### Example:

## Matrix of predictors
n.clust <- 30
n.id    <- 172
N.obs   <- n.id * n.clust

id      <- rep(1:n.id, each=n.clust)
day     <- rep(0:(n.clust-1), times=n.id)
dayofwk <- rep(rep(0:6, length.out=30), length.out=N.obs)

# randomized vector of covariate values
dmqsoc.cf <- c(-2.63, -2.19, -1.97, -1.75, -1.53, -1.31, -1.09, -0.87, -0.65,
               -0.43, -0.21, 0.01, 0.23, 0.45, 0.67, 0.88, 1.10, 1.32, 1.54)

dmqsoc.cf.prob <- c(0.006657019, 0.013892909, 0.020260492, 0.030969609, 0.017945007,
                    0.035600579, 0.067727931, 0.053835022, 0.064254703, 0.055571635,
                    0.042257598, 0.054124457, 0.102170767, 0.084515195, 0.096092619,
                    0.061360347, 0.083068017, 0.040520984, 0.069175109)

dmqsoc <- sample(dmqsoc.cf,
                 size = n.id,
                 replace = TRUE,
                 prob = dmqsoc.cf.prob)

# cyclical variables
cycAmp <- cos(2*pi*(dayofwk/7))
cycPha <- sin(2*pi*(dayofwk/7))

# full dummy variables
weekdayMon <- as.numeric(dayofwk==1)
weekdayTue <- as.numeric(dayofwk==2)
weekdayWed <- as.numeric(dayofwk==3)
weekdayThu <- as.numeric(dayofwk==4)
weekdayFri <- as.numeric(dayofwk==5)
weekdaySat <- as.numeric(dayofwk==6)

pred.df <- data.frame(id=id, day=day, dayofwk=dayofwk, dmqsoc=dmqsoc, cycPha=cycPha, cycAmp=cycAmp,
                      weekdayMon=weekdayMon, weekdayTue=weekdayTue, weekdayWed=weekdayWed,
                      weekdayThu=weekdayThu, weekdayFri=weekdayFri, weekdaySat=weekdaySat)

## Load fitted model estimates
load("modelfit.Rdata")

y.sim <- simulate.hurdle(data=pred.df,
                         formula.fe.bin= ~dmqsoc*(weekdayMon + weekdayTue + weekdayWed + weekdayThu + weekdayFri + weekdaySat),
                         formula.fe.cnt= ~dmqsoc*(cycPha + cycAmp),
                         formula.re.bin= ~1,
                         formula.re.cnt= ~cycPha + cycAmp,
                         coef.fe.bin=beta.bin,
                         coef.fe.cnt=beta.cnt,
                         vcov.fe.bin=vcov.beta.bin,
                         vcov.fe.cnt=vcov.beta.cnt,
                         vcov.re.bin=vcov.id.bin,
                         vcov.re.cnt=vcov.id.cnt,
                         dist.cnt="xtnb",
                         od.cnt=alpha.mu.cnt, link.bin="logit", link.cnt="log", idvar="id")

dash.sim.df <- cbind(pred.df, y.sim)


## incorporate missing data (assuming MCAR)
n.obs.orig <- 3455
n.obs.keep <- sort(sample(seq(1, N.obs), size=n.obs.orig, replace=FALSE))

dash.sim.df <- dash.sim.df[n.obs.keep, ]

## create variables necessary for analysis in glmmADMB
dash.sim.df <- within(dash.sim.df, {
  id <- factor(id)
  drinks <- y.sim
  
  dv.bin <- ifelse(y.sim==0, 0, 1)
  dv.cnt <- ifelse(y.sim==0, NA, y.sim)
})

dash.sim.bin.df <- dash.sim.df
dash.sim.cnt.df <- dash.sim.df[!is.na(dash.sim.df$dv.cnt), ]

##### Fit model to simulated data

bin.sim.model <- glmmadmb(dv.bin ~ dmqsoc*(weekdayMon + weekdayTue + weekdayWed + weekdayThu + weekdayFri + weekdaySat) + (1|id),
                      data=dash.sim.bin.df,
                      family="binomial",
                      verbose = FALSE)

summary(bin.sim.model)

cnt.sim.model <- glmmadmb(dv.cnt ~ dmqsoc*(cycPha + cycAmp) + (cycPha + cycAmp|id),
                      data=dash.sim.cnt.df,
                      family="truncnbinom",
                      #admb.opts=admbControl(shess=FALSE, noinit=FALSE),
                      verbose = FALSE)

summary(cnt.sim.model)
