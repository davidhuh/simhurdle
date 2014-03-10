### Example R code for using the simhurdle() function to zero-altered data
### from fitted hurdle model estimates

library(devtools)
library(glmmADMB)

## Load the simhurdle() function from github 
source_url("https://raw.github.com/davidhuh/simhurdle/master/simhurdle.R")


### Create artificial covariate data based on the DASH study

## Matrix of predictors
n.clust <- 30
n.id    <- 172
N.obs   <- n.id * n.clust

id      <- rep(1:n.id, each=n.clust)
day     <- rep(0:(n.clust-1), times=n.id)
dayofwk <- rep(rep(0:6, length.out=30), length.out=N.obs)

## randomly sample covariate values
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


#### Simulate outcome data

## Load fitted model estimates from which the base the simulation

beta.bin.names <- c("(Intercept)", "dmqsoc", "weekdayMon", "weekdayTue", "weekdayWed",
                    "weekdayThu", "weekdayFri", "weekdaySat", "dmqsoc:weekdayMon", "dmqsoc:weekdayTue",
                    "dmqsoc:weekdayWed", "dmqsoc:weekdayThu", "dmqsoc:weekdayFri", "dmqsoc:weekdaySat")

beta.bin <- c(-1.63053538, -0.15014852, -0.19634891, 0.26110985, -0.03715277, 0.68183510, 1.86443855, 1.65213433, -0.13456237, 0.15731037, 0.03097067, 0.20967669, 0.24922289, 0.13883190)
names(beta.bin) <- beta.bin.names

vcov.beta.bin <- matrix(c( 0.0216354681,  2.534447e-04,  1.809076e-03,  6.294931e-04,  1.497850e-03,  2.365040e-04, -1.357588e-03, -1.830763e-03,  8.557039e-04, -1.118421e-04,  0.0002916545, -1.814282e-04, -1.638488e-04, -0.0003221624,
                           0.0002534447,  2.168845e-02,  7.575438e-04, -1.414371e-04,  1.884152e-04, -1.739706e-04, -9.521889e-06, -2.779456e-04,  1.453559e-03,  5.275453e-04,  0.0011299571,  1.259441e-04, -1.188763e-03, -0.0014670354,
                           0.0018090759,  7.575438e-04,  3.168044e-02, -8.996861e-05, -6.369947e-04,  2.540645e-04,  3.423666e-04,  6.776400e-04,  4.837519e-03,  8.721741e-05, -0.0001258250,  1.083073e-04,  1.724082e-05,  0.0001790380,
                           0.0006294931, -1.414371e-04, -8.996861e-05,  2.838888e-02, -1.028297e-03,  2.211532e-04,  6.263986e-04,  1.167808e-03,  1.304910e-04,  3.046264e-04, -0.0002062624,  1.995110e-04,  1.332850e-04,  0.0002706246,
                           0.0014978500,  1.884152e-04, -6.369947e-04, -1.028297e-03,  2.989095e-02,  3.744320e-04,  6.120161e-04,  1.158923e-03, -8.096197e-05, -2.307845e-04,  0.0018303152,  1.620709e-04,  3.070250e-05,  0.0002636669,
                           0.0002365040, -1.739706e-04,  2.540645e-04,  2.211532e-04,  3.744320e-04,  2.691896e-02,  1.238497e-03,  2.098163e-03,  1.595735e-04,  2.162387e-04,  0.0002178249,  1.079318e-05,  2.489821e-04,  0.0005164031,
                          -0.0013575883, -9.521889e-06,  3.423666e-04,  6.263986e-04,  6.120161e-04,  1.238497e-03,  2.612749e-02,  3.984272e-03,  2.911298e-05,  1.775300e-04,  0.0001003319,  3.136831e-04,  5.010271e-04,  0.0007447751,
                          -0.0018307630, -2.779456e-04,  6.776400e-04,  1.167808e-03,  1.158923e-03,  2.098163e-03,  3.984272e-03,  2.647129e-02,  7.912052e-05,  2.721650e-04,  0.0002103955,  5.538832e-04,  8.457745e-04,  0.0012142366,
                           0.0008557039,  1.453559e-03,  4.837519e-03,  1.304910e-04, -8.096197e-05,  1.595735e-04,  2.911298e-05,  7.912052e-05,  3.243961e-02, -2.799853e-04, -0.0008105242,  7.997640e-05,  5.757238e-04,  0.0009292293,
                          -0.0001118421,  5.275453e-04,  8.721741e-05,  3.046264e-04, -2.307845e-04,  2.162387e-04,  1.775300e-04,  2.721650e-04, -2.799853e-04,  2.855086e-02, -0.0011158248,  8.336642e-05,  7.556149e-04,  0.0013076331,
                           0.0002916545,  1.129957e-03, -1.258250e-04, -2.062624e-04,  1.830315e-03,  2.178249e-04,  1.003319e-04,  2.103955e-04, -8.105242e-04, -1.115825e-03,  0.0297286564,  1.673015e-04,  9.018697e-04,  0.0014126578,
                          -0.0001814282,  1.259441e-04,  1.083073e-04,  1.995110e-04,  1.620709e-04,  1.079318e-05,  3.136831e-04,  5.538832e-04,  7.997640e-05,  8.336642e-05,  0.0001673015,  2.704709e-02,  1.385932e-03,  0.0022012642,
                          -0.0001638488, -1.188763e-03,  1.724082e-05,  1.332850e-04,  3.070250e-05,  2.489821e-04,  5.010271e-04,  8.457745e-04,  5.757238e-04,  7.556149e-04,  0.0009018697,  1.385932e-03,  2.606287e-02,  0.0036799796,
                          -0.0003221624, -1.467035e-03,  1.790380e-04,  2.706246e-04,  2.636669e-04,  5.164031e-04,  7.447751e-04,  1.214237e-03,  9.292293e-04,  1.307633e-03,  0.0014126578,  2.201264e-03,  3.679980e-03,  0.0263218176),
                        nrow=14,
                        byrow=TRUE,
                        dimnames=list(beta.bin.names, beta.bin.names))

beta.cnt.names <- c("(Intercept)", "dmqsoc", "cycPha", "cycAmp", "dmqsoc:cycPha", "dmqsoc:cycAmp")

beta.cnt <- c(1.152691622, 0.203528654, -0.206173030, -0.206259375, 0.113605291, -0.003125211)
names(beta.cnt) <- beta.cnt.names

vcov.beta.cnt <- matrix(c( 2.065793e-03, -1.919220e-04,  1.603239e-04,  9.402623e-05, -8.099336e-05,  1.858301e-05,
                          -1.919220e-04,  2.052724e-03, -8.969975e-05,  0.000000e+00,  1.864548e-04,  1.006906e-04,
                           1.603239e-04, -8.969975e-05,  1.489420e-03,  9.377653e-05, -2.436334e-04, -6.249133e-07,
                           9.402623e-05,  0.000000e+00,  9.377653e-05,  1.540955e-03,  9.881131e-06, -1.914843e-04,
                          -8.099336e-05,  1.864548e-04, -2.436334e-04,  9.881131e-06,  1.596402e-03,  1.159689e-04,
                           1.858301e-05,  1.006906e-04, -6.249133e-07, -1.914843e-04,  1.159689e-04,  1.638711e-03),
                        nrow=6,
                        byrow=TRUE,
                        dimnames=list(beta.cnt.names, beta.cnt.names))

re.bin.names <- "(Intercept)"
vcov.id.bin <- matrix(0.94573, dimnames=list(re.bin.names, re.bin.names))

re.cnt.names <- c("(Intercept)","cycPha","cycAmp")
vcov.id.cnt <- matrix(c(0.178290, 0.000000, 0.000000,
                        0.000000, 0.019143, 0.000000,
                        0.000000, 0.000000, 0.041581),
                      nrow=3,
                      byrow=TRUE,
                      dimnames=list(re.cnt.names, re.cnt.names))

alpha.mu.cnt <- 6.7429


## Simulate outcome values based based on artifical covariate values and fitted model estimates
y.sim <- simhurdle(data=pred.df,
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

## Combine artifical covariate data and simulated outcome values
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

##### Fit model in glmmADMB to simulated data

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
