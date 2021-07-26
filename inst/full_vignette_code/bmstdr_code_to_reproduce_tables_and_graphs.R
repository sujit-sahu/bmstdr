## The code in this file will re-produce the table results reported in the vignette 
## and paper: 
## bmstdr: Bayesian Modeling of Spatio-Temporal Data with R
## by Sujit Sahu, University of Southampton, 
## Email: S.K.Sahu@soton.ac.uk
## Binary versions of the package are available from: 
## https://www.southampton.ac.uk/~sks/bmbook/
## A second package databmstdr is also required to reproduce 
## some of the results. Please download and install it from the above webpage.

## The total time to run all the code in this vignette is about 2 hours on a fast PC. 
## Please start with a clear work space and restart R
rm(list=ls())
## Set a suitable working directory
yourpath <- "~/Dropbox/sks/bookproject/rbook/bmstdr/inst/full_vignette_code" 
setwd(yourpath)

## Please set a folder path for saving the tables. 
## The table numbers and file names refer to the tables in the 
## paper submitted to the JSS and also the vignette distributed 
## with the package 
tablepath <- "../txttables/"

## Please set a folder path for saving the figures. 
## The figure number and file name refers to the figures in the 
## paper submitted to the JSS only.
## The figures in the package vignette are drawn by the vignette Rmd file itself. 
figurepath <- "figures/"
# Set a file name to save the screen output if required 
# sink("sahu_jss-bmstdr.txt")
## Please load all the libraries. You may have to install these as required. 
library("bmstdr")
library("ggplot2")
require("databmstdr")
require("ggsn")
library("spTimer")
library("spTDyn")
library("tidyr")
library("xtable")
library("ggpubr")
library("akima")
library("tidyr")
library("doBy")


## Note the start time
start.time<-proc.time()[3]

## Code for drawing Figure 1. 
## This is a map of the state of New York where we also plot the 
## 28 air pollution monitoring sites. 


N <- 5000 # This is also the package default for number of iterations
Nstan <- 2000  # Stan runs take bit longer to run but requires 
## less sample size and hence we work with these choices 
burn.in <- 1000 # This is also the package default 

nymap <- map_data(database="state",regions="new york")
s <- c(1, 5, 10)
fcoords <- nyspatial[-s, c("Longitude", "Latitude")]
vcoords <- nyspatial[s,  c("Longitude", "Latitude")]

label <- tibble(
   long = -76.1,
   lat = 41.5,
   label = "25 fitted (circles) & 3  \n  validation (numbered) sites"
)
vsites3 <- ggplot() +
   geom_polygon(data=nymap, aes(x=long, y=lat, group=group),
                color="black", size = 0.6, fill=NA) +
   geom_point(data =fcoords, aes(x=Longitude,y=Latitude)) +
   geom_text(data=vcoords, aes(x=Longitude,y=Latitude, label=s), col=4) +
   labs(x="Longitude", y = "Latitude") +
   # geom_text(aes(label=label, x=long, y=lat), data = label, vjust = "top", hjust = "right")  +
   # geom_rect(mapping=aes(xmin=-78.7, xmax=-75.8, ymin=41, ymax=41.6), color="black", fill=NA) + 
   ggsn::scalebar(data =nymap, dist = 100, location = "bottomleft", transform=T, dist_unit = "km",
                  st.dist = .05, st.size = 5, height = .06, st.bottom=T, model="WGS84") +
   ggsn::north(data=nymap, location="topleft", symbol=12) 


## Draw the figure: fig.height=3, fig.width=4, fig.cap="Air pollution monitoring sites in New York."----
# 
vsites3 
## We save the figure. 
ggsave(filename = paste0(figurepath, "figure1.png"))

## Fitting independent error regression model 
f1 <- yo3 ~ xmaxtemp + xwdsp + xrh
M1 <- Bspatial(formula=f1, data=nyspatial, mchoice=T, N=N, burn.in = burn.in)
plot(M1)
## Plots are not saved 
summary(M1)
##  

## Fitting a spatial model without the nugget effect. 
## This run is very quick. 
M2 <- Bspatial(model="spat", formula=f1, data=nyspatial, 
        coordtype="utm", coords=4:5, phi=0.4, mchoice=T, N=N, burn.in = burn.in)
## 
summary(M2)
plot(M2)

## Illustrating a grid search method for choosing the decay parameter. 
# ?phichoice_sp
a <- phichoice_sp()
a

## Model fitting using spBayes 
M3 <- Bspatial(package="spBayes", formula=f1, data=nyspatial, coordtype="utm", 
            coords=4:5, prior.phi=c(0.005, 2), mchoice=T, 
            N=N, burn.in = burn.in, n.report=2)
summary(M3)

## Model fitting using code written in STAN
## Run with 2000 iterations and 1000 burn-in
## This takes 16 minutes to run. 
M4 <- Bspatial(package="stan", formula=f1, data=nyspatial, coordtype="utm", 
          coords=4:5,phi=0.4, mchoice=T, N=Nstan)
summary(M4)



## Model fitting using the INLA package 
M5  <- Bspatial(package="inla",formula=f1, data=nyspatial, 
            coordtype="utm", coords=4:5, mchoice=T, 
            N=N-burn.in)
summary(M5)

## The following command gets the model choice the results for fitting a baseline model with 
## without any covariates, i.e. with the intercept only.  
a3 <- Bmchoice(case="MC.sigma2.unknown", y=ydata)
## Now organize the all the results for forming Table 1. 
a5 <- rep(NA, 11)
a5[c(1, 3, 5, 7, 9:11)] <- unlist(M5$mchoice)
table1 <- cbind.data.frame(unlist(a3), M1$mchoice, M2$mchoice, M3$mchoice, M4$mchoice,  a5)
colnames(table1) <- paste("M", 0:5, sep="")
round(table1,  2)

dput(table1, file=paste0(tablepath, "table1.txt"))
## End reproducing Table 1. 

## Commands to reproduce model validation statistics in Table 2. 
## Select eight sites for validation 
s <- c(8,11,12,14,18,21,24,28)
M1.v <-  Bspatial(formula=f1, data=nyspatial, validrows=s, N=N, burn.in = burn.in )
M2.v <- Bspatial(package="none", model="spat", formula=f1, data=nyspatial,   
        coordtype="utm", coords=4:5,phi=0.4,  validrows=s, N=N, burn.in = burn.in )
M3.v <-  Bspatial(package="spBayes", prior.phi=c(0.005, 2), 
      formula=f1, data=nyspatial,   coordtype="utm", coords=4:5, 
      validrows=s, N=N, burn.in = burn.in)
# Takes 15 minutes to run 
M4.v  <- Bspatial(package="stan",formula=f1, data=nyspatial,   
      coordtype="utm", coords=4:5,phi=0.4 , validrows=s, N=Nstan)
M5.v  <- Bspatial(package="inla", formula=f1, data=nyspatial, 
        coordtype="utm", coords=4:5, validrows=s, N=N, burn.in = burn.in)

summary(M4.v)
## Gather results for Table 2
table2 <- cbind(unlist(M1.v$stats), unlist(M2.v$stats), unlist(M3.v$stats), 
          unlist(M4.v$stats), unlist(M5.v$stats))
colnames(table2) <- paste("M", 1:5, sep="")
round(table2, 3)
dput(table2, file=paste0(tablepath, "table2.txt"))

# Remove the last row as values are 100

## End reproducing Table 2. 


## Illustrate K-fold validation using M2 only. The runs are very fast. 
## Code below will reproduce Figure 2 and Table 3. 


set.seed(44)
x <- runif(n=28)
u <- order(x)
s1 <- u[1:7]
s2 <- u[8:14]
s3 <- u[15:21]
s4 <- u[22:28]

summary((1:28) - sort(c(s1, s2, s3, s4))) ## check

M2.v1 <- Bspatial(model="spat", formula=yo3~xmaxtemp+xwdsp+xrh, data=nyspatial, 
               coordtype="utm", coords=4:5,validrows= s1, phi=0.4, mchoice=F, 
               N=N, burn.in = burn.in)
M2.v2 <- Bspatial(model="spat", formula=yo3~xmaxtemp+xwdsp+xrh, data=nyspatial, 
               coordtype="utm", coords=4:5,validrows= s2,  phi=0.4, mchoice=F, 
               N=N, burn.in = burn.in)

M2.v4 <- Bspatial(model="spat", formula=yo3~xmaxtemp+xwdsp+xrh, data=nyspatial, 
               coordtype="utm", coords=4:5,validrows= s4,  phi=0.4, mchoice=F, 
               N=N, burn.in = burn.in)

M2.v3 <- Bspatial(model="spat", formula=yo3~xmaxtemp+xwdsp+xrh, data=nyspatial, 
                  coordtype="utm", coords=4:5, validrows= s3, phi=0.4, mchoice=F, 
                  N=N, burn.in = burn.in )
## Save this as Figure 2.  
ggsave(filename = paste0(figurepath, "figure2.png"))

## Now gather results to form Table 3. 

table3 <- cbind(unlist(M2.v1$stats), unlist(M2.v2$stats), unlist(M2.v3$stats), unlist(M2.v4$stats))
dimnames(table3)[[2]] <- paste("Fold", 1:4, sep="")
round(table3, 3)
dput(table3, file=paste0(tablepath, "table3.txt"))

## Finish reproducing Table 3. 



## Section 3:  Code for Modeling point reference spatio-temporal data 

f2 <- y8hrmax ~ xmaxtemp+xwdsp+xrh
M1 <- Bsptime(model="lm", formula=f2, data=nysptime, scale.transform = "SQRT",
              N=Ncar, burn.in = burn.in)
M2 <- Bsptime(model="separable", formula=f2, data=nysptime, 
              scale.transform = "SQRT", coordtype="utm", coords=4:5,  
              N=N, burn.in = burn.in)

names(M2)
summary(M2)
M2$phi.s
M2$phi.t

## Code for producing Figure 3. 
a <- residuals(M2)

ggsave(filename = paste0(figurepath, "figure3.png"))

valids <-  c(8,11,12,14,18,21,24,28)
vrows <-  which(nysptime$s.index%in% valids)
a2 <- residuals(M2)
a1 <- residuals(M1)
a1 <- residuals(M1, numbers=list(sn=28, tn=62))

## Demonstrate grid search method for selecting the spatial
## temporal decay parameters. 
# ?phichoicep
# Takes few minutes to run for 20 values 
a <- phichoicep()

valids <-  c(8,11,12,14,18,21,24,28)
vrows <-  which(nysptime$s.index%in% valids)
M2.1 <- Bsptime(model="separable",  formula=f2, data=nysptime, 
        validrows=vrows, coordtype="utm", coords=4:5,
        phi.s=0.005, phi.t=0.05, scale.transform = "SQRT", 
        N=N, burn.in = burn.in)

summary(M2.1)
plot(M2.1, segments = F)

M3 <- Bsptime(package="spTimer", formula=f2, data=nysptime, model="GP", 
              coordtype="utm", coords=4:5, scale.transform = "SQRT", 
              N=N, burn.in = burn.in, n.report=5)

# plot(M3)
print(M3)
summary(M3)
a3 <- residuals(M3)

tn <- 62; sn <- 28; valids <- c(1, 5, 10)
validt <- sort(sample(1:tn, size=31))
vrows <- getvalidrows(sn=sn, tn=tn, valids=valids, validt=validt)

M31 <- Bsptime(package="spTimer",formula=f2, data=nysptime, 
    coordtype="utm", coords=4:5, validrows=vrows, model="GP", 
    mchoice=F, scale.transform = "NONE", 
    N=N, burn.in = burn.in, n.report=5)

modfit <- M31$fit
fitall <- data.frame(modfit$fitted)
fitall$s.index <- rep(1:sn, each=tn)
vdat <- spT.subset(data=nysptime, var.name=c("s.index"), s=valids)
fitvalid <- spT.subset(data=fitall, var.name=c("s.index"), s=valids)
fitvalid$low <- fitvalid$Mean - 1.96 * fitvalid$SD
fitvalid$up <- fitvalid$Mean + 1.96 * fitvalid$SD
fitvalid$yobs <- vdat$y8hrmax
yobs <- matrix(fitvalid$yobs, byrow=T, ncol=tn)
y.valids.low <- matrix(fitvalid$low, byrow=T, ncol=tn)
y.valids.med <- matrix(fitvalid$Mean, byrow=T, ncol=tn)
y.valids.up <- matrix(fitvalid$up, byrow=T, ncol=tn)

p1 <- fig11.13.plot(yobs[1, ], y.valids.low[1, ], y.valids.med[1, ], 
                    y.valids.up[1, ], misst=validt)
p1 <- p1 + ggtitle("Validation for Site 1")

p1 
ggsave(filename = paste0(figurepath, "figure4.png"))

p2 <- fig11.13.plot(yobs[2, ], y.valids.low[2, ], y.valids.med[2, ],
                     y.valids.up[2, ], misst=validt)
p2 <- p2 + ggtitle("Validation for Site 5")
p2
p3 <- fig11.13.plot(yobs[3, ], y.valids.low[3, ], y.valids.med[3, ],
                     y.valids.up[3, ], misst=validt)
p3 <- p3 + ggtitle("Validation for Site 10")
p3


ggarrange(p1, p2, p3, common.legend = TRUE, legend = "top", nrow = 3, ncol = 1)
# This plot has not been  included in the JSS paper 


# function to calculate site-wise averages 
sitemeans <- function(a, sn, tn=62) { 
   u <- matrix(a, nrow=sn, ncol=tn, byrow=T)
   as.vector(apply(u, 1, mean))
}


post <- M3$fit
gpred <- predict(post, newdata=gridnysptime, newcoords=~Longitude+Latitude)
u <- gpred$pred.samples
v <- apply(u, 2, sitemeans, sn=100)
a <- get_parameter_estimates(t(v)) 
b <- data.frame(gridnyspatial[, 1:5], a) 


## Extracting the fitted values 
meanmat <- post$op
sig2eps <-  post$sig2ep
sige <- sqrt(sig2eps)
itmax <- ncol(meanmat)
nT <- nrow(nysptime)
sigemat <- matrix(rep(sige, each=nT), byrow=F, ncol=itmax)
a <- matrix(rnorm(nT*itmax), nrow=nT, ncol=itmax)
ypreds <- meanmat + a * sigemat
ypreds <-  (ypreds)^2
v <- apply(ypreds, 2, sitemeans, sn=28)
a <- get_parameter_estimates(t(v)) 
fits <- data.frame(nyspatial[, 1:5], a)

# Combine the fitted values and the predictions 
b <- rbind(b, fits)


coord <- nyspatial[, c("Longitude","Latitude")]
xo <- seq(from=min(coord$Longitude)-0.5, to = max(coord$Longitude)+0.8, length=200)
yo <- seq(from=min(coord$Latitude)-0.25, to = max(coord$Latitude)+0.8, length=200)
surf <- interp(b$Longitude, b$Latitude, b$mean,  xo=xo, yo=yo)
v <- fnc.delete.map.XYZ(xyz=surf)
interp1 <- data.frame(long = v$x, v$z )
names(interp1)[1:length(v$y)+1] <- v$y
interp1 <- gather(interp1,key = lat,value =Predicted,-long,convert = TRUE)

nymap <- map_data(database="state",regions="new york")
mappath <- cbind(nymap$long, nymap$lat)
zr <- range(interp1$Predicted, na.rm=T)
P <- ggplot() +  
geom_raster(data=interp1, aes(x = long, y = lat,fill = Predicted)) +
geom_polygon(data=nymap, aes(x=long, y=lat, group=group), color="black", size = 0.6, fill=NA) + 
geom_point(data=coord, aes(x=Longitude,y=Latitude))  +
stat_contour(data=na.omit(interp1), aes(x = long, y = lat,z = Predicted), colour = "black", binwidth =2) +
scale_fill_gradientn(colours=colpalette, na.value="gray95", limits=zr) +
theme(axis.text = element_blank(), axis.ticks = element_blank()) +
ggsn::scalebar(data =interp1, dist = 100, location = "bottomleft", transform=T, dist_unit = "km", st.dist = .05, st.size = 5, height = .06, st.bottom=T, model="WGS84") +
ggsn::north(data=interp1, location="topleft", symbol=12) +
labs(x="Longitude", y = "Latitude", size=2.5) 

# Repeat for sd 
surf <- interp(b$Longitude, b$Latitude, b$sd,  xo=xo, yo=yo)
v <- fnc.delete.map.XYZ(xyz=surf)
interp1 <- data.frame(long = v$x, v$z )
names(interp1)[1:length(v$y)+1] <- v$y
interp1 <- gather(interp1,key = lat,value =sd,-long,convert = TRUE)
nymap <- map_data(database="state",regions="new york")
mappath <- cbind(nymap$long, nymap$lat)
zr <- range(interp1$sd, na.rm=T)
Psd <- ggplot() +
    geom_raster(data=interp1, aes(x = long, y = lat,fill = sd)) +
    geom_polygon(data=nymap, aes(x=long, y=lat, group=group), color="black", size = 0.6, fill=NA) +
    geom_point(data=coord, aes(x=Longitude,y=Latitude))  +
    stat_contour(data=na.omit(interp1), aes(x = long, y = lat,z = sd), colour = "black", binwidth =0.1) +    scale_fill_gradientn(colours=colpalette, na.value="gray95", limits=zr) +
    theme(axis.text = element_blank(), axis.ticks = element_blank()) +
    ggsn::scalebar(data =interp1, dist = 100, location = "bottomleft", transform=T, dist_unit = "km",
                   st.dist = .05, st.size = 5, height = .06, st.bottom=T, model="WGS84") +
    ggsn::north(data=interp1, location="topleft", symbol=12) +
    labs(x="Longitude", y = "Latitude", size=2.5)

# Psd
ggpubr::ggarrange(P, Psd, common.legend = FALSE, nrow = 1, ncol = 2)
ggsave(filename = paste0(figurepath, "figure5.png"))

# Takes 2 minutes to run 
M4 <- Bsptime(package="stan",formula=f2, data=nysptime,   coordtype="utm", coords=4:5, 
              scale.transform = "SQRT", mchoice=T, verbose = F, N=Nstan)

M1 <- Bsptime(model="lm", formula=f2, data=nysptime, scale.transform = "SQRT", 
            mchoice=T, N=N, burn.in = burn.in)
M2 <- Bsptime(model="separable", formula=f2, data=nysptime, 
    scale.transform = "SQRT", coordtype="utm", coords=4:5, mchoice=T, 
    N=N, burn.in = burn.in)
M3 <- Bsptime(package="spTimer", formula=f2, data=nysptime, 
              model="GP", coordtype="utm", coords=4:5, 
        scale.transform = "SQRT", mchoice=T, N=N, burn.in = burn.in, 
        n.report=5)

table4 <- cbind(M1$mchoice, M2$mchoice, M3$mchoice, M4$mchoice)
round(table4, 2)
dput(table4, file=paste0(tablepath, "table4.txt"))


M5 <- Bsptime(package="spTimer", model="AR", formula=f2, data=nysptime, 
                coordtype="utm", coords=4:5, scale.transform = "SQRT", 
                mchoice=T,  N=N, burn.in = burn.in, n.report=5)

a <- residuals(M5)
# Takes 2 mins 23 sec 
M6 <- Bsptime(package="inla", model="AR", formula=f2, data=nysptime,
         coordtype="utm", coords=4:5, scale.transform = "SQRT",
         mchoice=T,  N=N-burn.in)

summary(M5)
summary(M6)

table5 <- rbind(unlist(M5$mchoice), unlist(M6$mchoice[5:7]))
rownames(table5) <- c("spTimer", "INLA")
table5
dput(table5, file=paste0(tablepath, "table5.txt"))

M5$params
M6$params

b1 <- rbind(rep(NA, 4), M6$params)
table6 <- cbind.data.frame(M5$params, b1)

round(table6, 3)
dput(table6, file=paste0(tablepath, "table6.txt"))


## Model fitting using the sptDyn package

f3 <- y8hrmax~ xmaxtemp + sp(xmaxtemp)+ tp(xwdsp) + xrh
M7 <- Bsptime(package="sptDyn", model="GP", formula=f3, data=nysptime, 
      coordtype="utm", coords=4:5, scale.transform = "SQRT", 
      mchoice=T,  n.report=5, N=N, burn.in = burn.in )

## Generating the plots 

out <- M7$fit
a <- out$betasp
u <- c(t(out$betasp))
sn <- nrow(a)
itmax <- ncol(a)
v <- rep(1:sn, each=itmax)
d <- data.frame(site=as.factor(v), sp = u)

ptmp <- ggplot(data=d, aes(x=site, y=sp)) + 
   geom_boxplot(outlier.colour="black", outlier.shape=1,
                outlier.size=0.5) +
   geom_abline(intercept=0, slope=0, color="blue") + 
   labs(title= "Spatial effects of maximum temperature", x="Site", y = "Effects", size=2.5) 


b <- out$betatp
tn <- nrow(b)
itmax <- ncol(b)
tids <- 1:tn 
stat <- apply(b[tids,], 1, quantile, prob=c(0.025,0.5,0.975))
tstat <- data.frame(tids, t(stat))
dimnames(tstat)[[2]] <- c("Days", "low", "median", "up")
# head(tstat)
yr <- c(min(c(stat)),max(c(stat)))
pwdsp <- ggplot(data=tstat, aes(x=Days, y=median)) + 
   geom_point(size=3) + 
   ylim(yr) + 
   geom_segment(data=tstat, aes(x=Days, y=median, xend=Days, yend=low), linetype=1) +
   geom_segment(data=tstat, aes(x=Days, y=median, xend=Days, yend=up), linetype=1) +
   geom_abline(intercept=0, slope=0, col="blue") +
   labs(title="Temporal effects of wind speed", x="Days", y="Temporal effects") 


library(ggpubr)
ggarrange(ptmp, pwdsp, common.legend = FALSE, nrow = 2, ncol = 1)
ggsave(filename = paste0(figurepath, "figure6.png"))

## Model fitting using spBayes 
M8 <- Bsptime(package="spBayes",  formula=f2, data=nysptime, 
prior.sigma2=c(2, 25), prior.tau2 =c(2, 25), prior.sigma.eta =c(2, 0.001),
coordtype="utm", coords=4:5, scale.transform = "SQRT", 
mchoice=T,  n.report=5, N=N,  burn.in = burn.in)


modfit <- M8$fit
tn <- 62
quant95 <- function(x){
   quantile(x, prob=c(0.025, 0.5, 0.975))
}
beta <- apply(modfit$p.beta.samples[burn.in:N,], 2, quant95)
theta <- apply(modfit$p.theta.samples[burn.in:N,], 2, quant95)
sigma.sq <- theta[,grep("sigma.sq", colnames(theta))]
tau.sq <- theta[,grep("tau.sq", colnames(theta))]
phi <- theta[,grep("phi", colnames(theta))]

adat <- data.frame(x=1:tn, low=sigma.sq[1, ], med=sigma.sq[2, ], up=sigma.sq[3, ])
head(adat)

psigma <- ggplot() + 
   geom_point(data = adat, aes(x =x, y = med, shape=19), shape=19, col="blue", size = 2) + 
   geom_ribbon(data = adat, aes(x =x, ymin =low, ymax = up), alpha = 0.2, color = "grey50") +
   theme(legend.position ="none") + 
   labs(y = "sigma2", x = "Days") 
psigma

adat <- data.frame(x=1:tn, low=tau.sq[1, ], med=tau.sq[2, ], up=tau.sq[3, ])

ptau <- ggplot() + 
   geom_point(data = adat, aes(x =x, y = med, shape=19), shape=19, col="blue", size = 2) + 
   geom_ribbon(data = adat, aes(x =x, ymin =low, ymax = up), alpha = 0.2, color = "grey50") +
   theme(legend.position ="none") + 
   labs(y = "tau2", x = "Days") 
ptau

adat <- data.frame(x=1:tn, low=3/phi[3, ], med=3/phi[2, ], up=3/phi[1, ])


prange <- ggplot() + 
   geom_point(data = adat, aes(x =x, y = med, shape=19), shape=19, col="blue", size = 2) + 
   geom_ribbon(data = adat, aes(x =x, ymin =low, ymax = up), alpha = 0.2, color = "grey50") +
   theme(legend.position ="none") + 
   labs(y = "Range", x = "Days") 
prange


ggarrange(psigma, prange, common.legend = TRUE, legend = "none", nrow = 2, ncol = 1)
ggsave(filename = paste0(figurepath, "figure7.png"))

vnames <- all.vars(f2)
xnames <- vnames[-1]
k <- 4
cat("\nOnly first 4 beta parameters are plotted\n")
beta.0 <- beta[,grep("Intercept", colnames(beta))]
adat <- data.frame(x=1:tn, low=beta.0[1,], med=beta.0[2,], up=beta.0[3,])
head(adat)
 pint <- ggplot() +
    geom_point(data = adat, aes(x =x, y = med, shape=19), shape=19, col="blue", size = 2) +
    geom_ribbon(data = adat, aes(x =x, ymin =low, ymax = up), alpha = 0.2, color = "grey50") +
    geom_hline(yintercept=0, linetype="dashed", color = "red", size=1)+
    theme(legend.position ="none") +
    labs(y = "Intercept", x = "Days")
pint
j <- 2
betaj <- beta[,grep(xnames[j-1], colnames(beta))]
adat <- data.frame(x=1:tn, low=betaj[1,], med=betaj[2,], up=betaj[3,])
ptmp <- ggplot() +
    geom_point(data = adat, aes(x =x, y = med, shape=19), shape=19, col="blue", size = 2) +
    geom_ribbon(data = adat, aes(x =x, ymin =low, ymax = up), alpha = 0.2, color = "grey50") +
    geom_hline(yintercept=0, linetype="dashed", color = "red", size=1)+
    theme(legend.position ="none") +
   labs(y = "Max temp", x = "Days")
ptmp
j <- 3
betaj <- beta[,grep(xnames[j-1], colnames(beta))]
 adat <- data.frame(x=1:tn, low=betaj[1,], med=betaj[2,], up=betaj[3,])
 head(adat)
pwdsp <- ggplot() +
    geom_point(data = adat, aes(x =x, y = med, shape=19), shape=19, col="blue", size = 2) +
    geom_ribbon(data = adat, aes(x =x, ymin =low, ymax = up), alpha = 0.2, color = "grey50") +
    theme(legend.position ="none") +
    geom_hline(yintercept=0, linetype="dashed", color = "red", size=1) +
    labs(y = "Wind speed", x = "Days")
pwdsp
 j <- 4
 betaj <- beta[,grep(xnames[j-1], colnames(beta))]
 adat <- data.frame(x=1:tn, low=betaj[1,], med=betaj[2,], up=betaj[3,])
 head(adat)
 prh <- ggplot() +
    geom_point(data = adat, aes(x =x, y = med, shape=19), shape=19, col="blue", size = 2) +
    geom_ribbon(data = adat, aes(x =x, ymin =low, ymax = up), alpha = 0.2, color = "grey50") +
    theme(legend.position ="none") +
    geom_hline(yintercept=0, linetype="dashed", color = "red", size=1)+
    labs(y = "Rel humidity", x = "Days")
 prh
# Plot not included 
ggarrange(pint, ptmp, pwdsp, prh, common.legend = TRUE, legend = "none", nrow = 2, ncol = 2)
 
###
M9 <-  Bsptime(package="spTimer", model="GPP", g_size=5, formula=f2, 
data=nysptime, coordtype="utm", coords=4:5, scale.transform = "SQRT", 
N=N, burn.in = burn.in, n.report=5)
a <- residuals(M9)
summary(M9)

valids <- c(8,11,12,14,18,21,24,28)
vrows <- getvalidrows(sn=28, tn=62, valids=valids, allt=T)

M1.v <- Bsptime(model="lm", formula=f2, data=nysptime, 
                scale.transform = "SQRT", validrows=vrows, mchoice=T, 
                N=N, burn.in = burn.in)

M2.v <- Bsptime(model="separable",  formula=f2, data=nysptime, 
                coordtype="utm", coords=4:5, 
                phi.s=0.005, phi.t=0.05, scale.transform = "SQRT", 
                validrows=vrows, mchoice=T, 
                N=N, burn.in = burn.in)

M3.v <- Bsptime(package="spTimer", model="GP", formula=f2, data=nysptime, 
                coordtype="utm", coords=4:5, scale.transform = "SQRT", 
                mchoice=T, validrows=vrows, n.report=2, 
                N=N, burn.in = burn.in)

M4.v <- Bsptime(package="stan",formula=f2, data=nysptime, 
                coordtype="utm", coords=4:5, scale.transform = "SQRT", 
                 mchoice=T, validrows=vrows, verbose = F, 
                N=Nstan)

# 8 mins 
M5.v <- Bsptime(package="spTimer", model="AR", formula=f2, data=nysptime, 
                coordtype="utm", coords=4:5, scale.transform = "SQRT", 
                validrows=vrows, mchoice=T,  n.report = 2, 
                N=N, burn.in = burn.in)

M6.v <- Bsptime(package="inla", model="AR", formula=f2, data=nysptime, 
                coordtype="utm", coords=4:5, scale.transform = "SQRT", 
                validrows=vrows,  mchoice=T, N=N-burn.in)
#6 mins 41s


M7.v <- Bsptime(package="sptDyn", model="GP", formula=f3, data=nysptime, 
                coordtype="utm", coords=4:5, scale.transform = "SQRT", 
                mchoice=T, validrows=vrows, n.report=5, 
                N=N, burn.in = burn.in)

M8.v <- Bsptime(package="spBayes",  formula=f2, data=nysptime, 
                coordtype="utm", coords=4:5, scale.transform = "SQRT", 
                prior.sigma2=c(2, 25),
                prior.tau2 =c(2, 25),
                prior.sigma.eta =c(2, 0.001),
                mchoice=T,  validrows=vrows,  n.report=2, 
                N=N, burn.in = burn.in)

M9.v <-   Bsptime(package="spTimer", model="GPP", g_size=5, 
                  formula=f2, data=nysptime, validrow=vrows, 
                  coordtype="utm", coords=4:5, scale.transform = "SQRT", mchoice=T, 
                  n.report=2, N=N, burn.in = burn.in)


results <- cbind.data.frame(lm=unlist(M1.v$stats), separable=unlist(M2.v$stats), 
                            spTimerGP=unlist(M3.v$stats), stan=unlist(M4.v$stats),inla=unlist(M6.v$stats),
                            spTimerAR=unlist(M5.v$stats), spTDyn=unlist(M7.v$stats), 
                            spBayes=unlist(M8.v$stats),  sptimerGPP=unlist(M9.v$stats))

mcres <- cbind.data.frame(lm=unlist(M1.v$mchoice)[9:11], separable=unlist(M2.v$mchoice)[9:11], 
                          spTimerGP=unlist(M3.v$mchoice)[9:11], stan=unlist(M4.v$mchoice)[9:11],
                          inla=unlist(M6.v$mchoice)[5:7],
                          spTimerAR=unlist(M5.v$mchoice), spTDyn=unlist(M7.v$mchoice), 
                          spBayes=unlist(M8.v$mchoice),  sptimerGPP=unlist(M9.v$mchoice))                   
results
round(results, 2)
mcres
xtable(results, digits=2)
table7 <- rbind.data.frame(results, mcres)
round(table7, 2)
table7 <- table7[, -8]
table7

dput(table7, file=paste0(tablepath, "table7.txt"))

## Ocean temperature 


deep <- argo_floats_atlantic_2003[argo_floats_atlantic_2003$depth==3, ]
deep$x2inter <- deep$xinter*deep$xinter
deep$month <- factor(deep$month)

deep$lat2 <- (deep$lat)^2
deep$sin1 <- round(sin(deep$time*2*pi/365), 4)
deep$cos1 <- round(cos(deep$time*2*pi/365), 4)
deep$sin2 <- round(sin(deep$time*4*pi/365), 4)
deep$cos2 <- round(cos(deep$time*4*pi/365), 4)

head(deep)
## scaling variables 
deep[, c( "xlat2", "xsin1", "xcos1", "xsin2", "xcos2")] <- 
  scale(deep[,c("lat2", "sin1", "cos1", "sin2", "cos2")])

# 15 minutes 
f2 <- temp ~ xlon + xlat + xlat2+ xinter + x2inter 
M2atl <- Bmoving_sptime(formula=f2, data = deep, coordtype="lonlat", coords = 1:2,
                     N=1100, burn.in=100, validrows =NULL, mchoice =T)
summary(M2atl)
plot(M2atl)
names(M2atl)

a <-residuals(M2atl)
a <- fitted(M2atl)
summary(a)


listofdraws <- rstan::extract(M2atl$fit)
names(listofdraws)
dim(listofdraws$xbmodel)
dim(listofdraws$bigS)
dat <- M2atl$datatostan
names(dat)
tempdata <- dat



## Extract the diagonal elements of all St 
v <- numeric()
x <- numeric()
m <- length(dat$nts)
for (i in 1:m) {
  k <- dat$nts[i]
  a1 <- 1:k^2
  a2 <- seq(from=1, to =k^2, length=k)
  b1 <- rep(0, k^2)
  b1[which(a1==a2)] <- 1
  cbind(a1, b1)
  v <- c(v, a1) 
  x <- c(x, b1) ## indicates if the corresponding index is a diagonal 
}



varsamps <- listofdraws$bigS[, which(x>0)]/365
dim(varsamps)
xbeta <- listofdraws$xbmodel
dim(xbeta)
sdsamps <- sqrt(varsamps)

ntot <- nrow(xbeta) * ncol(xbeta)

zsamp <- matrix(rnorm(ntot), nrow=nrow(xbeta), ncol=ncol(xbeta))
dim(zsamp)
zsamp <- xbeta + zsamp * sdsamps
ansummary <- get_parameter_estimates(zsamp)
head(ansummary)
summary(ansummary$mean)
summary(ansummary$sd)
summary(ansummary$low)
summary(ansummary$up)


pdata <- cbind(deep, ansummary)


atlmap <- map_data("world", xlim=c(-70, 10), ylim=c(15, 65))
head(atlmap)
atlmap <- atlmap[atlmap$long < 5, ]
atlmap <- atlmap[atlmap$long > -70, ]
atlmap <- atlmap[atlmap$lat < 65, ]
atlmap <- atlmap[atlmap$lat > 10, ]

xo <- seq(from=-70, to = 5,  length=200)
yo <- seq(from= 10, to = 65, length=200)

## Temperature 
surf <- interp(pdata$lon, pdata$lat, pdata$mean,  xo=xo, yo=yo)
names(surf)
surf <- list(x=surf$x, y=surf$y, z=surf$z)
interp1 <- data.frame(long = surf$x, surf$z )
names(interp1)[1:length(surf$y)+1] <- surf$y
interp1 <- gather(interp1,key = lat,value =Temp,-long,convert = TRUE)

head(interp1)
summary(interp1$Temp)
pcolors <- topo.colors(5)
P <- ggplot() +  
  geom_raster(data=interp1, aes(x = long, y = lat,fill = Temp)) + 
  # scale_fill_gradientn(colours=c("#0000FFFF","#FFFFFFFF","#FF0000FF")) + 
  scale_fill_gradientn(colours=pcolors) +  
  geom_contour(data=interp1, aes(x = long, y = lat,z = Temp)) + 
  geom_polygon(data=atlmap, aes(x=long, y=lat, group=group),
               color="black", size = 0.6, fill=NA) + 
  geom_point(data =deep, aes(x=lon, y=lat, colour=month), size=1) +
  labs( title= "Annual temperature in deep ocean in 2003", x="Longitude", y = "Latitude") +
  ggsn::scalebar(data =atlmap, dist = 1000, location = "bottomleft", transform=T, dist_unit = "km",
                 st.dist = .05, st.size =5, height = .05, st.bottom=T, model="WGS84") +
  ggsn::north(data=atlmap, location="topright", symbol=12) 

P 

ggsave(filename = paste0(figurepath, "temp_deep.png"),  width=8.27, height=5, dpi=300)

## sd of temperature 

surf <- interp(pdata$lon, pdata$lat, pdata$sd,  xo=xo, yo=yo)
names(surf)
surf <- list(x=surf$x, y=surf$y, z=surf$z)
interp1 <- data.frame(long = surf$x, surf$z )
names(interp1)[1:length(surf$y)+1] <- surf$y
interp1 <- gather(interp1,key = lat,value =sd,-long,convert = TRUE)

head(interp1)
summary(interp1$sd)

psd <- ggplot() +  
  geom_raster(data=interp1, aes(x = long, y = lat,fill = sd)) + 
  # scale_fill_gradientn(colours=c("#0000FFFF","#FFFFFFFF","#FF0000FF")) + 
  scale_fill_gradientn(colours=pcolors) + 
  geom_contour(data=interp1, aes(x = long, y = lat,z = sd)) + 
  geom_polygon(data=atlmap, aes(x=long, y=lat, group=group),
               color="black", size = 0.6, fill=NA) + 
  geom_point(data =deep, aes(x=lon, y=lat, colour=month), size=1) +
  labs( title= "sd of annual temperature in deep ocean in 2003", x="Longitude", y = "Latitude") +
  ggsn::scalebar(data = atlmap, dist = 1000, location = "bottomleft", transform=T, dist_unit = "km",
                 st.dist = .05, st.size =5, height = .05, st.bottom=T, model="WGS84") +
  ggsn::north(data=atlmap, location="topright", symbol=12) 

psd 
ggsave(filename = paste0(figurepath, "sd_temp_deep.png"),  width=8.27, height=5, dpi=300)


## Ocean temperature example complete. 

## Remove the large model output 
rm(list=ls(pattern="M"))

# Code for Section 4 results 

## fig.cap="Weekly Covid-19 death rate per 100,000"----------
engdeaths$covidrate <- 100000*engdeaths$covid/engdeaths$popn
ptime <- ggplot(data=engdeaths,  aes(x=factor(Weeknumber), y=covidrate)) +
  geom_boxplot() +
  labs(x = "Week", y = "Death rate per 100,000")  +
  stat_summary(fun=median, geom="line", aes(group=1, col="red")) +
  theme(legend.position = "none")
ptime


ggsave(filename = paste0(figurepath, "figure8.png"))

bdf <- merge(englamap, engtotals, by.x="id", by.y="mapid", all.y=T, all.x=F)
bdf$covidrate <- bdf$covid/bdf$popn*100000
plimits <- range(bdf$covidrate)
prate <-  ggplot(data=bdf, aes(x=long, y=lat, group = group, fill=covidrate)) +
   scale_fill_gradientn(colours=colpalette, na.value="black",limits=plimits)  +
   geom_polygon(colour='black',size=0.25) +
   geom_polygon(data=engregmap, aes(x=long, y=lat, group = group), fill=NA, colour='black',size=0.6)  +
   coord_equal() + guides(fill=guide_colorbar(title="Death rate")) +
   theme_bw()+theme(text=element_text(family="Times")) +
   labs(x="", y = "") +
   theme(axis.text.x = element_blank(), axis.text.y = element_blank(),axis.ticks = element_blank())   +
   theme(legend.position =c(0.2, 0.5))

prate 

ggsave(filename = paste0(figurepath, "figure9.png"))

Ncar <- 50000
burn.in.car <- 10000
thin <- 10
inlaN <- (Ncar - burn.in.car)/thin

nweek <- rep(1, nrow(engdeaths))
scol <- "spaceid"
tcol <-  "Weeknumber"

f10 <- highdeathsmr ~  jsa + log10(houseprice) + log(popdensity) 

M1st_linear <- Bcartime(formula=f10, data=engdeaths, scol=scol, tcol=tcol, trials=nweek, 
                        W=Weng, model="linear", family="binomial", package="CARBayesST", 
                        N=Ncar, burn.in=burn.in.car, thin=thin)


M1st_anova_nointer <- Bcartime(formula=f10, data=engdeaths, scol=scol, tcol=tcol, trials=nweek, 
                               W=Weng, model="anova", interaction=F, family="binomial", 
                               package="CARBayesST", N=Ncar, burn.in=burn.in.car, thin=thin)
summary(M1st_anova_nointer)


M1st_sepspat <- Bcartime(formula=f10, data=engdeaths, scol=scol, tcol=tcol, trials=nweek, 
                         W=Weng, model="sepspatial", family="binomial", 
                         package="CARBayesST", N=Ncar, burn.in=burn.in.car, thin=thin)
summary(M1st_sepspat)


M1st_ar <- Bcartime(formula=f10, data=engdeaths, scol=scol, tcol=tcol, trials=nweek, 
                    W=Weng, model="ar", AR=1, family="binomial", package="CARBayesST", 
                    N=Ncar, burn.in=burn.in.car, thin=thin)
summary(M1st_ar)

M1st_ar2 <- Bcartime(formula=f10, data=engdeaths, scol=scol, tcol=tcol, trials=nweek, 
                     W=Weng, model="ar", AR=2, family="binomial", package="CARBayesST", 
                     N=Ncar, burn.in=burn.in.car, thin=thin)
summary(M1st_ar2)


model <- c("bym", "ar1")
M1st_inla.bym <- Bcartime(data=engdeaths, formula=f10, W=Weng, 
        scol =scol, tcol=tcol,  model=model,   trials=nweek, 
        family="binomial", package="inla", N=inlaN) 
summary(M1st_inla.bym)
names(M1st_inla.bym)

a <- rbind(M1st_linear$mchoice,  M1st_anova_nointer$mchoice, 
           M1st_sepspat$mchoice, M1st_ar$mchoice, M1st_ar2$mchoice)
a
a <- a[, -(5:6)]
a <- a[, c(2, 1, 4, 3)]
a

b <- M1st_inla.bym$mchoice[1:4]
table8 <- rbind(a, b)

rownames(table8) <- c("Linear", "Anova", "Separable", "AR (1)", "AR (2)", "INLA-BYM")
colnames(table8) <- c("pDIC", "DIC", "pWAIC", "WAIC")
round(table8, 2)
dput(table8, file=paste0(tablepath, "table8.txt"))

## ---- echo=TRUE--------------------------------------------------------------------------------------------------
f20 <-  covid ~ offset(logEdeaths) + jsa + log10(houseprice) + log(popdensity) + n0

# 2 mins 51 sec
M2st_linear <- Bcartime(formula=f20, data=engdeaths, scol=scol, tcol=tcol,  
                        W=Weng, model="linear", family="poisson", package="CARBayesST", 
                        N=Ncar, burn.in=burn.in.car, thin=thin)
summary(M2st_linear)


M2st_anova <- Bcartime(formula=f20, data=engdeaths, scol=scol, tcol=tcol,  
                       W=Weng, model="anova", family="poisson", package="CARBayesST", 
                       N=Ncar, burn.in=burn.in.car, thin=thin)
summary(M2st_anova)


M2st_anova_nointer <- Bcartime(formula=f20, data=engdeaths, scol=scol, tcol=tcol,  
                               W=Weng, model="anova",interaction=F,  family="poisson", 
                               package="CARBayesST", N=Ncar, burn.in=burn.in.car, thin=thin)
summary(M2st_anova_nointer)


M2st_sepspat <- Bcartime(formula=f20, data=engdeaths, scol=scol, tcol=tcol, 
                         W=Weng, model="sepspatial",family="poisson", package="CARBayesST", 
                         N=Ncar, burn.in=burn.in.car, thin=thin)
summary(M2st_sepspat)

M2st_ar <- Bcartime(formula=f20, data=engdeaths, scol=scol, tcol=tcol,  
                    W=Weng, model="ar", family="poisson", package="CARBayesST", 
                    N=Ncar, burn.in=burn.in.car, thin=thin)
summary(M2st_ar)


M2st_ar2 <- Bcartime(formula=f20, data=engdeaths, scol=scol, tcol=tcol,  
                     W=Weng, model="ar", AR=2, family="poisson", package="CARBayesST", 
                     N=Ncar, burn.in=burn.in.car, thin=thin)
summary(M2st_ar2)

model <- c("bym", "ar1")
f2inla <-  covid ~  jsa + log10(houseprice) + log(popdensity) + n0 

M2stinla <- Bcartime(data=engdeaths, formula=f2inla, W=Weng, scol =scol, tcol=tcol,  
                     offsetcol="logEdeaths",  model=model,  link="log", family="poisson", 
                     package="inla", N=inlaN) 


a <- rbind(M2st_linear$mchoice, M2st_anova_nointer$mchoice,  M2st_anova$mchoice,
           M2st_sepspat$mchoice, M2st_ar$mchoice, M2st_ar2$mchoice)
a
a <- a[, -(5:6)]
a <- a[, c(2, 1, 4, 3)]
a

b <- M2stinla$mchoice[1:4]
table9 <- rbind(a, b)

rownames(table9) <- c("Linear", "Anova (without interaction)", "Anova (with interaction)", "Separable", "AR (1)", "AR (2)", "INLA-BYM")
colnames(table9) <- c("pDIC", "DIC", "pWAIC", "WAIC")
table9
dput(table9, file=paste0(tablepath, "table9.txt"))

## Code for Figure 10

allpost <- M2st_ar2$fit
names(allpost)

a <- allpost$samples$fitted ## nsamp by nla*nweek
class(a)
b <- as.matrix(a)
class(b)

# head(b)

dim(b)

calwkmeans <- function(x, wkpattern=engdeaths$Weeknumber) {
  u <- data.frame(x=100000*x/engdeaths$popn, wkpattern=wkpattern)
  # Transform to per 100,000
  v <- summaryBy(x~wkpattern, data=u)
  as.vector(v$x.mean)
}
u <- calwkmeans(b[1,])
length(u) 
summary(u)

a <- apply(b, 1, calwkmeans)
dim(a)

fits <- apply(a, 1, mean) 
# fits <- apply(a, 1, quantile, probs=c(0.025, 0.5, 0.975)) 
fits
lims <-  apply(a, 1, FUN=quantile, probs=c(0.025, 0.975))  
lims
head(lims)
dim(lims)

library(doBy)

head(engdeaths)
engdeaths$rcovid <- 100000*engdeaths$covid/engdeaths$popn 
wkmeans <- summaryBy(rcovid ~ Weeknumber, data =engdeaths, FUN =mean) # Weekly mean 
colnames(wkmeans)[2] <- "covid"
head(wkmeans)

# dataandfits <- cbind(wkmeans, fits=fits[2, ], lower=fits[1,],  upper=fits[3,])
dataandfits <- cbind(wkmeans, fits=fits, lower=lims[1,],  upper=lims[2,])
head(dataandfits)


k <- nrow(dataandfits)
k
adata <- data.frame(Weeknumber=rep(dataandfits$Weeknumber, 3), 
                    fits=c(rep("fits", k), rep("lower", k), rep("upper", k)),  
                    vfits=c(dataandfits$fits, dataandfits$lower, dataandfits$upper))
adata$fits <- as.factor(adata$fits)
adata$fits <- factor(adata$fits, levels=rev(levels(adata$fits)))
head(adata)
table(adata$fits)
adata <- adata[adata$fits !="fits", ]

adata <- data.frame(Weeknumber=rep(dataandfits$Weeknumber, 2), 
                    intervals=c(rep("lower", k), rep("upper", k)),  
                    vfits=c(dataandfits$lower, dataandfits$upper)) 
adata$intervals <- as.factor(adata$intervals)
adata$intervals <- factor(adata$intervals, levels=levels(adata$intervals))
dim(adata)
head(adata)
table(adata$fits)
pdata <- data.frame(Weeknumber=rep(dataandfits$Weeknumber, 2), 
                    fits=c(rep("fitted", k), rep("observed", k)),  
                    vfits=c(dataandfits$fits, dataandfits$covid))
pdata$fits <- as.factor(pdata$fits)
pdata$fits <- factor(pdata$fits, levels=levels(pdata$fits))
levels(adata$fits)
levels(pdata$fits)

head(pdata)
pwkfit <- ggplot() + 
  geom_line(data=adata, aes(x=factor(Weeknumber), y=vfits, group=intervals, color=intervals)) +
  geom_point(data=pdata, aes(x=factor(Weeknumber), y=vfits,  shape=fits)) + 
  labs(x ="Week number", y = "Average number of covid deaths per 100,000")+
  theme(legend.position=c(0.65, 0.5))
pwkfit


ggsave(filename = paste0(figurepath, "figure10.png"))

###

set.seed(5)
vs <- sample(nrow(engdeaths), 0.1*nrow(engdeaths))

M2st_anova.0 <- Bcartime(formula=f20, data=engdeaths, scol=scol, tcol=tcol,  
                         W=Weng, model="anova", family="poisson", package="CARBayesST", 
                         N=Ncar, burn.in=burn.in.car, thin=thin, validrows=vs)


M2st_ar.0 <- Bcartime(formula=f20, data=engdeaths, scol=scol, tcol=tcol,  
                      W=Weng, model="ar", AR=1,  family="poisson", package="CARBayesST", 
                      N=Ncar, burn.in=burn.in.car, thin=thin, 
                      validrows=vs, verbose=T)

M2st_ar2.0 <- Bcartime(formula=f20, data=engdeaths, scol=scol, tcol=tcol,  
                       W=Weng, model="ar", AR=2, family="poisson", package="CARBayesST", 
                       N=Ncar, burn.in=burn.in.car, thin=thin, 
                       validrows=vs, verbose=T)


M2stinla.0  <- Bcartime(data=engdeaths, formula=f2inla, W=Weng, scol =scol, tcol=tcol,  
                        offsetcol="logEdeaths",  model=model,  link="log", family="poisson", 
                        package="inla", validrow=vs, N=inlaN) 

table10 <- rbind(unlist(M2st_anova.0$stats), unlist(M2st_ar.0$stats), unlist(M2st_ar2.0$stats), 
                 unlist(M2stinla.0$stats)) 

table10
table10 <- as.data.frame(table10)
rownames(table10) <- c("Anova", "AR (1)", "AR (2)", "INLA-BYM")
dput(table10, file=paste0(tablepath, "table10.txt"))

summary(M2st_ar2.0)
yobspred <- M2st_ar2.0$yobs_preds
names(yobspred)
yobs <- yobspred$covid
predsums <- get_validation_summaries(t(M2st_ar2.0$valpreds))
dim(predsums)
b <- obs_v_pred_plot(yobs, predsums, segments=T) 
b
names(M2stinla.0)
inlayobspred <- M2stinla.0$yobs_preds
names(inlayobspred)
inlapredsums <- get_validation_summaries(t(M2stinla.0$valpreds))
dim(inlapredsums)
a <- obs_v_pred_plot(yobs, inlapredsums, segments=T) 
a

inlavalid <- a$pwithseg
ar2valid <- b$pwithseg

ggarrange(ar2valid, inlavalid, common.legend = TRUE, legend = "top", nrow = 1, ncol = 2)
ggsave(filename = paste0(figurepath, "figure11.png"))

f2s <-  covid ~ offset(logEdeaths) + jsa + log10(houseprice) + log(popdensity)
M2s <- Bcartime(formula=f2s, data=engtotals, scol=scol, W=Weng, model="bym",
   AR=2, family="poisson", N=Ncar, burn.in=burn.in.car, thin=thin)
summary(M2s)


colnames(engdeaths)

f3 <-  sqrt(no2) ~  jsa + log10(houseprice) + log(popdensity) 

scol <- "spaceid"
tcol <-  "Weeknumber"

M3st_linear <- Bcartime(formula=f3, data=engdeaths, scol=scol, tcol=tcol, 
                        W=Weng, model="linear", family="gaussian", package="CARBayesST", 
                        N=Ncar, burn.in=burn.in.car, thin=thin)
summary(M3st_linear)


M3st_anova <- Bcartime(formula=f3, data=engdeaths, scol=scol, tcol=tcol, 
                       W=Weng, model="anova", family="gaussian", package="CARBayesST", 
                       interaction=F, N=Ncar, burn.in=burn.in.car, thin=thin)
summary(M3st_anova)



M3st_ar <- Bcartime(formula=f3, data=engdeaths, scol=scol, tcol=tcol, 
                    W=Weng, model="ar", AR=1, family="gaussian", package="CARBayesST", 
                    N=Ncar, burn.in=burn.in.car, thin=thin)
summary(M3st_ar)


M3st_ar2 <- Bcartime(formula=f3, data=engdeaths, scol=scol, tcol=tcol, 
                     W=Weng, model="ar", AR=2, family="gaussian", package="CARBayesST", 
                     N=Ncar, burn.in=burn.in.car, thin=thin)
summary(M3st_ar2)


model <- c("bym", "ar1")
M3inla.bym.ar1 <- Bcartime(data=engdeaths, formula=f3, W=Weng, scol =scol, tcol=tcol,  
                           model=model,  family="gaussian", package="inla", 
                            N=inlaN) 
summary(M3inla.bym.ar1)



a <- rbind(M3st_linear$mchoice,  M3st_anova$mchoice,
           M3st_ar$mchoice, M3st_ar2$mchoice)
a <- a[, 1:4]
a <- a[, c(2, 1, 4, 3)]
rownames(a) <- c("Linear",  "Anova",  "AR (1)", "AR (2)")
colnames(a) <- c("pDIC", "DIC", "pWAIC", "WAIC")
table11 <- a
table11
dput(table11, file=paste0(tablepath, "table11.txt"))

##
M3s <- Bcartime(formula=f3, data=engtotals, scol=scol, W=Weng, model="leroux", family="gaussian", N=Ncar, burn.in=burn.in.car, thin=thin)

summary(M3s)

end.time <- proc.time()[3]
comp.time<-end.time-start.time
# comp.time<-fancy.time(comp.time)
print(comp.time)

