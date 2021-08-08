## ----style, echo = FALSE, results = 'asis'------------------------------------
  BiocStyle::markdown()

## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>")

## ----setup,  eval=TRUE, echo=FALSE, include=FALSE-----------------------------
library(bmstdr)
library(ggplot2)
require(ggsn)
library(tidyr)
library(huxtable)
library(RColorBrewer)
library(akima)
knitr::opts_chunk$set(eval = F)
colpalette <- c("dodgerblue4",  "dodgerblue2",  "firebrick2",   "firebrick4",   "purple")     
figpath <- system.file("figs", package = "bmstdr") 
tablepath <- system.file('last3tables', package = 'bmstdr')

tabs <- lapply(list.files(system.file('txttables', package = 'bmstdr'), full.names = TRUE), dget)
# print(tabs)
# fns <- list.files(system.file('txttables', package = 'bmstdr')) 
# how the table data and file names match 
table1 <-  tabs[[1]]
table10 <- tabs[[2]]
table11 <- tabs[[3]]
table2 <-  tabs[[4]]
table3 <-  tabs[[5]]
table4 <- tabs[[6]]
table5 <- tabs[[7]]
table6 <- tabs[[8]]
table7 <- tabs[[9]]
table8 <- tabs[[10]]
table9 <- tabs[[11]]

## ----vsites3, echo=FALSE, eval=TRUE, fig.cap="25 fitted sites and 3 validation sites (numbered) in  New York", fig.width=6, fig.height=4----
nymap <- map_data(database="state",regions="new york")
s <- c(1, 5, 10)
fcoords <- nyspatial[-s, c("Longitude", "Latitude")]
vcoords <- nyspatial[s,  c("Longitude", "Latitude")]
library(tidyr)
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
   labs( title= "28 air pollution monitoring sites in New York", x="Longitude", y = "Latitude") +
   geom_text(aes(label=label, x=long, y=lat), data = label, vjust = "top", hjust = "right")  +
   # geom_rect(mapping=aes(xmin=-80.2, xmax=-77.3, ymin=41, ymax=41.6), color="black", fill=NA) + 
   geom_rect(mapping=aes(xmin=-78.7, xmax=-75.8, ymin=41, ymax=41.6), color="black", fill=NA) + 
   ggsn::scalebar(data =nymap, dist = 100, location = "bottomleft", transform=T, dist_unit = "km",
                  st.dist = .05, st.size = 5, height = .06, st.bottom=T, model="WGS84") +
   ggsn::north(data=nymap, location="topleft", symbol=12) 
vsites3

## ----chunk1, echo=TRUE, eval=FALSE--------------------------------------------
#  M1 <- Bspatial(formula=yo3~xmaxtemp+xwdsp+xrh, data=nyspatial, mchoice=T)

## ----m2spatial, echo=TRUE, eval=FALSE-----------------------------------------
#  M2 <- Bspatial(model="spat", formula=yo3~xmaxtemp+xwdsp+xrh, data=nyspatial,
#        coordtype="utm", coords=4:5, phi=0.4, mchoice=T)

## ---- echo=TRUE, eval=FALSE---------------------------------------------------
#  asave <- phichoice_sp()
#  asave

## ---- echo=TRUE, eval=FALSE---------------------------------------------------
#  M3 <- Bspatial(package="spBayes", formula=yo3~xmaxtemp+xwdsp+xrh, data=nyspatial,
#                 coordtype="utm", coords=4:5, prior.phi=c(0.005, 2), mchoice=T)
#  M4 <- Bspatial(package="stan", formula=yo3~xmaxtemp+xwdsp+xrh, data=nyspatial,
#                 coordtype="utm", coords=4:5,phi=0.4, mchoice=T)
#  M5  <- Bspatial(package="inla",formula=yo3~xmaxtemp+xwdsp+xrh, data=nyspatial,
#                 coordtype="utm", coords=4:5, mchoice=T)

## ---- echo=FALSE, eval=FALSE--------------------------------------------------
#  a3 <- Bmchoice(case="MC.sigma2.unknown", y=ydata)
#  # Now organize the all the results for forming Table 1.
#  a5 <- rep(NA, 11)
#  a5[c(1, 3, 5, 7, 9:11)] <- unlist(M5$mchoice)
#  table1 <- cbind.data.frame(unlist(a3), M1$mchoice, M2$mchoice, M3$mchoice, M4$mchoice,  a5)
#  colnames(table1) <- paste("M", 0:5, sep="")
#  round(table1,  2)
#  dput(table1, file=paste0(path, "table1.txt"))

## ----tab, echo=FALSE, eval=TRUE-----------------------------------------------
   table1 %>%
   as_hux(add_colnames = FALSE) %>%
   set_number_format(2)     %>%
   map_text_color(by_cols("darkred", "blue", "darkgreen")) %>%
   add_colnames("Criteria") %>%
   set_header_rows(1, TRUE) %>% 
   add_rownames() %>%
   set_all_borders(1)  %>%
   set_bold(1, everywhere) %>% 
   style_headers(bold = TRUE, text_color = "red") %>% 
   set_caption('(#tab:mchoicenyspatial) Model choice criteria for various models fitted to the nyspatial data set.')


## ---- echo=TRUE, eval=FALSE---------------------------------------------------
#  Bmchoice(case="MC.sigma2.unknown", y=ydata).

## ---- echo=TRUE, eval=FALSE---------------------------------------------------
#  s <- c(8,11,12,14,18,21,24,28)
#  f1 <- yo3~xmaxtemp+xwdsp+xrh
#  M1.v <-  Bspatial(package="none", model="lm", formula=f1,
#                    data=nyspatial, validrows=s)
#  M2.v <- Bspatial(package="none", model="spat", formula=f1,
#          data=nyspatial,   coordtype="utm", coords=4:5,phi=0.4,  validrows=s)
#  M3.v <-  Bspatial(package="spBayes", prior.phi=c(0.005, 2), formula=f1,
#          data=nyspatial,   coordtype="utm", coords=4:5, validrows=s)
#  M4.v  <- Bspatial(package="stan",formula=f1,
#      data=nyspatial,   coordtype="utm", coords=4:5,phi=0.4 , validrows=s)
#  M5.v  <- Bspatial(package="inla", formula=f1, data=nyspatial,
#          coordtype="utm", coords=4:5, validrows=s)

## ---- echo=FALSE, eval=TRUE---------------------------------------------------
   table2[-4, ] %>%
   as_hux(add_colnames = FALSE) %>%
   set_number_format(3)     %>%
   map_text_color(by_cols("darkred", "blue", "darkgreen", "purple")) %>%
   add_colnames("Criteria") %>%
   set_header_rows(1, TRUE) %>% 
   add_rownames() %>%
   set_all_borders(1)  %>%
   set_bold(1, everywhere) %>% 
   style_headers(bold = TRUE, text_color = "red")%>% 
   set_caption('(#tab:mvalidnyspatial) Model validation statistics for the five models fitted to the nyspatial data set.')

## ---- echo=TRUE, eval=TRUE----------------------------------------------------
set.seed(44)
x <- runif(n=28)
u <- order(x)
s1 <- u[1:7]
s2 <- u[8:14]
s3 <- u[15:21]
s4 <- u[22:28]

## ---- echo=FALSE, eval=TRUE---------------------------------------------------
   table3 %>%
   as_hux(add_colnames = FALSE) %>%
   set_number_format(3)     %>%
   map_text_color(by_cols("darkred", "blue", "darkgreen", "purple")) %>%
   add_colnames("Criteria") %>%
   set_header_rows(1, TRUE) %>% 
   add_rownames() %>%
   set_all_borders(1)  %>%
   set_bold(1, everywhere) %>% 
   style_headers(bold = TRUE, text_color = "red") %>%
   set_caption('(#tab:m2-4-fold-validnyspatial) 4-fold cross-validation statistics for model M2 fitted to the nyspatial data set.')


## ----valplot1, echo=T, eval=T, message=FALSE, results='hide', fig.cap="Prediction against observation plot with the prediction intervals included. The `in/out' symbol in the plot indicates whether or not a prediction interval incudes the 45 degree line."----
M2.v3 <- Bspatial(model="spat", formula=yo3~xmaxtemp+xwdsp+xrh, data=nyspatial, 
               coordtype="utm", coords=4:5, validrows= s3, phi=0.4, verbose = FALSE)

## ---- echo=TRUE, eval=FALSE, message=FALSE, results='hide'--------------------
#  names(M2.v3)
#  psums <- get_validation_summaries(M2.v3$valpreds)
#  names(psums)
#  a <- obs_v_pred_plot(yobs=M2.v3$yobs_preds$yo3, predsums=psums, segments=FALSE, summarystat = "mean" )

## ---- echo=TRUE, eval=TRUE, message=FALSE, results='hide'---------------------
f2 <- y8hrmax~xmaxtemp+xwdsp+xrh
M1 <- Bsptime(model="lm", formula=f2, data=nysptime, scale.transform = "SQRT")
M2 <- Bsptime(model="separable", formula=f2, data=nysptime, scale.transform = "SQRT",
              coordtype="utm", coords=4:5)

## ----residM2, echo=TRUE, eval=TRUE, fig.cap="A multiple time series plot of residuals"----
a <- residuals(M2)

## ---- echo=TRUE, eval=FALSE, message=FALSE, results='hide'--------------------
#  a <- residuals(M1, numbers=list(sn=28, tn=62))

## ---- echo=TRUE, eval=FALSE---------------------------------------------------
#  M2$phi.s; M2$phi.t

## ----valplot21, echo=TRUE, eval=TRUE, message=FALSE, results='hide', fig.show='hide'----
valids <-  c(1, 5, 10)
vrows <-  which(nysptime$s.index%in% valids)
M2.1 <- Bsptime(model="separable",  formula=f2, data=nysptime, 
             validrows=vrows, coordtype="utm", coords=4:5,
             phi.s=0.005, phi.t=0.05, scale.transform = "SQRT")

## ---- echo=TRUE, eval=TRUE, message=FALSE-------------------------------------
summary(M2.1)

## ---- echo=TRUE, eval=TRUE,message=FALSE, results='hide'----------------------
M3 <- Bsptime(package="spTimer", formula=f2, data=nysptime, n.report=5, 
              coordtype="utm", coords=4:5, scale.transform = "SQRT")

## ---- echo=T, eval=T----------------------------------------------------------
set.seed(44)
tn <- 62
sn <- 28
valids <- c(1, 5, 10)
validt <- sort(sample(1:tn, size=31))
vrows <- getvalidrows(sn=sn, tn=tn, valids=valids, validt=validt)

## ---- echo=TRUE, eval=TRUE, message=FALSE, results='hide', fig.show='hide'----
M31 <- Bsptime(package="spTimer",formula=f2, data=nysptime, 
               coordtype="utm", coords=4:5,
               validrows=vrows, model="GP", report=5, 
               mchoice=F, scale.transform = "NONE")

## ---- echo=T, eval=T----------------------------------------------------------
modfit <- M31$fit
fitall <- data.frame(modfit$fitted)
fitall$s.index <- rep(1:sn, each=tn)
library(spTimer)
vdat <- spT.subset(data=nysptime, var.name=c("s.index"), s=valids)
fitvalid <- spT.subset(data=fitall, var.name=c("s.index"), s=valids)
fitvalid$low <- fitvalid$Mean - 1.96 * fitvalid$SD
fitvalid$up <- fitvalid$Mean + 1.96 * fitvalid$SD
fitvalid$yobs <- vdat$y8hrmax
yobs <- matrix(fitvalid$yobs, byrow=T, ncol=tn)
y.valids.low <- matrix(fitvalid$low, byrow=T, ncol=tn)
y.valids.med <- matrix(fitvalid$Mean, byrow=T, ncol=tn)
y.valids.up <- matrix(fitvalid$up, byrow=T, ncol=tn)

## ----valid3sites, echo=T, eval=T, fig.height=8, fig.cap="Time series of observed and predicted values at three sites."----
p1 <- fig11.13.plot(yobs[1, ], y.valids.low[1, ], y.valids.med[1, ], 
                    y.valids.up[1, ], misst=validt)
p1 <- p1 + ggtitle("Validation for Site 1")
p2 <- fig11.13.plot(yobs[2, ], y.valids.low[2, ], y.valids.med[2, ], 
                    y.valids.up[2, ], misst=validt)
p2 <- p2 + ggtitle("Validation for Site 5")
p3 <- fig11.13.plot(yobs[3, ], y.valids.low[3, ], y.valids.med[3, ], 
                    y.valids.up[3, ], misst=validt)
p3 <- p3 + ggtitle("Validation for Site 10")
library(ggpubr)
ggarrange(p1, p2, p3, common.legend = TRUE, legend = "top", nrow = 3, ncol = 1)

## ---- echo=T, eval=T----------------------------------------------------------
sitemeans <- function(a, sn, tn=62) { 
   u <- matrix(a, nrow=sn, ncol=tn, byrow=T)
   b <- apply(u, 1, mean)
   as.vector(b)
}

## ---- echo=T, eval=T, message=FALSE, results='hide'---------------------------
post <- M3$fit
gpred <- predict(post, newdata=gridnysptime, newcoords=~Longitude+Latitude)
u <- gpred$pred.samples
v <- apply(u, 2, sitemeans, sn=100)
a <- get_parameter_estimates(t(v)) 
b <- data.frame(gridnyspatial[, 1:5], a) 

## ---- echo=TRUE, eval=TRUE, message=FALSE, results='hide'---------------------
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

## ---- echo=TRUE, eval=TRUE, message=FALSE, results='hide'---------------------
b <- rbind(b, fits)

## ---- echo=TRUE, eval=TRUE, message=FALSE, results='hide'---------------------
coord <- nyspatial[, c("Longitude","Latitude")]
library(akima)
xo <- seq(from=min(coord$Longitude)-0.5, to = max(coord$Longitude)+0.8, length=200)
yo <- seq(from=min(coord$Latitude)-0.25, to = max(coord$Latitude)+0.8, length=200)
surf <- interp(b$Longitude, b$Latitude, b$mean,  xo=xo, yo=yo)
v <- fnc.delete.map.XYZ(xyz=surf)

interp1 <- data.frame(long = v$x, v$z )
names(interp1)[1:length(v$y)+1] <- v$y
library(tidyr)
interp1 <- gather(interp1,key = lat,value =Predicted,-long,convert = TRUE)
library(ggplot2)
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

## ---- echo=F, eval=T, fig.height=6, fig.width=8, fig.cap="Predicted map of average ozone air pollution in New York"----
P

## ---- echo=TRUE, eval=FALSE---------------------------------------------------
#  M1.c <- Bsptime(model="lm", formula=f2, data=nysptime,
#            scale.transform = "SQRT", mchoice=T)
#  M2.c <- Bsptime(model="separable",  formula=f2, data=nysptime,
#            coordtype="utm", coords=4:5, phi.s=0.005, phi.t=0.05,
#            scale.transform = "SQRT", mchoice=T)
#  M3.c <- Bsptime(package="spTimer", model="GP",
#          formula=f2, data=nysptime, coordtype="utm",
#          coords=4:5, scale.transform = "SQRT",
#          mchoice=T, N=5000)
#  M4.c <- Bsptime(package="stan",formula=f2, data=nysptime,
#          coordtype="utm", coords=4:5, scale.transform = "SQRT",
#          N=1500, burn.in=500, mchoice=T, verbose = F)

## ---- echo=FALSE, eval=TRUE---------------------------------------------------
   colnames(table4) <- c("M1", "M2", "M3", "M4")
   table4 %>%
   as_hux(add_colnames = FALSE) %>%
   set_number_format(2)     %>%
   map_text_color(by_cols("darkred", "blue", "darkgreen", "purple")) %>%
   add_colnames(c("Criteria")) %>%
   set_header_rows(1, TRUE) %>% 
   add_rownames() %>%
   set_all_borders(1)  %>%
   set_bold(1, everywhere) %>% 
   style_headers(bold = TRUE, text_color = "red") %>% 
   set_caption('(#tab:mchoice-m1-m4) Model choice  criteria for the four spatio-temporal models M1 to M4.')

## ---- echo=TRUE, eval=FALSE, message=FALSE, results='hide', fig.show='hide'----
#  M5 <- Bsptime(package="spTimer", model="AR", formula=f2, data=nysptime,
#                  coordtype="utm", coords=4:5, scale.transform = "SQRT",
#                  mchoice=T,  validrows = vrows)
#  

## ---- echo=TRUE, eval=FALSE, message=FALSE, results='hide'--------------------
#  M6 <- Bsptime(package="inla", model="AR", formula=f2, data=nysptime,
#          coordtype="utm", coords=4:5, scale.transform = "SQRT",
#          mchoice=T, validrows=vrows)

## ---- echo=FALSE, eval=TRUE---------------------------------------------------
   table5 %>%
   as_hux(add_colnames = FALSE) %>%
   set_number_format(2)     %>%
   map_text_color(by_cols("darkred", "blue", "darkgreen")) %>%
   add_colnames("Criteria") %>%
   set_header_rows(1, TRUE) %>% 
   add_rownames() %>%
   set_all_borders(1)  %>%
   set_bold(1, everywhere) %>% 
   style_headers(bold = TRUE, text_color = "red") %>%
   set_caption('(#tab:mcvm5m6) Model choice and validation statistics for the two AR models M5 and M6.')

## ---- echo=FALSE, eval=TRUE---------------------------------------------------
   table6 %>%
   as_hux(add_colnames = FALSE) %>%
   set_number_format(3)     %>%
   map_text_color(by_cols("black", "darkred", "blue", "darkgreen", "purple", "darkred", "blue", "darkgreen", "purple")) %>%
   add_colnames("Criteria") %>%
   set_header_rows(1, TRUE) %>% 
   add_rownames() %>%
   set_all_borders(1)  %>%
   set_bold(1, everywhere) %>% 
   style_headers(bold = TRUE, text_color = "red") %>% 
   set_caption('(#tab:mparm-m5-m6) Parameter estimates for the two AR models M5 (first four columns) and M6 (last four columns).')

## ---- echo=TRUE, eval=TRUE, message=FALSE, results='hide'---------------------
library(spTDyn)
f3 <- y8hrmax~ xmaxtemp + sp(xmaxtemp)+ tp(xwdsp) + xrh
M7 <- Bsptime(package="sptDyn", model="GP", formula=f3, data=nysptime, 
      coordtype="utm", coords=4:5, scale.transform = "SQRT", n.report=2)

## ----speffects,  echo=TRUE, eval=TRUE, message=FALSE,  fig.cap="Spatial effects of maximum temperature from model M7"----
out <- M7$fit
dim(out$betasp)
a <- out$betasp
u <- c(t(out$betasp))
sn <- nrow(a)
itmax <- ncol(a)
v <- rep(1:sn, each=itmax)
d <- data.frame(site=as.factor(v), sp = u)
p <- ggplot(data=d, aes(x=site, y=sp)) + 
   geom_boxplot(outlier.colour="black", outlier.shape=1,
                outlier.size=0.5) +
   geom_abline(intercept=0, slope=0, color="blue") + 
   labs(title= "Spatial effects of maximum temperature", x="Site", y = "Effects", size=2.5) 
p 

## ----temporaleffects, echo=TRUE, eval=TRUE, message=FALSE, fig.cap="Temporal effects  of wind spped from model M7"----
b <- out$betatp
tn <- nrow(b)
itmax <- ncol(b)
tids <- 1:tn 
stat <- apply(b[tids,], 1, quantile, prob=c(0.025,0.5,0.975))
tstat <- data.frame(tids, t(stat))
dimnames(tstat)[[2]] <- c("Days", "low", "median", "up")
# head(tstat)
yr <- c(min(c(stat)),max(c(stat)))
p <- ggplot(data=tstat, aes(x=Days, y=median)) + 
   geom_point(size=3) + 
   ylim(yr) + 
   geom_segment(data=tstat, aes(x=Days, y=median, xend=Days, yend=low), linetype=1) +
   geom_segment(data=tstat, aes(x=Days, y=median, xend=Days, yend=up), linetype=1) +
   geom_abline(intercept=0, slope=0, col="blue") +
   labs(title="Temporal effects of wind speed", x="Days", y="Temporal effects") 
p 

