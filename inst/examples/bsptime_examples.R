
# Set the total number of iterations 
N <- 45
# Set the total number of burn-in iterations 
burn.in <- 5
# How many times to report progress 
n.report <- 2
# Model formula used in most model fitting  
f2 <- y8hrmax ~ xmaxtemp + xwdsp + xrh
# Check out the data set 
head(nysptime)
## Fit linear regression model 
M1 <- Bsptime(model = "lm", data = nysptime, formula = f2,
    scale.transform = "SQRT", N = N, burn.in = burn.in, mchoice = TRUE)
names(M1)
plot(M1)
print(M1)
summary(M1)
a <- residuals(M1, numbers = list(sn = 28, tn = 62))
M2 <- Bsptime(model = "separable", data = nysptime, formula = f2,
    coordtype = "utm", coords = 4:5, mchoice = TRUE, scale.transform = "SQRT",
    N = N, burn.in = burn.in)
names(M2)
plot(M2)
print(M2)
summary(M2)
b <- residuals(M2)
# Spatio-temporal model fitting and validation
valids <- c(8, 11)
vrows <- which(nysptime$s.index %in% valids)
## Fit separable spatio-temporal model 
M2.1 <- Bsptime(model = "separable", formula = f2, data = nysptime,
    validrows = vrows, coordtype = "utm", coords = 4:5, phi.s = 0.005,
    phi.t = 0.05, scale.transform = "SQRT", N = N)
summary(M2.1)
plot(M2.1)
# Use spTimer to fit independent GP model 
M3 <- Bsptime(package = "spTimer", formula = f2, data = nysptime,
    coordtype = "utm", coords = 4:5, scale.transform = "SQRT", mchoice = TRUE,
    N = N, burn.in = burn.in, n.report = 2)
summary(M3)
\donttest{
valids <- c(1, 5, 10)
validt <- sort(sample(1:62, size = 31))
vrows <- getvalidrows(sn = 28, tn = 62, valids = valids, validt = validt)
ymat <- matrix(nysptime$y8hrmax, byrow = TRUE, ncol = 62)
yholdout <- ymat[valids, validt]
# Perform validation 
M31 <- Bsptime(package = "spTimer", formula = f2, data = nysptime,
    coordtype = "utm", coords = 4:5, validrows = vrows, model = "GP",
    scale.transform = "NONE", N = N, burn.in = burn.in, n.report = 2)
summary(M31)
modfit <- M31$fit
## Extract the fits for the validation sites
fitall <- data.frame(modfit$fitted)
head(fitall)
tn <- 62
fitall$s.index <- rep(1:28, each = tn)
library(spTimer)
vdat <- spT.subset(data = nysptime, var.name = c("s.index"), s = valids)
fitvalid <- spT.subset(data = fitall, var.name = c("s.index"), s = valids)
head(fitvalid)
fitvalid$low <- fitvalid$Mean - 1.96 * fitvalid$SD
fitvalid$up <- fitvalid$Mean + 1.96 * fitvalid$SD
fitvalid$yobs <- sqrt(vdat$y8hrmax)
fitvalid$yobs <- vdat$y8hrmax
yobs <- matrix(fitvalid$yobs, byrow = TRUE, ncol = tn)
y.valids.low <- matrix(fitvalid$low, byrow = TRUE, ncol = tn)
y.valids.med <- matrix(fitvalid$Mean, byrow = TRUE, ncol = tn)
y.valids.up <- matrix(fitvalid$up, byrow = TRUE, ncol = tn)
library(ggplot2)
p1 <- fig11.13.plot(yobs[1, ], y.valids.low[1, ], y.valids.med[1, ],
    y.valids.up[1, ], misst = validt)
p1 <- p1 + ggtitle("Validation for Site 1")
p1
p2 <- fig11.13.plot(yobs[2, ], y.valids.low[2, ], y.valids.med[2, ],
    y.valids.up[2, ], misst = validt)
p2 <- p2 + ggtitle("Validation for Site 5")
p2
p3 <- fig11.13.plot(yobs[3, ], y.valids.low[3, ], y.valids.med[3, ],
    y.valids.up[3, ], misst = validt)
p3 <- p3 + ggtitle("Validation for Site 10")
p3

## Independent marginal GP model fitting using rstan

M4 <- Bsptime(package = "stan", formula = f2, data = nysptime,
    coordtype = "utm", coords = 4:5, N = N, burn.in = burn.in,
    verbose = FALSE)
summary(M4)

# Spatio-temporal hierarchical auto-regressive modeling useing spTimer 
M5 <- Bsptime(package = "spTimer", model = "AR", formula = f2, data = nysptime,
    coordtype = "utm", coords = 4:5, scale.transform = "SQRT", mchoice = TRUE,
    n.report = n.report, N = N, burn.in = burn.in)
summary(M5)
a <- residuals(M5)

## Spatio-temporal dynamic model fitting using spTDyn
library(spTDyn)

f3 <- y8hrmax ~ xmaxtemp + sp(xmaxtemp) + tp(xwdsp) + xrh
M7 <- Bsptime(package = "sptDyn", model = "GP", formula = f3, data = nysptime,
    coordtype = "utm", coords = 4:5, scale.transform = "SQRT", mchoice = TRUE,
    N = N, burn.in = burn.in, n.report = n.report)
summary(M7)

# Dynamic Model fitting using spBayes 
M8 <- Bsptime(package = "spBayes", formula = f2, data = nysptime,
    prior.sigma2 = c(2, 25), prior.tau2 = c(2, 25), prior.sigma.eta = c(2,
        0.001), coordtype = "utm", coords = 4:5, scale.transform = "SQRT",
    N = N, burn.in = burn.in, n.report = n.report)
summary(M8)

## Gussian Predictive Process based model fitting using spTimer 
M9 <- Bsptime(package = "spTimer", model = "GPP", g_size = 5, formula = f2,
    data = nysptime, coordtype = "utm", coords = 4:5, scale.transform = "SQRT",
    N = N, burn.in = burn.in, n.report = n.report)
summary(M9)

# This INLA run may take a long time
if (require(INLA) & require(inlabru)) {
    f2 <- y8hrmax ~ xmaxtemp + xwdsp + xrh
    M6 <- Bsptime(package = "inla", model = "AR", formula = f2, data = nysptime,
        coordtype = "utm", coords = 4:5, scale.transform = "SQRT", 
        offset = c(100, 200), max.edge = c(500, 10000),
        mchoice = TRUE, plotit=TRUE)
    # Takes a minute
    summary(M6)
}
}

