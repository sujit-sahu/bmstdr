
N <- 10
burn.in <- 5
n.report <- 2
a <- Bspatial(formula = mpg ~ wt, data = mtcars, package = "none", model = "lm",
    N = N)
summary(a)
plot(a)
print(a)
b <- Bspatial(formula = mpg ~ disp + wt + qsec + drat, data = mtcars,
    validrows = c(8, 11, 12, 14, 18, 21, 24, 28), N = N)
#' print(b)
summary(b)
## Illustration with the nyspatial data set
head(nyspatial)
## Linear regression model fitting
M1 <- Bspatial(formula = yo3 ~ xmaxtemp + xwdsp + xrh, data = nyspatial,
    mchoice = TRUE, N = N)
print(M1)
plot(M1)
a <- residuals(M1)
summary(M1)
## Spatial model fitting
M2 <- Bspatial(model = "spat", formula = yo3 ~ xmaxtemp + xwdsp +
    xrh, data = nyspatial, coordtype = "utm", coords = 4:5, phi = 0.4,
    mchoice = TRUE, N = N)
names(M2)
print(M2)
plot(M2)
a <- residuals(M2)
summary(M2)
## Fit model 2 on the square root scale
M2root <- Bspatial(model = "spat", formula = yo3 ~ xmaxtemp + xwdsp + xrh,
    data = nyspatial, coordtype = "utm", coords = 4:5, scale.transform = "SQRT")
summary(M2root)
## Spatial model fitting using spBayes
M3 <- Bspatial(package = "spBayes", formula = yo3 ~ xmaxtemp + xwdsp + xrh,
    data = nyspatial, coordtype = "utm", coords = 4:5, prior.phi = c(0.005,
        2), mchoice = TRUE, N = N, burn.in = burn.in, n.report = n.report)
summary(M3)
\donttest{
# Spatial model fitting using stan (with a small number of iterations)
M4 <- Bspatial(package = "stan", formula = yo3 ~ xmaxtemp + xwdsp + xrh,
    data = nyspatial, coordtype = "utm", coords = 4:5, phi = 0.4, N = N,
    burn.in = burn.in, mchoice = TRUE)
summary(M4)


## K fold cross-validation for M2 only
set.seed(44)
x <- runif(n = 28)
u <- order(x)
# Here are the four folds
s1 <- u[1:7]
s2 <- u[8:14]
s3 <- u[15:21]
s4 <- u[22:28]
summary((1:28) - sort(c(s1, s2, s3, s4)))  ## check
v1 <- Bspatial(model = "spat", formula = yo3 ~ xmaxtemp + xwdsp + xrh,
    data = nyspatial, coordtype = "utm", coords = 4:5, validrows = s1,
    phi = 0.4, N = N)
v2 <- Bspatial(model = "spat", formula = yo3 ~ xmaxtemp + xwdsp + xrh,
    data = nyspatial, coordtype = "utm", coords = 4:5, validrows = s2,
    phi = 0.4, N = N)
v3 <- Bspatial(model = "spat", formula = yo3 ~ xmaxtemp + xwdsp + xrh,
    data = nyspatial, coordtype = "utm", coords = 4:5, validrows = s3,
    phi = 0.4, N = N)
v4 <- Bspatial(model = "spat", formula = yo3 ~ xmaxtemp + xwdsp + xrh,
    data = nyspatial, coordtype = "utm", coords = 4:5, validrows = s4,
    phi = 0.4, N = N)
M2.val.table <- cbind(unlist(v1$stats), unlist(v2$stats), unlist(v3$stats),
    unlist(v4$stats))
dimnames(M2.val.table)[[2]] <- paste("Fold", 1:4, sep = "")
round(M2.val.table, 3)

## Model validation
s <- c(1, 5, 10)
M1.v <- Bspatial(model = "lm", formula = yo3 ~ xmaxtemp + xwdsp + xrh,
    data = nyspatial, coordtype = "utm", coords = 4:5, validrows = s, N = N,
    burn.in = burn.in)
M2.v <- Bspatial(model = "spat", formula = yo3 ~ xmaxtemp + xwdsp + xrh,
    data = nyspatial, coordtype = "utm", coords = 4:5, validrows = s, phi = 0.4,
    N = N, burn.in = burn.in)
M3.v <- Bspatial(package = "spBayes", formula = yo3 ~ xmaxtemp + xwdsp +
    xrh, data = nyspatial, coordtype = "utm", coords = 4:5, validrows = s,
    prior.phi = c(0.005, 2), n.report = 2, N = N, burn.in = burn.in)
# Collect all the results
Mall.table <- cbind(unlist(M1.v$stats), unlist(M2.v$stats), unlist(M3.v$stats))
colnames(Mall.table) <- paste("M", c(1:3), sep = "")
round(Mall.table, 3)


if (require(INLA)) {
    N <- 10
    burn.in <- 5
    # Spatial model fitting using INLA
    M5 <- Bspatial(package = "inla", formula = yo3 ~ xmaxtemp + xwdsp + xrh,
        data = nyspatial, coordtype = "utm", coords = 4:5, mchoice = TRUE,
        N = N, burn.in = burn.in)
    summary(M5)
}
}
