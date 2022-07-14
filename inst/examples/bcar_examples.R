

# Set the validation row numbers 
vs <- sample(nrow(engtotals), 5)
# Total number of iterations 
N <- 60
# Number of burn in iterations 
burn.in <- 10
# Number of thinning iterations
thin <- 1

# Set the model formula for binomial distribution based modeling 
f1 <- noofhighweeks ~ jsa + log10(houseprice) + log(popdensity) + sqrt(no2)
## Independent error logistic regression
M1 <- Bcartime(formula = f1, data = engtotals, family = "binomial",
    trials = "nweek", N = N, burn.in = burn.in, thin = thin,
    verbose = TRUE)
summary(M1)
# Leroux model
M1.leroux <- Bcartime(formula = f1, data = engtotals, scol = "spaceid",
    model = "leroux", W = Weng, family = "binomial", trials = "nweek",
    N = N, burn.in = burn.in, thin = thin)
summary(M1.leroux)
# BYM model
M1.bym <- Bcartime(formula = f1, data = engtotals, scol = "spaceid",
    model = "bym", W = Weng, family = "binomial", trials = "nweek",
    N = N, burn.in = burn.in, thin = thin, verbose = FALSE)
summary(M1.bym)

# Validation for the Leroux model
M1.leroux.v <- Bcartime(formula = f1, data = engtotals, scol = "spaceid",
    model = "leroux", W = Weng, family = "binomial", trials = "nweek",
    validrows = vs, N = N, burn.in = burn.in, thin = thin, verbose = FALSE)
summary(M1.leroux.v)


## Poisson Distribution based models ####################################
# Model formula
f2 <- covid ~ offset(logEdeaths) + jsa + log10(houseprice) + log(popdensity) +
    sqrt(no2)

# Independent error Poisson regression
M2 <- Bcartime(formula = f2, data = engtotals, family = "poisson", N = N,
    burn.in = burn.in, thin = thin, verbose = FALSE)
summary(M2)
## Poisson regression with Leroux Model
M2.leroux <- Bcartime(formula = f2, data = engtotals, scol = "spaceid",
    model = "leroux", family = "poisson", W = Weng, N = N, burn.in = burn.in,
    thin = thin, verbose = FALSE)
summary(M2.leroux)
# Poisson regression with BYM Model
M2.bym <- Bcartime(formula = f2, data = engtotals, scol = "spaceid",
    model = "bym", family = "poisson", W = Weng, N = N, burn.in = burn.in,
    thin = thin)
summary(M2.bym)


## Gaussian distribution based models  ###############
f3 <- sqrt(no2) ~ jsa + log10(houseprice) + log(popdensity)

# Independent error model 
M3 <- Bcartime(formula = f3, data = engtotals, family = "gaussian", N = N,
    burn.in = burn.in, thin = thin, verbose = FALSE)
summary(M3)
# Leroux model 
M3.leroux <- Bcartime(formula = f3, data = engtotals, scol = "spaceid",
    model = "leroux", family = "gaussian", W = Weng, N = N, burn.in = burn.in,
    thin = thin, verbose = FALSE)
summary(M3.leroux)

## Validation
M3.leroux.v <- Bcartime(formula = f3, data = engtotals, scol = "spaceid",
    model = "leroux", family = "gaussian", W = Weng, N = N, burn.in = burn.in,
    thin = thin, validrows = vs, verbose = FALSE)
summary(M3.leroux.v)

\donttest{

## Spatio-temporal modeling ##################################################
head(engdeaths)
dim(engdeaths)
colnames(engdeaths)
vs <- sample(nrow(engdeaths), 5)


## Binomial distribution
engdeaths$nweek <- rep(1, nrow(engdeaths))
f1 <- highdeathsmr ~ jsa + log10(houseprice) + log(popdensity)

M1st_linear <- Bcartime(formula = f1, data = engdeaths, scol = "spaceid",
    tcol = "Weeknumber", trials = "nweek", W = Weng, model = "linear",
    family = "binomial", package = "CARBayesST", N = N, burn.in = burn.in,
    thin = thin, verbose = TRUE)
summary(M1st_linear)
M1st_sepspat <- Bcartime(formula = f1, data = engdeaths, scol = "spaceid",
    tcol = "Weeknumber", trials = "nweek", W = Weng, model = "sepspatial",
    family = "binomial", package = "CARBayesST", N = N, burn.in = burn.in,
    thin = thin, verbose = FALSE)
summary(M1st_sepspat)
M1st_ar <- Bcartime(formula = f1, data = engdeaths, scol = "spaceid",
    tcol = "Weeknumber", trials = "nweek", W = Weng, model = "ar", AR = 1,
    family = "binomial", package = "CARBayesST", N = N, burn.in = burn.in,
    thin = thin, verbose = FALSE)
summary(M1st_ar)
# Model validation
M1st_ar.v <- Bcartime(formula = f1, data = engdeaths, scol = "spaceid",
    tcol = "Weeknumber", trials = "nweek", W = Weng, model = "ar", AR = 1,
    family = "binomial", package = "CARBayesST", N = N, burn.in = burn.in,
    thin = thin, validrows = vs, verbose = FALSE)
summary(M1st_ar.v)


## Spatio temporal Poisson models###################################
colnames(engdeaths)
f2 <- covid ~ offset(logEdeaths) + jsa + log10(houseprice) + log(popdensity) +
    n0

M2st_linear <- Bcartime(formula = f2, data = engdeaths, scol = "spaceid",
    tcol = "Weeknumber", W = Weng, model = "linear", family = "poisson",
    package = "CARBayesST", N = N, burn.in = burn.in, thin = thin,
    verbose = FALSE)
summary(M2st_linear)
M2st_anova <- Bcartime(formula = f2, data = engdeaths, scol = "spaceid",
    tcol = "Weeknumber", W = Weng, model = "anova", family = "poisson",
    package = "CARBayesST", N = N, burn.in = burn.in, thin = thin,
    verbose = FALSE)
summary(M2st_anova)
M2st_anova_nointer <- Bcartime(formula = f2, data = engdeaths, scol = "spaceid",
    tcol = "Weeknumber", W = Weng, model = "anova", interaction = FALSE,
    family = "poisson", package = "CARBayesST", N = N, burn.in = burn.in,
    thin = thin, verbose = FALSE)
summary(M2st_anova_nointer)
M2st_sepspat <- Bcartime(formula = f2, data = engdeaths, scol = "spaceid",
    tcol = "Weeknumber", W = Weng, model = "sepspatial", family = "poisson",
    package = "CARBayesST", N = N, burn.in = burn.in, thin = thin,
    verbose = FALSE)
summary(M2st_sepspat)
M2st_ar <- Bcartime(formula = f2, data = engdeaths, scol = "spaceid",
    tcol = "Weeknumber", W = Weng, model = "ar", AR = 1, family = "poisson",
    package = "CARBayesST", N = N, burn.in = burn.in, thin = thin,
    verbose = FALSE)
summary(M2st_ar)
M2st_ar.v <- Bcartime(formula = f2, data = engdeaths, scol = "spaceid",
    tcol = "Weeknumber", W = Weng, model = "ar", family = "poisson",
    package = "CARBayesST", N = N, burn.in = burn.in, thin = thin,
    validrows = vs, verbose = FALSE)
M2st_anova.v <- Bcartime(formula = f2, data = engdeaths, scol = "spaceid",
    tcol = "Weeknumber", W = Weng, model = "anova", family = "poisson",
    package = "CARBayesST", N = N, burn.in = burn.in, thin = thin,
    validrows = vs, verbose = FALSE)
summary(M2st_ar.v)
summary(M2st_anova.v)


## Spatio-temporal Normal models ###############################
colnames(engdeaths)
f3 <- sqrt(no2) ~ jsa + log10(houseprice) + log(popdensity)

M3st_linear <- Bcartime(formula = f3, data = engdeaths, scol = "spaceid",
    tcol = "Weeknumber", W = Weng, model = "linear", family = "gaussian",
    package = "CARBayesST", N = N, burn.in = burn.in, thin = thin,
    verbose = FALSE)
summary(M3st_linear)
M3st_anova <- Bcartime(formula = f3, data = engdeaths, scol = "spaceid",
    tcol = "Weeknumber", W = Weng, model = "anova", family = "gaussian",
    package = "CARBayesST", N = N, burn.in = burn.in, thin = thin,
    verbose = FALSE)
summary(M3st_anova)
M3st_anova_nointer <- Bcartime(formula = f3, data = engdeaths, scol = "spaceid",
    tcol = "Weeknumber", W = Weng, model = "anova", interaction = FALSE,
    family = "gaussian", package = "CARBayesST", N = N, burn.in = burn.in,
    thin = thin, verbose = FALSE)
summary(M3st_anova_nointer)
M3st_ar <- Bcartime(formula = f3, data = engdeaths, scol = "spaceid",
    tcol = "Weeknumber", W = Weng, model = "ar", AR = 2, family = "gaussian",
    package = "CARBayesST", N = N, burn.in = burn.in, thin = thin,
    verbose = FALSE)
summary(M3st_ar)

# Execute the following examples if INLA is available
if (require(INLA)) {
    N <- 55
    burn.in <- 5
    thin <- 1
    vs <- sample(nrow(engtotals), 5)
    

    # Spatial Binomial GLM

    f1 <- noofhighweeks ~ jsa + log10(houseprice) + log(popdensity) +  sqrt(no2)

     M1.inla.bym <- Bcartime(data = engtotals, formula = f1,
        W = Weng, scol = "spaceid", model = c("bym"), package = "inla",
        family = "binomial", link="logit",  trials = "nweek", N = N, burn.in = burn.in,
        thin = thin)
     summary(M1.inla.bym)



    ## Spatial only Poisson

    f2inla <- covid ~ jsa + log10(houseprice) + log(popdensity) +
        sqrt(no2)

    M2.inla.bym <- Bcartime(data = engtotals, formula = f2inla,
        W = Weng, scol = "spaceid", offsetcol = "logEdeaths",
        model = c("bym"), package = "inla", link = "log",
        family = "poisson", N = N, burn.in = burn.in, thin = thin)
    summary(M2.inla.bym)



    ## Normal models

    f3 <- sqrt(no2) ~ jsa + log10(houseprice) + log(popdensity)

    ## Fit BYM with iid random effect
    M3.inla.bym <- Bcartime(data = engtotals, formula = f3,
        W = Weng, scol = "spaceid", model = c("bym"), package = "inla",
        family = "gaussian", N = N, burn.in = burn.in, thin = thin)
    summary(M3.inla.bym)

    # validation
    M3.inla.bym.v <- Bcartime(data = engtotals, formula = f3,
        W = Weng, scol = "spaceid", model = c("bym"), package = "inla",
        family = "gaussian",  validrows = vs, N = N, burn.in = burn.in,
        thin = thin)

    summary(M3.inla.bym.v)


    ###### Spatio-temporal INLA models

   
    f1 <- highdeathsmr ~ jsa + log10(houseprice) + log(popdensity)
    nweek <- rep(1, nrow(engdeaths))
    engdeaths$nweek <- rep(1, nrow(engdeaths))

    ## INLA Binomial
  
    model <- c("bym", "ar1")
    M1st_inla.bym <- Bcartime(data = engdeaths, formula = f1,
        W = Weng, scol = "spaceid", tcol = "Weeknumber",
        model = model, trials = "nweek", family = "binomial", link="logit", 
        package = "inla", N = N, burn.in = burn.in, thin = thin)
    summary(M1st_inla.bym)

    M1st_inla_v <- Bcartime(data = engdeaths, formula = f1,
        W = Weng, scol = "spaceid", tcol = "Weeknumber",
        offsetcol = "logEdeaths", model = model, trials = "nweek",
        family = "binomial",  link="logit", package = "inla", validrow = vs,
        N = N, burn.in = burn.in, thin = thin)
    summary(M1st_inla_v)


    model <- c("bym", "none")
    M1st_inla.bym.none <- Bcartime(data = engdeaths, formula = f1,
        W = Weng, scol = "spaceid", tcol = "Weeknumber",
        model = model, trials = "nweek", family = "binomial", link="logit", 
        package = "inla", N = N, burn.in = burn.in, thin = thin)
    summary(M1st_inla.bym.none)


    model <- c("bym")
    M1st_inla.bym.none <- Bcartime(data = engdeaths, formula = f1,
        W = Weng, scol = "spaceid", tcol = "Weeknumber",
        model = model, trials = "nweek", family = "binomial", link="logit", 
        package = "inla", N = N, burn.in = burn.in, thin = thin)
    summary(M1st_inla.bym.none)


    ## Poisson models
    f2inla <- covid ~ jsa + log10(houseprice) + log(popdensity) +
        n0 + n1 + n2 + n3

    model <- c("bym", "ar1")
    M2stinla <- Bcartime(data = engdeaths, formula = f2inla,
        W = Weng, scol = "spaceid", tcol = "Weeknumber",
        offsetcol = "logEdeaths", model = model, link = "log",
        family = "poisson", package = "inla", N = N, burn.in = burn.in,
        thin = thin)

    summary(M2stinla)

    M2stinla.v <- Bcartime(data = engdeaths, formula = f2inla,
        W = Weng, scol = "spaceid", tcol = "Weeknumber",
        offsetcol = "logEdeaths", model = model, link = "log",
        family = "poisson", package = "inla", validrows = vs,
        N = N, burn.in = burn.in, thin = thin)

    summary(M2stinla.v)

    ## Normal models

    f3 <- sqrt(no2) ~ jsa + log10(houseprice) + log(popdensity)

    model <- c("bym", "iid")
    M3inla.bym.iid <- Bcartime(data = engdeaths, formula = f3,
        W = Weng, scol = "spaceid", tcol = "Weeknumber",
        model = model, family = "gaussian", package = "inla",
        validrows = vs, N = N, burn.in = burn.in, thin = thin)
    summary(M3inla.bym.iid)

    model <- c("bym", "ar1")
    M3inla.bym.ar1 <- Bcartime(data = engdeaths, formula = f3,
        W = Weng, scol = "spaceid", tcol = "Weeknumber",
        model = model, family = "gaussian", package = "inla",
        validrows = vs, N = N, burn.in = burn.in, thin = thin)
    summary(M3inla.bym.ar1)

}
}

