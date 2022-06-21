# Print method for bmstdr objects. 
#' Provides basic information regarding the fitted model. 
#' @param x A bmstdr model fit object. 
#' @param digits How many significant digits after the decimal to print, defaults to 3.
#' @param ... Any other additional arguments 
## #' @seealso \link{\code{summary}} for summary of model fitting and   \link{\code{plot}} for plotting
## #' @seealso \link{\code{fitted}} for fitted values and \link{\code{resid}} for residuals. 
#' @return No return value, called for side effects. 
#' @method print bmstdr
#' @export
print.bmstdr <- function(x, digits=3, ...)
{
  if (x$package == "none") message("\n The ",  x$model, " model has been fitted using bmstdr code in R. \n")
  else message("\n The ",  x$model, " model has been fitted using the ", x$package, ".\n")
  
  message("Call:\n")
  print(x$call)
  
  message("\nModel formula\n")
  print(x$formula)
  message("\n")
  
  message("\nCoefficients:\n")
  print(round(x$params[, 1], digits = digits))
  message("\n")
  
  if (exists("computation.time", x))  { #  Validation has been performed  
    cat(x$computation.time, "\n")
  }
  
}
# Coefficient method for bmstdr objects. 
#' Prints and returns the estimates of the coefficients
#' @param object A bmstdr model fit object. 
#' @param digits How many significant digits after the decimal to print, defaults to 3.
#' @param ... Any other additional arguments 
#' @return Coefficients are returned as a data frame preserving the names of 
#' the covariates  
#' @method coef bmstdr
#' @export
coef.bmstdr <- function(object, digits=3, ...)
{
  # message("\nCoefficients:\n")
  # print(round(x$params[, 1], digits = digits))
  # message("\n")
  coefs <- object$params[, 1, drop=FALSE] 
  round(coefs, digits = digits)
}
# terms method for bmstdr objects. 
#' Prints the terms
#' @param x A bmstdr model fit object. 
#' @param ... Any other additional arguments 
#' @return Terms in the model formula
#' @method terms bmstdr
#' @export
terms.bmstdr <- function(x,...)
{
  terms(x$formula)
}
## Summary method for bmstdr objects. 
#' Provides basic summaries of model fitting. 
#' @param object A bmstdr model fit object. 
#' @param digits How many significant digits after the decimal to print, defaults to 3. 
#' @param ... Any other additional arguments 
## #' @seealso \link{\code{print}} for basic information regarding the fitted model. 
## #' @seealso \link{\code{plot}} for plotting and  \link{\code{fitted}} for 
## #' fitted values and \link{\code{resid}} for residuals.
#' @return No return value, called for side effects. 
#' @method summary bmstdr 
#' @rdname summary
#' @export
summary.bmstdr <- function(object, digits=3, ...)
{
 
  if (object$package == "none") message("\n The ",  object$model, " model has been fitted using bmstdr code in R. \n")
  else cat("\n The ",  object$model, " model has been fitted using the ", object$package, " package.\n")
  
  cat("Call:\n")
  print(object$call)

  
  if (exists("computation.time", object))  { #  Validation has been performed  
    cat(object$computation.time)
  }
  
  cat("\nModel formula\n")
  print(object$formula)
  cat("\n")
  
  cat("\nParameter Estimates:\n")
  print(round(object$params, digits = digits))
  
  if (exists("mchoice", object))  { # Model choice has been performed  
    cat("\nModel Choice Statistics:\n")
    print(round(object$mchoice, digits = digits))
  }
  
  if (exists("stats", object))  { #  Validation has been performed  
    cat("\nValidation Statistics:\n")
    # print(round(unlist(object$stats), digits = digits))
    print(object$stats)
  }
  
}
#' Extract fitted values from bmstdr objects. 
#' @param object A bmstdr model fit object. 
#' @param ... Any other additional arguments. 
#' @return Returns a vector of fitted values.
#' @method fitted bmstdr
#' @rdname fitted.bmstdr
#' @export
fitted.bmstdr <- function(object, ...)
{
if (object$scale.transform !="NONE") {  
  message("\n Note that the residuals are provided on the transformed scale. 
    Please see the scale.transform argument.\n")
} 
if (exists("stats", object)) { 
  message("Validation has been performed\n")
message("The fitted values include the validation observations as well. 
Expect the return value to be of the same length as the supplied data frame. \n")
}
  
object$fitteds
} 
#' Plot method for  bmstdr objects. 
#' @param x A bmstdr model fit object. 
#' @param segments TRUE or FALSE. It decides whether to draw the prediction intervals
#' as line segments. 
#' @param ... Any other additional arguments. 
#' @return  It plots the observed values on the original scale 
#' against the predictions and the 95\% prediction intervals if validation has been 
#' performed. It then plots the residuals against fitted values. It then applies 
#' plotting method to the model fitted object as returned by the chosen named package. 
#' There is no return value.   
#' @method plot bmstdr
#' @rdname plot.bmstdr 
#' @export
plot.bmstdr <- function(x, segments=TRUE, ...) { 
  old.par <- par(no.readonly = TRUE)
  on.exit(par(old.par))
  if (x$package == "spTimer")  {   
    plot(x$fit) 
    par(mfrow=c(1,1))
    par(ask=F)
  } else { 
    message("\n No other plots implemented for this model fitting method.\n")
  }  
  v <- residuals(x)
  u <- fitted(x)
  adf <- data.frame(residvals=v, fitvals=u)
  ndf <- na.omit(adf)
  p <- ggplot() + 
    geom_point(data=ndf, aes(x=fitvals, y=residvals), size=1) + 
    geom_abline(intercept=0, slope=0, col="blue") +
    labs(x="Fitted values", y="Residuals") 
  plot(p)              
  if (exists("stats", x))  { #  Validation has been performed  
    message("\nValidation Statistics:\n")
    print(round(unlist(x$stats), digits=3))
    # df <- x$yobs_preds
    # k <- ncol(df)
    # predsums <- df[, (k-4):k]
    # u <- all.vars(x$formula)
    # b <- as.character(u[1])
    # yobs <- df[, b]
    # allplots <- obs_v_pred_plot(yobs, predsums, segments = segments, summarystat = summarystat) 
    if (segments) plot(x$validationplots$pwithseg)
    else plot(x$validationplots$pwithoutseg)
    }   
}
#' Extract residuals from a bmstdr  fitted object. 
#' @param object A bmstdr model fit object. 
#' @param numbers a list with two components: sn=number of spatial locations 
#' tn=number of time points. Residuals will be assumed to follow the arrangement 
#' of the data frame - sorted by space and then time within space.
#' @return Returns a vector of residuals. If appropriate, it draws a 
#' time series plot of residuals. Otherwise, it draws a plot of residuals 
#' against observation numbers.     
#' @param ... Any other additional arguments. 
## #' @seealso \link{\code{print}} for basic information regarding the fitted model, 
## #' \link{\code{summary}} for summaries of model fitting,  
## #' \link{\code{fitted}} for extracting the fitted values, 
## #' \link{\code{plot}} for plotting. 
#' @method residuals bmstdr
#' @rdname residuals.bmstdr 
#' @export
residuals.bmstdr <- function(object, numbers=NULL, ...)
{
if (object$scale.transform !="NONE") {  
message("\n Note that the residuals are provided on the transformed scale. Please see the scale.transform argument.\n")
} 
  if (exists("stats", object)) { 
    message("Validation has been performed. The residuals include the validation observations as well. 
    Expect the return value to be of the same length as the supplied data frame. \n")
  }
  
 # object <- M5 
 a <- object$residuals
 message("\nSummary of the residuals\n")
 print(summary(a))
 sn <- object$sn 
 tn <- object$tn 
 if (tn==0) {
  # No spatio-temporal model has been fitted
  # read tn from the user as the supplied argument of this function
  sn <- numbers$sn
  tn <- numbers$tn
  if (length(tn)==0) { 
    message("tn has not been supplied in residuals and it is 
    not possible to figure this out. Hence a time series plot of the residuals
    has not been provided here.\n")
    tn <- 0
    }
 }
 if (tn>1) { 
 rdata <- data.frame(s.index=rep(1:sn, each=tn), Time=rep(1:tn, sn), residuals=a)
 rdata$s.index <- as.factor(rdata$s.index)
 rdata <- na.omit(rdata)
 rplot <- ggplot(data=rdata, aes(x=Time, y=residuals, group=s.index)) +
   geom_line() + 
   geom_abline(intercept = 0, slope = 0, col="red") + 
   labs(title= "Time series plot of residuals for each location", x="Time", 
        y = "Residuals", size=2.5) 
 plot(rplot) 
 } else if(tn==1) {
   rdata <- data.frame(s.index=1:sn,  residuals=a)
   # rdata$s.index <- as.factor(rdata$s.index)
   rdata <- na.omit(rdata)
   rplot <- ggplot(data=rdata, aes(x=s.index, y=residuals)) +
     geom_line() + 
     geom_point() + 
     geom_abline(intercept = 0, slope = 0, col="red") + 
     labs(title= "Plot of residuals against observation numbers", x="Observation number", 
          y = "Residuals", size=2.5) 
   plot(rplot)  
 }
 a 
} 
## Extract residuals from a bmstdr  fitted object. 
## @param x A bmstdr model fit object. 
## ## @seealso \link{\code{print}} for basic information regarding the fitted model, 
## ## \link{\code{summary}} for summaries of model fitting,  
## ## \link{\code{fitted}} for extracting the fitted values, 
## ## \link{\code{plot}} for plotting.  
# ## @method resid bmstdr
resid.bmstdr <- function(x)
{
 if (x$scale.transform !="NONE") {  
    message("\n Note that the residuals are provided on the transformed scale. 
    See the scale.transform argument.\n")
 }  
  if (exists("stats", x)) { 
   message("Validation has been performed.\n")
   message("The residuals include the validation observations as well. 
    Expect the return value to be of the same length as the supplied data frame. \n")
 }
  
  a <- x$residuals
  message("\nSummary of the residuals\n")
  print(summary(a))
  a 
} 
#' Is it a bmstdr model fitted object?
#' @param x Any R object. 
#' @return A TRUE/FALSE logical output.
#' @export
is.bmstdr <- function(x) inherits(x, "bmstdr")

