#' Average air pollution values from 28 sites in New York.
#'
#' @source See the NYdata set in spTimer package.
#'  Each data point is the mean of the available daily 8-hour maximum average ozone concentrations
#'  in parts per billion (ppb) at each of the
#'  28 sites. The daily values are for the month of July and August in 2006.
#' @format A vector with 28 real values.
## library(spTimer)
## attach(NYdata)
## ydata <- as.vector(tapply(o8hrmax, INDEX=s.index,  FUN=mean, na.rm=T))
ydata <- c(45.08065,  48.34758,  54.65677,  44.95774,  43.92903,  55.37726,  51.08589,
44.13661,  44.82355,  41.24903,  45.86645,  45.88557,  48.93226,  46.31767,  49.96500,
48.18387,  45.64532,  39.87500,  54.25150,  57.59550,  48.31339,  46.14145,  52.99596,
48.16083,  44.35355,  39.18065,  52.80557,  52.39814)
# usethis::use_data(ydata, overwrite = TRUE)
#' @examples
#' \dontrun{
#'  ydata
#' }
"ydata"
