#' Temperature and salinity data from Argo floats in the North Atlantic Ocean at 
#' three layers of depth: surface (less than 50 meters),  mid-layer (between 475-525 meters) 
#' and deep (975 to 1025 meters) during 2003. 
#' @source \insertCite{SahuandChallenor2008;textual}{bmstdr}
#'  @format A data frame with 6978 rows and 11 columns:
#' \describe{
#'  \item{lon}{Longitude of the argo float}
#'  \item{lat}{Latitude of the argo float}
#'  \item{time}{Cumulative day of the year in 2003}
#'  \item{day}{Day within each month in 2003}
#'  \item{month}{Month in 2003}
#'  \item{temp}{Temperature recorded by the Argo float in degree Celsius}
#'   \item{sali}{Salinity in practical salinity units}
#'   \item{xlon}{Centered and scaled values of longitude at each depth}
#'   \item{xlat}{Centered and scaled values of latitude at each depth}
#'   \item{xinter}{Centered and scaled values of longitude times latitude at each depth}
#' }
# usethis::use_data(argo_floats_atlantic_2003, overwrite = TRUE)
#' @references
#' \insertAllCited{}
#' @examples
#'  head(argo_floats_atlantic_2003)
#'  # Data for the surface layer 
#'  surface <- argo_floats_atlantic_2003[argo_floats_atlantic_2003$depth==1, ] 
#'  # Data for the mid-layer 
#'  midlayer <- argo_floats_atlantic_2003[argo_floats_atlantic_2003$depth==2, ] 
#'  # Data at the deep ocean 
#'  deep <- argo_floats_atlantic_2003[argo_floats_atlantic_2003$depth==3, ]  
"argo_floats_atlantic_2003"
## argo_floats_atlantic_2003 <- argo_floats_atlatic_2003
