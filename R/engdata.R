#' Total number of weekly Covid-19 deaths and cases in the 313 local 
#' Local Authority Districts, Counties  and Unitary Authorities (LADCUA) in 
#' England during the first peak from March 13 to July 31, 2020. 
#' @source \insertCite{SahuBohning2021;textual}{bmstdr}. 
#'  @format A data frame with 313 rows and 19 columns:
#' \describe{
#'   \item{Areacode}{Areacode identifier of the 313 
#'   Local Authority Districts, Counties  and Unitary Authorities (LADCUA)}
#'   \item{mapid}{A numeric column identifying the map area needed for plotting}
#'   \item{spaceid}{A numeric variable taking value 1 to 313 identifying the LADCUA's}
#'   \item{Region}{Identifies one of the 9 English regions}
#'   \item{popn}{Population number in mid-2019}
#'   \item{jsa}{Percentage of the working age  population receiving job-seekers allowance 
#'    during January 2020}
#'   \item{houseprice}{Median house price  in March 2020}
#'   \item{popdensity}{Population density in mid-2019}    
#'   \item{startdate}{Start date of the week}
#'   \item{Weeknumber}{Week numbers 11 to 30}
#'   \item{no2}{Estimated average value of NO2 at the centroid of the LADCUA in that week}
#'   \item{covid}{Number of Covid-19 deaths within 28 days of a positive test}
#'   \item{allcause}{Total number deaths}
#'   \item{noofcases}{Total number of cases}
#'   \item{Edeaths}{Expected number of Covid-19 deaths. See Sahu and 
#'   Bohning (2021) for methodology.} 
#'   \item{Ecases}{Expected number of cases.} 
#'   \item{logEdeaths}{Log of the column \code{Edeaths}}
#'   \item{logEcases}{Log of the column \code{Ecases}}
#'   \item{casesmr}{Standaridised morbidity rate for the number of cases, 
#'   \code{noofcases}/\code{Ecases}}
#'   \item{nweek}{Number of weeks during March 13 to July 31, 2020. All values are 20.}
#'   \item{noofhighweeks}{Number of  weeks out of 20 when the \code{casesmr} 
#'      was greater than 1}      
#' }
# usethis::use_data(engtotals, overwrite = TRUE)
#' @references
#' \insertAllCited{}
#' @examples
#'  colnames(engtotals)
#'  dim(engtotals)
#'  summary(engtotals[, 5:14])
"engtotals"


#' Number of weekly Covid-19 deaths and cases in the 313 local 
#' Local Authority Districts, Counties  and Unitary Authorities (LADCUA) in 
#' England during the 20 peaks in the first peak from March 13 to July 31, 2020. 
#' @source \insertCite{SahuBohning2021;textual}{bmstdr}. 
#'  @format A data frame with 6260 rows and 24 columns:
#' \describe{
#'   \item{Areacode}{Areacode identifier of the 313 
#'   Local Authority Districts, Counties  and Unitary Authorities (LADCUA)}
#'   \item{mapid}{A numeric column identifying the map area needed for plotting}
#'   \item{spaceid}{A numeric variable taking value 1 to 313 identifying the LADCUA's}
#'   \item{Region}{Identifies one of the 9 English regions}
#'   \item{popn}{Population number in mid-2019}
#'   \item{jsa}{Percentage of the working age  population receiving job-seekers allowance 
#'    during January 2020}
#'   \item{houseprice}{Median house price  in March 2020}
#'   \item{popdensity}{Population density in mid-2019}    
#'   \item{no2}{Estimated average value of NO2 at the centroid of the LADCUA}
#'   \item{covid}{Number of Covid-19 deaths within 28 days of a positive test}
#'   \item{allcause}{Number deaths}
#'   \item{noofcases}{Number of cases}
#'   \item{n0}{Log of the standardized case morbidity during the current week}
#'   \item{n1}{Log of the standardized case morbidity during the week before}
#'   \item{n2}{Log of the standardized case morbidity during the second week before}
#'   \item{n3}{Log of the standardized case morbidity during the third week before}
#'   \item{n4}{Log of the standardized case morbidity during the fourth week before}
#'   \item{Edeaths}{Expected number of Covid-19 deaths. See Sahu and 
#'   Bohning (2021) for methodology. }
#'   \item{Ecases}{Expected number of cases.} 
#'   \item{logEdeaths}{Log of the column \code{Edeaths}}
#'   \item{logEcases}{Log of the column {Ecases}}
#'   \item{highdeathsmr}{A binary (0-1) random variable taking the value 1 if 
#'   the SMR of Covid-19 death is higher than 1}       
#' }
# usethis::use_data(engdeaths, overwrite = TRUE)
#' @references
#' \insertAllCited{}
#' @examples
#'  colnames(engdeaths)
#'  dim(engdeaths)
#'  summary(engdeaths[, 11:24])
"engdeaths"


#' A 313 by 313 proximity matrix for the 313 LADCUAS in England. Each entry is either 0 or 1 
#' and is 1 if the corresponding row and column LADCUAs share a common boundary.  
#' @source \insertCite{SahuBohning2021;textual}{bmstdr}. 
# usethis::use_data(Weng, overwrite = TRUE)
#' @references
#' \insertAllCited{}
#' @examples
#'  dim(Weng)
#'  summary(apply(Weng, 1, sum))
"Weng"


