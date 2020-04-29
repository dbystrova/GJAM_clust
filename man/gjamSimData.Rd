\name{gjamSimData}
\alias{gjamSimData}
\title{Simulated data for gjam analysis}
\description{
  Simulates data for analysis by \code{\link{gjam}}.
}
\usage{
  gjamSimData(n = 1000, S = 10, Q = 5, x = NULL, nmiss = 0, typeNames, effort = NULL)
}
\arguments{
  \item{n}{Sample size}

  \item{S}{Number of response variables (columns) in \code{y}, typically less than
   \code{n}}

  \item{Q}{Number of predictors (columns) in design matrix \code{x << n}}
  
  \item{x}{design \code{matrix}, if supplied \code{n} and \code{Q} will be set to \code{nrow(x)} and \code{ncol(x)}, respectively}

  \item{nmiss}{Number of missing values to in \code{x << n}}

  \item{typeNames}{Character vector of data types, see \bold{Details}}

  \item{effort}{List containing '\code{columns}' specifying columns to which
  effort applies, and '\code{values}', a length-\code{n} vector of effort per 
  observation.}
}
\value{

  \item{formula}{R formula for model, e.g., \code{~ x1 + x2}}
  
  \item{xdata}{\code{data.frame} includes columns for predictors in the
  design matrix}
  
    \item{ydata}{\code{data.frame} for the simulated response}
    
    \item{y}{response as a  \code{n} by \code{S} \code{matrix} as assembled in \code{\link{gjam}}.}
  
    \item{w}{\code{n} by \code{S} latent states}
  
  \item{typeY}{vector of data types corresponding to columns in \code{y}, see \bold{Details}} 
  
 \item{typeNames}{vector of data types corresponding to columns in \code{ydata}} 
  
  \item{trueValues}{list containing true parameter values \code{beta} (regression
  coefficients), \code{sigma} (covariance matrix), \code{corSpec} (correlation 
  matrix corresponding to \code{sigma}), and \code{cuts} (partition matrix for 
  ordinal data).}
  
  \item{effort}{ see \bold{Arguments}.}
  
}
\details{
  Generates simulated data and parameters for analysis by \code{\link{gjam}}.  Because both parameters and data are stochastic, not all simulations will give good results.   
  
  \code{typeNames} can be '\code{PA}' (presenceAbsence), '\code{CA}'  
  (continuous), '\code{DA}' (discrete),  '\code{FC}' (fractional composition),
  '\code{CC}' (count composition), '\code{OC}' (ordinal counts), and '\code{CAT}' (categorical levels).  If more than one \code{'CAT'} is included, each defines a multilevel categorical reponse.
  One additional type, '\code{CON}' (continuous), is not censored at zero by default.  
  
  If defined as a single \code{character} value \code{typeNames} applies to all columns in \code{y}.  If not, \code{typeNames} is length-\code{S} \code{character vector}, identifying each response by column in \code{y}.  If a column \code{'CAT'} is included, a random number of levels will be generated, \code{a, b, c, ...}.
  
  A more detailed vignette is can be obtained with:

\code{browseVignettes('gjam')}

website 'http://sites.nicholas.duke.edu/clarklab/code/'.

}
\author{
James S Clark, \email{jimclark@duke.edu}
}
\references{
Clark, J.S., D. Nemergut, B. Seyednasrollah, P. Turner, and S. Zhang. 2016.  Generalized joint attribute modeling for biodiversity analysis: Median-zero, multivariate, multifarious data.  Ecological Monographs 87, 34-56.
}
\seealso{
\code{\link{gjam}}
}
\examples{
\dontrun{
## ordinal data, show true parameter values
sim <- gjamSimData(S = 5, typeNames = 'OC')  
sim$ydata[1:5,]                              # example data
sim$trueValues$cuts                          # simulated partition
sim$trueValues$beta                          # coefficient matrix

## continuous data censored at zero, note latent w for obs y = 0
sim <- gjamSimData(n = 5, S = 5, typeNames = 'CA')  
sim$w
sim$y

## continuous and discrete data
types <- c(rep('DA',5), rep('CA',4))
sim   <- gjamSimData(n = 10, S = length(types), Q = 4, typeNames = types)
sim$typeNames
sim$ydata
                             
## composition count data  
sim <- gjamSimData(n = 10, S = 8, typeNames = 'CC')
totalCount <- rowSums(sim$ydata)
cbind(sim$ydata, totalCount)  # data with sample effort

## multiple categorical responses - compare matrix y and data.frqme ydata
types <- rep('CAT',2)
sim   <- gjamSimData(S = length(types), typeNames = types)
head(sim$ydata)
head(sim$y)

## discrete abundance, heterogeneous effort 
S   <- 5
n   <- 1000
ef  <- list( columns = 1:S, values = round(runif(n,.5,5),1) )
sim <- gjamSimData(n, S, typeNames = 'DA', effort = ef)
sim$effort$values[1:20]

## combinations of scales, partition only for 'OC' columns
types <- c('OC','OC','OC','CC','CC','CC','CC','CC','CA','CA','PA','PA')
sim   <- gjamSimData(S = length(types), typeNames = types)
sim$typeNames                           
head(sim$ydata)
sim$trueValues$cuts
}
}
