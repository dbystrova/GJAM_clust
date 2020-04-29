\name{gjamPlot}
\alias{gjamPlot}
\title{Plot gjam analysis}
\description{
  Constructs plots of posterior distributions, predictive distributions, 
  and additional analysis from output of \code{\link{gjam}}.
}
\usage{
  gjamPlot(output, plotPars)
}
\arguments{
  \item{output}{object of \code{\link{class}} "gjam"}
  
  \item{plotPars}{\code{list} having default values described in \bold{Details}}
}
\value{
  Summary tables of parameter estimates are:
  
  \item{\code{betaEstimates}}{Posterior summary of beta coefficients.}
  
  \item{\code{clusterIndex}}{cluster index for responses in grid/cluster plots.}
    
  \item{\code{clusterOrder}}{order for responses in grid/cluster plots.}
  
  \item{\code{eComs}}{groups based on clustering \code{ematrix}.}

  \item{\code{ematrix}}{\code{S X S} response correlation matrix for E.}
  
  \item{\code{eValues}}{eigenvalues of \code{ematrix}.}
  
   \item{\code{eVecs}}{eigenvectors of \code{ematrix}.}
    
  \item{\code{fit}}{\code{list} containing DIC, score, and rmspe.}
  
}

\details{

  \code{plotPars} a \code{list} that can contain the following, listed with default values:
  
  \tabular{ll}{
  
    \code{PLOTY = T} \tab plot predicted \code{y}.\cr
  
    \code{PLOTX = T} \tab plot inverse predicted \code{x}.\cr
    
    \code{PREDICTX = T} \tab inverse prediction of \code{x}; does not work if \code{PREDICTX = F} in \code{link{gjam}}.\cr
  
    \code{ncluster} \tab number of clusters to highlight in 
    cluster diagrams, default based on \code{S}.\cr
  
    \code{CORLINES = T} \tab draw grid lines on grid plots of R and E.\cr
  
    \code{cex = 1} \tab text size for grid plots, see \code{\link{par}}.\cr
  
    \code{BETAGRID = T} \tab draw grid of beta coefficients.\cr
  
    \code{PLOTALLY = F} \tab an individual plot for each column in \code{y}.\cr
    
    \code{SMALLPLOTS = T} \tab avoids plot margin error on some devices, better appearance if \code{FALSE}.\cr
    
    \code{GRIDPLOTS = F} \tab cluster and grid plots derived from parameters; matrices R and E are discussed in Clark et al. (2016).\cr
    
    \code{SAVEPLOTS = F} \tab plots saved in pdf format.\cr
    
    \code{outfolder = 'gjamOutput'} \tab folder for plot files if \code{SAVEPLOTS = T}.\cr
  
    \code{width, height = 4} \tab can be small values, in inches, to avoid 
    plot margin error on some devices.\cr
    
    \code{specColor = 'black'} \tab color for posterior box-and-whisker plots.\cr
    
    \code{ematAlpha = .95} \tab prob threshold used to infer that a covariance value in \code{Emat} is not zero.\cr
    
    \code{ncluster = 4} \tab number of clusters to identify in \code{ematrix}.\cr
    
  }
  
  The 'plot margin' errors mentioned above are device-dependent.  They can be avoided by specifying small \code{width, height} (in inches) and by omitting the grid plots (\code{GRIDPLOTS = F}).  If plotting does not produce a 'plot margin error', better appearance is obtained with \code{SMALLPLOTS = F}.  
  
  
Names will not be legible for large numbers of species.  Specify \code{specLabs = F} and use a character vector for \code{specColor} to identify species groups (see the gjam vignette on dimension reduction).
  
Box and whisker plots bound 0.68 and 0.95 credible and predictive intervals. 
  
}
\examples{
\dontrun{
## ordinal data
f   <- gjamSimData(S = 15, Q = 3, typeNames = 'OC') 
ml  <- list(ng = 200, burnin = 50, typeNames = f$typeNames, holdoutN = 10)
out <- gjam(f$formula, f$xdata, f$ydata, modelList = ml)

# repeat with ng = 2000, burnin = 500, then plot data here:
pl  <- list(trueValues = f$trueValues, width=3, height=2)
fit <- gjamPlot(output = out, plotPars = pl)
}
}
\author{
James S Clark, \email{jimclark@duke.edu}
}
\references{
Clark, J.S., D. Nemergut, B. Seyednasrollah, P. Turner, and S. Zhang. 2017.  Generalized joint attribute modeling for biodiversity analysis: Median-zero, multivariate, multifarious data.  Ecological Monographs 87, 34-56.
}
\seealso{
\code{\link{gjam}}
A more detailed vignette is can be obtained with:

\code{browseVignettes('gjam')}

website 'http://sites.nicholas.duke.edu/clarklab/code/'.
}

