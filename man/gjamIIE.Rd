\name{gjamIIE}
\alias{gjamIIE}
\title{Indirect effects and interactions for gjam data}
\description{
  Evaluates direct, indirect, and interactions from a \code{\link{gjam}} object.  Returns a \code{list} of objects that can be plotted by \code{\link{gjamIIEplot}}.
}
\usage{
  gjamIIE(output, xvector, MEAN = T, keepNames = NULL, omitY = NULL, 
          sdScaleX = T, sdScaleY = F)
          
}
\arguments{
  \item{output}{object of \code{\link{class}} inheriting from "gjam".}

  \item{xvector}{vector of predictor values, with names, corresponding to columns in \code{output$x}.}

  \item{MEAN}{\code{logical}, if false, then median used.}
  
  \item{omitY}{\code{character} vector of columns in \code{output$y} to omit from calculations.}

  \item{keepNames}{\code{character} vector of columns in \code{output$y}. If omitted, all columns used.}
 
 \item{sdScaleX}{standardize coefficients to X scale.}
 
  \item{sdScaleY}{standardize coefficients to correlation scale.}
}
\details{
For plotting or recovering effects.  The \code{list fit$IIE} has matrices for main effects (\code{mainEffect}), interactions (\code{intEffect}), direct effects (\code{dirEffect}), indirect effects (\code{indEffectTo}), and standard deviations for each.  The direct effects are the sum of main effects and interactions.  The indirect effects include main effects and interactions that come through other species, determined by covariance matrix \code{sigma}.  

If \code{sdScaleX = T} effects are standandardized from the Y/X to Y scale.  This is the typical standardization for predictor variables.  If \code{sdScaleY = T} effects are given on the correlation scale.  If both are true effects are dimensionless.  See the gjam vignette on dimension reduction.
}
\value{
  A list of objects for plotting by \code{\link{gjamIIEplot}}.
}
\examples{
\dontrun{
sim <- gjamSimData(S = 12, Q = 5, typeNames = 'CA')
ml  <- list(ng = 50, burnin = 5, typeNames = sim$typeNames)
out <- gjam(sim$formula, sim$xdata, sim$ydata, modelList = ml)

xvector <- colMeans(out$inputs$x)  #predict at mean values for data
xvector[1] <- 1

fit <- gjamIIE(output = out, xvector)

gjamIIEplot(fit, response = 'S1', effectMu = c('main','ind'), 
            effectSd = c('main','ind'), legLoc = 'topleft')
}
}
\author{
James S Clark, \email{jimclark@duke.edu}
}
\references{
Clark, J.S., D. Nemergut, B. Seyednasrollah, P. Turner, and S. Zhang. 2016.  Generalized joint attribute modeling for biodiversity analysis: Median-zero, multivariate, multifarious data.  Ecological Monographs 87, 34-56.
}
\seealso{
\code{\link{gjamIIEplot}} plots output from \code{gjamIIE}

A more detailed vignette is can be obtained with:

\code{browseVignettes('gjam')}

web site 'http://sites.nicholas.duke.edu/clarklab/code/'.
}

