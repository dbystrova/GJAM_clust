\name{gjamPriorTemplate}
\alias{gjamPriorTemplate}
\title{Prior coefficients for gjam analysis}
\description{
  Constructs coefficient matrices for low and high limits on the uniform prior distribution for \code{beta}.
}
\usage{
  gjamPriorTemplate(formula, xdata, ydata, lo = NULL, hi = NULL)
}
\arguments{
  \item{formula}{object of class \code{formula}, starting with \code{~}, matches the \code{formula} passed to \code{gjam}}

  \item{xdata}{\code{n x Q} observation by predictor \code{data.frame}}

  \item{ydata}{\code{n x Q} observation by response \code{data.frame}}
  
  \item{lo}{\code{list} of lower limits}
  
  \item{hi}{\code{list} of upper limits}
}
\value{
A \code{list} containing two matrices.  \code{lo} is a \code{Q x S matrix} of lower coefficient limits.  \code{hi} is a \code{Q x S matrix} of upper coefficient limits.  Unless specied in \code{lo}, all values in \code{lo = -Inf}.  Likewise, unless specied in \code{hi}, all values in \code{hiBeta = -Inf}.
}
\details{
The prior distribution for a coefficient \code{beta[q,s]} for predictor \code{q} and response \code{s}, is \code{dunif(lo[q,s], hi[q,s])}.  \code{gjamPriorTemplate} generates these matrices.  The default values are (\code{-Inf, Inf}), i.e., all values in \code{lo} equal to \code{-Inf} and \code{hi} equal to \code{Inf}.  These templates can be modified by changing specific values in \code{lo} and/or \code{hi}.  

Alternatively, desired lower limits can be passed as the \code{list lo}, assigned to names in \code{xdata} (same limit for all species in \code{ydata}), in \code{ydata} (same limit for all predictors in \code{xdata}), or both, separating names in \code{xdata} and \code{ydata} by \code{"_"}.  The same convention is used for upper limits in \code{hi}. 

These matrices are supplied in as list \code{betaPrior}, which is included in \code{modelList} passed to \code{gjam}.  See examples and \code{browseVignettes('gjam')}.

Note that the informative prior slows computation.
}
\author{
James S Clark, \email{jimclark@duke.edu}
}
\references{
Clark, J.S., D. Nemergut, B. Seyednasrollah, P. Turner, and S. Zhang. 2017.  Generalized joint attribute modeling for biodiversity analysis: Median-zero, multivariate, multifarious data.  Ecological Monographs 87, 34-56.
}
\seealso{
\code{\link{gjam}}
}
\examples{
\dontrun{
library(repmis)
source_data("https://github.com/jimclarkatduke/gjam/blob/master/forestTraits.RData?raw=True")

xdata       <- forestTraits$xdata
plotByTree  <- gjamReZero(forestTraits$treesDeZero) # re-zero
traitTypes  <- forestTraits$traitTypes
specByTrait <- forestTraits$specByTrait

tmp <- gjamSpec2Trait(pbys = plotByTree, sbyt = specByTrait, 
                      tTypes = traitTypes)
tTypes <- tmp$traitTypes
traity <- tmp$plotByCWM
censor <- tmp$censor

formula <- as.formula(~ temp + deficit)
lo <- list(temp_ring = 0, deficit_drought = 0) # positive combinations
b  <- gjamPriorTemplate(formula, xdata, ydata = traity, lo = lo, hi = hi)
bp <- list(loBeta = b$lo, hiBeta = b$hi)

ml <- list(ng=4000, burnin=1000, typeNames = tTypes, censor = censor, 
            betaPrior = bp)
out <- gjam(formula, xdata, ydata = traity, modelList = ml)

S   <- ncol(traity)
sc  <- rep('black',S)
sc[colnames(traity) %in% c('ring','drought')] <- 'blue'
pl  <- list(SMALLPLOTS=F, specColor=sc)           
gjamPlot(output = out, plotPars = pl)         
}
}
