\name{gjamSpec2Trait}
\alias{gjamSpec2Trait}
\title{Ecological traits for gjam analysis}
\description{
  Constructs community-weighted mean-mode (CWMM) trait matrix for analysis with \code{\link{gjam}} for \code{n} observations, \code{S} species, \code{P} traits, and \code{M} total trait levels.
}
\usage{
  gjamSpec2Trait(pbys, sbyt, tTypes)
}
\arguments{
  \item{pbys}{\code{n x S} plot by species matrix (presence-absence, abundance)}

  \item{sbyt}{\code{S x P} species by trait matrix}

  \item{tTypes}{\code{P} data types for trait columns}
}
\value{
  \item{plotByCWM}{\code{n x M matrix} of community-weight means (numeric) or modes (ordinal)}
  
  \item{traitTypes}{\code{character vector} of data types for traits}
  
  \item{specByTrait}{\code{S x M matrix} translates species to traits}
  
  \item{censor}{\code{list} of censored columns, values, and intervals; see \code{\link{gjamCensorY}}}
}
\details{
  Generates the objects needed for a trait response model (TRM).  As inputs the \code{sbyt} \code{data.frame} has \code{P} columns containing numeric values, ordinal scores, and categorical variables, identified by data type in \code{tTypes}.  Additional trait columns can appear in the \code{n x M} output \code{matrix plotByCWMM}, because each level of a category becomes a new \code{'FC'} column as a CWMM.  Thus, \code{M} can exceed \code{P}, depending on the number of factors in \code{sbyt}. The exception is for categorical traits with only two levels, which can be treated as (0, 1) censored \code{'CA'} data.
  
  As output, the CWMM data types are given in \code{traitTypes}.  
  
  The list \code{censor = NULL} unless some data types are censored.  In the example below there are two censored columns.

  A detailed vignette on trait analysis is obtained with:

\code{browseVignettes('gjam')}

}
\author{
James S Clark, \email{jimclark@duke.edu}
}
\references{
Clark, J.S. 2016. Why species tell us more about traits than traits tell us about species: Predictive models. Ecology 97, 1979-1993.

Clark, J.S., D. Nemergut, B. Seyednasrollah, P. Turner, and S. Zhang. 2017.  Generalized joint attribute modeling for biodiversity analysis: Median-zero, multivariate, multifarious data.  Ecological Monographs 87, 34-56.
}
\seealso{
\code{\link{gjam}}, \code{\link{gjamCensorY}}
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

ml  <- list(ng=2000, burnin=500, typeNames = tTypes, censor = censor)
out <- gjam(~ temp + stdage + deficit, xdata, ydata = traity, modelList = ml)
pl  <- list(SMALLPLOTS=F)           
gjamPlot(output = out, plotPars = pl)         
}
}
