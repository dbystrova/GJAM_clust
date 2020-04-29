\name{gjamOrdination}
\alias{gjamOrdination}
\title{Ordinate gjam data}
\description{
  Ordinate data from a gjam object using correlation corresponding to reponse matrix \code{E}.
}
\usage{
  gjamOrdination(output, specLabs = NULL, col = NULL, cex = 1, 
                 PLOT=T, method = 'PCA')
}
\arguments{
  \item{output}{object of \code{\link{class}} "gjam".}

  \item{specLabs}{\code{character vector} of variable names in \code{colnames(output$y)}.}
  
  \item{col}{\code{character vector} of columns in \code{output$y} to label in plots.}
  
  \item{cex}{text size in plot.}
  
  \item{PLOT}{\code{logical}, if true, draw plots.}
  
  \item{method}{\code{character} variable can specify \code{'NMDS'}.}
}
\details{
 Ordinates the response correlation \code{ematrix} contained in \code{output$parameterTables}.  If \code{method = 'PCA'} returns eigenvalues and eigenvectors.  If \code{method = 'PCA'} returns three NMDS dimensions.  If \code{PLOT}, then plots will be generated.  Uses principle components analysis or non-metric multidimensional scale (NMDS).
}
\value{

 \item{\code{eVecs}}{ \code{S x S} or, if there is an \code{other} response variable to be excluded,  \code{S-1 x S-1} matrix of eigenvectors for species (rows) by eigenvectors (columns).}
 
 \item{\code{eValues}}{If \code{method = 'PCA'} returns length-\code{S} or, there is an \code{other} response variable to be excluded,  length-\code{S-1} vector of eigenvalues. If \code{method = 'NMDS'} this variable is \code{NULL}.}
 
}
\examples{
\dontrun{
f      <- gjamSimData(S = 30, typeNames = 'CA') 
ml     <- list(ng = 30, burnin = 5, typeNames = f$typeNames, holdoutN = 10)
output <- gjam(f$formula, f$xdata, f$ydata, modelList = ml)
ePCA   <- gjamOrdination(output, PLOT=FALSE)
eNMDS  <- gjamOrdination(output, PLOT=FALSE, method='NMDS')
}
}

\author{
James S Clark, \email{jimclark@duke.edu}
}
\references{
Clark, J.S., D. Nemergut, B. Seyednasrollah, P. Turner, and S. Zhang. 2017.  Generalized joint attribute modeling for biodiversity analysis: Median-zero, multivariate, multifarious data.  Ecological Monographs 87, 34-56.
}
\seealso{
\code{\link{gjam}} fits the data

A more detailed vignette is can be obtained with:

\code{browseVignettes('gjam')}

website 'http://sites.nicholas.duke.edu/clarklab/code/'.
}

