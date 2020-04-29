\name{gjamTrimY}
\alias{gjamTrimY}
\title{Trim gjam response data}
\description{
  Returns a \code{list} that includes a subset of columns in \code{y}.  Rare species can be aggregated into a single class.
}
\usage{
  gjamTrimY(y, minObs = 2, maxCols = NULL, OTHER = TRUE)
}
\arguments{
  \item{y}{\code{n} by \code{S} numeric response \code{matrix}}

  \item{minObs}{minimum number of non-zero observations}

  \item{maxCols}{maximum number of response variables}

  \item{OTHER}{\code{logical} or \code{character} string.  If \code{OTHER = TRUE}, rare species are aggregated in a new column \code{'other'}.  A \code{character} vector contains the names of columns in \code{y} to be aggregated with rare species in the new column \code{'other'}.}
}
\value{
  Returns a \code{list} containing three elements.

  \item{\code{y}}{trimmed version of \code{y}.}

  \item{\code{colIndex}}{length-\code{S vector} of indices for new columns in \code{y}.}

  \item{\code{nobs}}{number of non-zero observations by column in \code{y}.}
}
\details{
  Data sets commonly have many responses that are mostly zeros, large numbers of rare species, even singletons.  Response matrix \code{y} can be trimmed to include only taxa having \code{> minObs} non-zero observations or to \code{<= maxCol} total columns.  The option \code{OTHER} is recommended for composition data ('CC', 'FC'), where the \code{'other'} column is taken as the reference class.  If there are unidentified species they might be included in this class.  [See \code{\link{gjamSimData}} for \code{typeName} codes].  
}
\examples{
\dontrun{
library(repmis)
source_data("https://github.com/jimclarkatduke/gjam/blob/master/forestTraits.RData?raw=True")

y   <- gjamReZero(fungEnd$yDeZero)     # re-zero data
dim(y)
y   <- gjamTrimY(y, minObs = 200)$y    # species in >= 200 observations
dim(y)
tail(colnames(y))    # last column is 'other'
}
}
\author{
James S Clark, \email{jimclark@duke.edu}
}
\references{
Clark, J.S., D. Nemergut, B. Seyednasrollah, P. Turner, and S. Zhang. 2017. 
Generalized joint attribute modeling for biodiversity analysis: Median-zero, multivariate, multifarious data. Ecological Monographs 87, 34-56.
}
\seealso{
  \code{\link{gjamSimData}} simulates data
  \code{\link{gjam}} analyzes data
  
  A more detailed vignette is can be obtained with:

\code{browseVignettes('gjam')}

web site 'http://sites.nicholas.duke.edu/clarklab/code/'.
}

