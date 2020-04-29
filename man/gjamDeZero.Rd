\name{gjamDeZero}
\alias{gjamDeZero}
\title{Compress (de-zero) gjam data}
\description{
  Returns a de-zeroed (sparse matrix) version of \code{matrix ymat} with objects needed to re-zero it.
}
\usage{
  gjamDeZero(ymat)
}
\arguments{
  \item{ymat}{\code{n} by \code{S} response \code{matrix}}
}
\value{
  Returns a list containing the de-zeroed \code{ymat} as a vector \code{yvec}.

  \item{\code{yvec}}{non-zero elements of \code{ymat}}

  \item{\code{n}}{no. rows of \code{ymat} }
  
  \item{\code{S}}{no. cols of \code{ymat} }

  \item{\code{index}}{index for non-zeros}
  
  \item{\code{ynames}}{column names of \code{ymat}}
}
\details{
Many abundance data sets are mostly zeros.  \code{\link{gjamDeZero}} extacts non-zero elements for storage.  
}
\examples{
\dontrun{
library(repmis)
source_data("https://github.com/jimclarkatduke/gjam/blob/master/fungEnd.RData?raw=True")

ymat <- gjamReZero(fungEnd$yDeZero)  # OTUs stored without zeros
length(fungEnd$yDeZero$yvec)         # size of stored version
length(ymat)                         # full size
yDeZero <- gjamDeZero(ymat)
length(yDeZero$yvec)                 # recover de-zeroed vector
}
}
\author{
James S Clark, \email{jimclark@duke.edu}
}
\references{
Clark, J.S., D. Nemergut, B. Seyednasrollah, P. Turner, and S. Zhang. 2016.  Generalized joint attribute modeling for biodiversity analysis: Median-zero, multivariate, multifarious data.  Ecological Monographs 87, 34-56.
}
\seealso{
  \code{\link{gjamReZero}} to recover \code{ymat}
  
  \code{browseVignettes('gjam')}

website 'http://sites.nicholas.duke.edu/clarklab/code/'.
}

