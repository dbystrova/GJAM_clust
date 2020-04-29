\name{gjamReZero}
\alias{gjamReZero}
\title{Expand (re-zero) gjam data}
\description{
  Returns a re-zeroed \code{matrix y} from the de-zeroed vector, a sparse matrix.
}
\usage{
  gjamReZero( yDeZero )
}
\arguments{
  \item{yDeZero}{\code{list} created by \code{\link{gjamReZero}} containing number of rows \code{n}, number of columns \code{S}, index for non-zeros \code{index}, the vector of non-zero values \code{yvec}, and the column names \code{ynames}.}
}
\value{
  \item{\code{ymat}}{re-zeroed \code{n} by \code{S} matrix.}
}
\details{
Many abundance data sets are mostly zeros.  \code{gjamReZero} recovers the full matrix from de-zeroed \code{list yDeZero} written by \code{\link{gjamDeZero}}
}
\examples{
\dontrun{
library(repmis)
source_data("https://github.com/jimclarkatduke/gjam/blob/master/fungEnd.RData?raw=True")
ymat <- gjamReZero(fungEnd$yDeZero)  # OTUs stored without zeros
length(fungEnd$yDeZero$yvec)         # size of stored version
length(ymat)                         # full size
}
}
\author{
James S Clark, \email{jimclark@duke.edu}
}
\references{
Clark, J.S., D. Nemergut, B. Seyednasrollah, P. Turner, and S. Zhang. 2016.  Generalized joint attribute modeling for biodiversity analysis: Median-zero, multivariate, multifarious data.  Ecological Monographs 87, 34-56.
}
\seealso{
  \code{\link{gjamDeZero}} to de-zero \code{ymat}
  
  \code{browseVignettes('gjam')}

website: 'http://sites.nicholas.duke.edu/clarklab/code/'.
}

