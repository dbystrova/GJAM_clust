\name{gjamFillMissingTimes}
\alias{gjamFillMissingTimes}
\title{Fill out data for time series (state-space) gjam}
\description{
  Fills in predictor, response, and effort matrics for time series data where there are multiple multivariate time series.  Time series gjam is still under development.
}
\usage{
  gjamFillMissingTimes(xdata, ydata, edata, groups, times, 
                       sequences=NULL, fillNA=T, fillTimes=T)
}
\arguments{
  \item{xdata}{\code{n by Q data.frame} holding predictor variables}
  
  \item{ydata}{\code{n by S matrix} holding response variables}
  
  \item{edata}{\code{n by S matrix} holding effort}
  
  \item{groups}{column name in \code{xdata} holding group indicator}
  
  \item{times}{column name in \code{xdata} holding time indicator}
  
  \item{sequences}{column name in \code{xdata} holding sequence indicator}
  
  \item{fillNA}{fill new rows in ydata with NA; otherwise interpolated value}
  
  \item{fillTimes}{insert rows for missing times: integers values in column \code{times}}
  

}
\details{
Fills missing times in \code{xdata, ydata, edata}.  If \code{fillTimes = TRUE} insert rows for missing times, including a \code{times = 0}, which serves as a prior mean for \code{times = 1}.  \code{groups} must be uniquely defined.  Within \code{groups} can be \code{sequences}, as when there is a time series for \code{groups} for each year.
}
\value{
A list containing the following:

  \item{xdata}{filled version of \code{xdata}}
  
  \item{ydata}{filled version of \code{ydata}}
  
  \item{edata}{filled version of \code{edata}}
  
  \item{timeZero}{row numbers in new data where each time series begins, with times = 0}
  
  \item{timeLast}{row numbers in new data where each time series ends}
  
  \item{rowInserts}{row numbers for all inserted rows}
  
  \item{noEffort}{rows for which effort is zero, except times = 0}
 
}
\examples{
\dontrun{
# under construction
}
}

\author{
James S Clark, \email{jimclark@duke.edu}
}
\references{
Clark, J.S., D. Nemergut, B. Seyednasrollah, P. Turner, and S. Zhang. 2016.  Generalized joint attribute modeling for biodiversity analysis: Median-zero, multivariate, multifarious data.  Ecological Monographs, 87, 34-56.
}
\seealso{
\code{\link{gjam}} for more on \code{xdata}, \code{ydata}, and effort.

A more detailed vignette is can be obtained with:

\code{browseVignettes('gjam')}

web site 'http://sites.nicholas.duke.edu/clarklab/code/'.
}

