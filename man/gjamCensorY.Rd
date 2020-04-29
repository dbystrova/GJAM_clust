\name{gjamCensorY}
\alias{gjamCensorY}
\title{Censor gjam response data}
\description{
  Returns a \code{list} with censored values, intervals, and censored response \code{matrix y}.
}
\usage{
  gjamCensorY(values, intervals, y, type='CA', whichcol = c(1:ncol(y)))
}
\arguments{
  \item{values}{Values in \code{y} that are censored, specified by \code{intervals}}

  \item{intervals}{\code{matrix} having two rows and one column for each value 
  in \code{values}. The first row holds lower bounds.  The second row holds
  upper bounds.  See \bold{Examples}.}

  \item{y}{Response \code{matrix}, \code{n} rows by \code{S} columns.  All values within \code{intervals} will be replaced with \code{values}}
  
  \item{type}{Response type, see \code{typeNames} in \code{\link{gjam}}}

  \item{whichcol}{Columns in \code{y} that are censored (often not all responses are
  censored)}
}
\value{
  Returns a list containing two elements.

  \item{\code{y}}{n by S matrix updated with censored {values} substituted for
  those falling within \code{intervals}.}

  \item{\code{censor}}{\code{list} containing \code{$columns} that are censored and
  \code{$partition}, a matrix with 3 rows used in \code{\link{gjam}} 
  and \code{\link{gjamPlot}}, one column per censor interval.  Rows 
  are values, followed by lower and upper bounds.}
}
\details{
Any values in \code{y} that fall within censored \code{intervals} are replaced with censored \code{values}.  The example below simulates data collected on an 'octave scale': 0, 1, 2, 4, 8, ..., an approach to accelerate data collection with approximate bins.
}
\examples{
\dontrun{
# data in octaves
v  <- up <- c(0, 2^c(0:4), Inf)         
dn <- c(-Inf, v[-length(v)])
i  <- rbind( dn, up )  # intervals

f  <- gjamSimData(n = 2000, S = 15, Q = 3, typeNames='CA')
y  <- f$y
cc <- c(3:6)                                                # censored columns
g  <- gjamCensorY(values = v, intervals = i, y = y, whichcol = cc)
y[,cc] <- g$y                                               # replace columns
ml     <- list(ng = 50, burnin = 10, censor = g$censor, typeNames = f$typeNames)
output <- gjam(f$formula, xdata = f$xdata, ydata = y, modelList = ml)

#repeat with ng = 2000, burnin = 500, then:
pl  <- list(trueValues = f$trueValues, width = 3, height = 3)   
gjamPlot(output, pl)

# upper detection limit
up <- 5
v  <- up
i  <- matrix(c(up,Inf),2)
rownames(i) <- c('down','up')

f   <- gjamSimData(typeNames='CA')   
g   <- gjamCensorY(values = v, intervals = i, y = f$y)
ml  <- list(ng = 50, burnin = 10, censor = g$censor, typeNames = f$typeNames)
out <- gjam(f$formula, xdata = f$xdata, ydata = g$y, modelList = ml)

#repeat with ng = 2000, burnin = 500, then:
pl  <- list(trueValues = f$trueValues, width = 3, height = 3)           
gjamPlot(out, pl)

# lower detection limit
lo        <- .001
values    <- upper <- lo
intervals <- matrix(c(-Inf,lo),2)
rownames(intervals) <- c('lower','upper')
}
}
\author{
James S Clark, \email{jimclark@duke.edu}
}
\references{
Clark, J.S., D. Nemergut, B. Seyednasrollah, P. Turner, and S. Zhang. 2017.  Generalized joint attribute modeling for biodiversity analysis: Median-zero, multivariate, multifarious data.  Ecological Monographs 87, 34-56.
}
\seealso{
  \code{\link{gjamSimData}} simulates data
  \code{\link{gjam}} analyzes data
  
  A more detailed vignette is can be obtained with:

\code{browseVignettes('gjam')}

website 'http://sites.nicholas.duke.edu/clarklab/code/'.
}

