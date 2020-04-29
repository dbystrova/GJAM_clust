\name{gjamPredict}
\alias{gjamPredict}
\title{Predict gjam data}
\description{
  Predicts data from a gjam object, including conditional and out-of-sample prediction.
}
\usage{
  gjamPredict(output, newdata = NULL, y2plot = NULL, ylim = NULL, 
              FULL = FALSE)
}
\arguments{
  \item{output}{object of \code{\link{class}} "gjam".}

  \item{newdata}{a \code{list} of data for prediction, see \bold{Details}.}
  
  \item{y2plot}{\code{character} vector of columns in \code{output$y} to plot.}
  
  \item{ylim}{vector of lower and upper bounds for prediction plot}
  
  \item{FULL}{will return full chains for predictions as \code{output$ychains}}
}
\details{
 If \code{newdata} is not specified, the response is predicted from \code{xdata}  as an in-sample prediction. If \code{newdata} is specified, prediction is either conditional or out-of-sample. 
 
Conditional prediction on a new set of \code{y} values is done if \code{newdata} includes the \code{matrix ycondData}, which holds columns to condition on.  \code{ycondData} must be a \code{matrix} and have column names matching those in \code{y} that it will replace.  \code{ycondData} must have at least one column, but fewer than \code{ncol(y)} columns.  Columns not included in \code{ycondData} will be predicted conditionally.
 
Alternatively, the \code{list newdata} can include a new version of \code{xdata} for out-of-sample prediction.  The version of \code{xdata} passed in \code{newdata} has the columns with the same names and variable types as \code{xdata} passed to \code{\link{gjam}}.  Note that factor levels must also match those included when fitting the model.  All columns in \code{y} will be predicted out-of-sample.  

For count composition data the effort (total count) is 1000.

Because there is no out-of-sample effort for \code{'CC'} data, values are predicted on the [0, 1] scale.

See examples below.
}
\value{

 \item{\code{x}}{design matrix.} 
 
 \item{\code{sdList}}{\code{list} of predictive means and standard errors includes \code{yMu, yPe} (predictive mean, SE), \code{wMu, wSe} (mean latent states and SEs)}
 
 \item{\code{piList}}{predictive intervals, only generated if \code{length(y) < 10000}, includes \code{yLo, yHi} (0.025, 0.975) prediction interval, \code{wLo, wHi} (0.025, 0.975) for latent states } 
 
  \item{\code{prPresent}}{\code{n x S matrix} of probabilities of presence}
  
  \item{\code{ematrix}}{effort}
  
  \item{\code{ychains}}{full prediction chains if \code{FULL = T}}
 
}
\examples{
\dontrun{
S   <- 5
f   <- gjamSimData(n = 200, S = S, Q = 3, typeNames = 'CC') 
ml  <- list(ng = 50, burnin = 5, typeNames = f$typeNames, holdoutN = 10)
out <- gjam(f$formula, f$xdata, f$ydata, modelList = ml)

# predict data
par(mfrow=c(1,3),bty='n')             
gjamPredict(out, y2plot = colnames(f$ydata)) #predict the data in-sample
title('full sample')

# out-of-sample prediction
xdata     <- f$xdata[1:20,]
xdata[,3] <- mean(f$xdata[,3])     # mean for x[,3]
xdata[,2] <- seq(-2,2,length=20)   # gradient x[,2]
newdata   <- list(xdata = xdata, nsim = 50 )
p1 <- gjamPredict(out, newdata = newdata)

# plus/minus 1 prediction SE, default effort = 1000
x2   <- p1$x[,2]
ylim <- c(0, max(p1$sdList$yMu[,1] + p1$sdList$yPe[,1]))
plot(x2, p1$sdList$yMu[,1],type='l',lwd=2, ylim=ylim, xlab='x2',
     ylab = 'Predicted')
lines(x2, p1$sdList$yMu[,1] + p1$sdList$yPe[,1], lty=2)
lines(x2, p1$sdList$yMu[,1] - p1$sdList$yPe[,1], lty=2)

# .95 prediction error
lines(x2, p1$piList$yLo[,1], lty=3)
lines(x2, p1$piList$yHi[,1], lty=3)
title('SE and prediction, Sp 1')

# conditional prediction
ydataCond <- out$inputs$y[,1,drop=FALSE]*0          #set first column to zero
newdata   <- list(ydataCond = ydataCond, nsim=50)
p0        <- gjamPredict(output = out, newdata = newdata)

ydataCond <- ydataCond + 20                  #first column is 20
newdata   <- list(ydataCond = ydataCond, nsim=50)
p1        <- gjamPredict(output = out, newdata = newdata)

plot(out$inputs$y[,4],p0$sdList$yMu[,4], cex=.4,col='orange'); abline(0,1,lty=2)
points(out$inputs$y[,4],p1$sdList$yMu[,4], cex=.4,col='blue')
title('Cond. on 1st Sp')
}
}

\author{
James S Clark, \email{jimclark@duke.edu}
}
\references{
Clark, J.S., D. Nemergut, B. Seyednasrollah, P. Turner, and S. Zhang. 2016.  Generalized joint attribute modeling for biodiversity analysis: Median-zero, multivariate, multifarious data.  Ecological Monographs, 87, 34-56.
}
\seealso{
\code{\link{gjamSimData}} simulates data

A more detailed vignette is can be obtained with:

\code{browseVignettes('gjam')}

web site 'http://sites.nicholas.duke.edu/clarklab/code/'.
}

