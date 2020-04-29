\name{gjamSensitivity}
\alias{gjamSensitivity}
\title{Sensitivity coefficients for gjam}
\description{
  Evaluates sensitivity coefficients for full response matrix or subsets of it.
  Uses output from \code{\link{gjam}}.  Returns a \code{matrix} of samples by predictors.
}
\usage{

  gjamSensitivity(output, group=NULL, nsim=100)
  
}

\arguments{
  \item{output}{object fitted with \code{gjam}.}

  \item{group}{\code{character vector} of response-variable names from \code{output$inputs$y}.}

  \item{nsim}{number of samples from posterior distribution.}

}

\details{

Sensitivity to predictors of entire reponse matrix or a subset of it, identified by the \code{character string group}.  The equations for sensitivity are given here:


\code{browseVignettes('gjam')}
}

\value{

Returns a \code{nsim} by predictor matrix of sensitivities to predictor variables.

}
\examples{
\dontrun{
## combinations of scales
types <- c('DA','DA','OC','OC','OC','OC','CC','CC','CC','CC','CC','CA','CA','PA','PA')         
f    <- gjamSimData(S = length(types), typeNames = types)
ml   <- list(ng = 50, burnin = 5, typeNames = f$typeNames)
out  <- gjam(f$formula, f$xdata, f$ydata, modelList = ml)

ynames <- colnames(f$y)
group  <- ynames[types == 'OC']

full <- gjamSensitivity(out)
cc   <- gjamSensitivity(out, group)

nt <- ncol(full)

boxplot( full, boxwex = 0.25,  at = 1:nt - .21, col='blue', log='y',
         xaxt = 'n', xlab = 'Predictors', ylab='Sensitivity')
boxplot( cc, boxwex = 0.25, at = 1:nt + .2, col='forestgreen', add=T,
         xaxt = 'n')
axis(1,at=1:nt,labels=colnames(full))
legend('bottomright',c('full response','CC data'),
       text.col=c('blue','forestgreen'))
}
}
\author{
James S Clark, \email{jimclark@duke.edu}
}
\references{
Clark, J.S., D. Nemergut, B. Seyednasrollah, P. Turner, and S. Zhang. 2017.  Generalized joint attribute modeling for biodiversity analysis: Median-zero, multivariate, multifarious data.  Ecological Monographs, 87, 34-56.
}
\seealso{
\code{\link{gjamSimData}} simulates data

A more detailed vignette is can be obtained with:

\code{browseVignettes('gjam')}

website 'http://sites.nicholas.duke.edu/clarklab/code/'.
}

