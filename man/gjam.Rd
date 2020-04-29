\name{gjam}
\alias{gjam}
\alias{print.gjam}
\alias{summary.gjam}
\title{Gibbs sampler for gjam data}
\description{
  Analyzes joint attribute data (e.g., species abundance) with Gibbs sampling.  
  Input can be output from \code{\link{gjamSimData}}.  Returns a list of objects from Gibbs sampling that can be plotted by \code{\link{gjamPlot}}.
}
\usage{
  gjam(formula, xdata, ydata, modelList)
  
  \method{print}{gjam}(x, ...)
  
  \method{summary}{gjam}(object, ...)
}
\arguments{
  \item{formula}{R formula for model, e.g., \code{~ x1 + x2}.}

  \item{xdata}{\code{data.frame} containing predictors in \code{formula}.  If not found in \code{xdata} variables, they must be available from the user's workspace.}

  \item{ydata}{\code{n} by \code{S} response \code{matrix} or \code{data.frame}.  Column names are unique labels, e.g., species names.  All columns will be included in the analysis.}

  \item{modelList}{\code{list} specifying inputs, including \code{ng} (number of Gibbs steps), \code{burnin}, and \code{typeNames}. Can include the number of holdouts for out-of-sample prediction, \code{holdoutN}.  See \bold{Details}.}
  
  \item{x}{object of \code{class \link{gjam}}.}
  
  \item{object}{currently, also an object of \code{class \link{gjam}}.}
  
  \item{...}{further arguments not used here.}
}

\details{


Note that \code{formula} begins with \code{~}, not \code{y ~}.  The response \code{matrix} is passed in the form of a \code{n} by \code{S} matrix or \code{data.frame} \code{ydata}.  
  
Both predictors in \code{xdata} and responses in \code{ydata} can include missing values as \code{NA}. Factors in \code{xdata} should be declared using \code{factor}.   For computational stability variables that are not factors are standardized by mean and variance, then transformed back to original scales in \code{output}.  To retain a variable in its original scale during computation include it in the character string \code{notStandard} as part of the \code{list modelList}. (example shown in the \code{vignette} on traits).  
  
  
\code{modelList} has these defaults and provides these options:
  
  \code{ng = 2000}, number of Gibbs steps.  
    
  \code{burnin = 500}, no. initial steps, must be less than \code{ng}.  
    
  \code{typeNames} can be \code{'PA'} (presenceAbsence), \code{'CON'} (continuous on \code{(-Inf, Inf)}), \code{'CA'} (continuous abundance, zero censoring), \code{'DA'} (discrete abundance),  \code{'FC'} (fractional composition),
  \code{'CC'} (count composition), \code{'OC'} (ordinal counts), \code{'CAT'} (categorical classes).  \code{typeNames} can be a single value that applies to all columns in \code{ydata}, or there can be one value for each column.
    
  \code{holdoutN = 0}, number of observations to hold out for out-of-sample
    prediction.  
    
  \code{holdoutIndex = numeric(0)}, \code{numeric vector} of observations (row numbers) to holdout for out-of-sample prediction.  
    
  \code{censor = NULL}, \code{list} specifying columns, values, and intervals for
    censoring, see \code{\link{gjamCensorY}}.
    
  \code{effort = NULL}, \code{list} containing \code{'columns'}, a vector of length \code{<= S} giving the names of columns in in \code{y}, and \code{'values'}, a length-\code{n} vector of effort or a \code{n} by \code{S} matrix (see Examples).  \code{effort} can be plot area, search time, etc. for discrete count data \code{'DA'}.
  
  \code{FULL = F} in \code{modelList} will save full prediction chains in \code{$chains$ygibbs}.
     
  \code{notStandard = NULL}, \code{character vector} of column names in \code{xdata} that should not be standardized.
  
  \code{reductList = list(N = 20, r = 3)}, \code{list} of dimension reduction parameters, invoked when \code{reductList} is included in \code{modelList} or automatically when \code{ydata} has too many columns.  See \code{vignette} on Dimension Reduction.
  
  \code{random}, \code{character} string giving the name of a column in \code{xdata} that will be used to specify random effects.  The random group column should be declared as a \code{factor}.  There should be replication, i.e., each group level occurs multiple times.
  
  \code{REDUCT = F} in \code{modelList} overrides automatic dimension reduction. 
  
  \code{FCgroups, CCgroups}, are length-\code{S vectors} assigned to columns in \code{ydata} indicating composition \code{'FC'} or \code{'CC'} group membership.  For example, if there are two 'CA' columns in \code{ydata} followed by two groups of fractional composition data, each having three columns, then \code{typeNames = c('CA','CA','FC','FC','FC','FC','FC','FC')} and \code{FCgroups = c(0,0,1,1,1,2,2,2)}.  note: \code{gjamSimData} is not currently set up to simulate multiple composition groups, but \code{gjam} will model it.
  
  \code{PREDICTX = T} executes inverse prediction of \code{x}.  Speed-up by setting \code{PREDICTX = F}.
  
    \code{ematAlpha = .5} is the probability assigned for conditional and marginal independence in the \code{ematrix}.
  
  \code{traitList = list(plotByTrait, traitTypes, specByTrait)}, list of trait objects.  See vignette on Trait analysis.
  
  
More detailed vignettes can be obtained with:

\code{browseVignettes('gjam')}
}

\value{

Returns an object of \code{\link{class} "gjam"}, which is a list containing the following components:

\item{call}{function call}

\item{chains}{\code{list} of MCMC matrices, each with \code{ng} rows; includes coefficients \code{bgibbs}(\code{Q*S} columns), \code{bgibbsUn} (unstandardized for \code{x}), sensitivity \code{fgibbs} (\code{Q1} columns), and \code{fbgibbs} (\code{Q1} columns, where \code{Q1 = Q - 1}, unless there are multilevel factors); covariance \code{sgibbs} has \code{S*(S + 1)/2} columns (\code{REDUCT == F}) or \code{N*r} columns (\code{REDUCT == T}).}

\item{fit}{\code{list} of diagnostics (\code{DIC, rmspeAll, rmspeBySpec, xscore, yscore}).}

\item{inputs}{\code{list} of input summaries, including \code{breakMat} (partition matrix), \code{classBySpec} (interval assignment),  \code{designTable} (summary of design matrix), [\code{factorBeta, interBeta, intMat, linFactor}] (factor and interaction information), \code{other, notOther} (response columns to exclude and not), [\code{standMat, standRows, standX}] means and variances to standardize \code{x}, [\code{x, xdata, y}] cleaned versions of data.}

\item{missing}{\code{list} of missing objects, including locations for predictors \code{xmiss} and responses \code{ymiss} in \code{xdata} and \code{ydata}, respectively, predictor means \code{xmissMu} and standard errors \code{xmissSe}, response means \code{ymissMu} and standard errors \code{ymissSe} .}

\item{modelList}{\code{list} of model specifications from input \code{modelList}.}

\item{parameters}{\code{list} of parameter estimates, including coefficient matrices on standardized (\code{betaMu, betaSe}), unstandardized  (\code{betaMuUn, betaSeUn}), and dimensionless (\code{fBetaMu, fBetaSd}) scales; correlation (\code{corMu, corSe}) and covariance (\code{sigMu, sigSe}) matrices;  sensitivities to predictors (\code{fmatrix, fMu, fSe}); environmental response matrix (\code{ematrix}), with locations of zero elements, conditionally (\code{whConZero}) and marginally (\code{whichZero}), set at probability level \code{modelList$ematAlpha}); and latent variables (\code{wMu, wSd}).}

\item{prediction}{\code{list} of predicted values, including species \code{richness} (responses predicted > 0); inverse predicted \code{x} (\code{xpredMu, xpredSd}) and predicted \code{y} (\code{ypredMu, ypredSd}) matrices. 
    }

  
If \bold{traits} are modeled, then \code{parameters} will additionally include \code{betaTraitMu}, \code{betaTraitSe} (coefficients), \code{sigmaTraitMu}, \code{sigmaTraitSe} (covariance).  \code{prediction} will additionally include \code{tMuOrd} (ordinal trait means), \code{tMu}, \code{tSe} (trait predictions).


}
\examples{
\dontrun{
## combinations of scales
types <- c('DA','DA','OC','OC','OC','OC','CC','CC','CC','CC','CC','CA','CA','PA','PA')         
f    <- gjamSimData(S = length(types), typeNames = types)
ml   <- list(ng = 50, burnin = 5, typeNames = f$typeNames)
out  <- gjam(f$formula, f$xdata, f$ydata, modelList = ml)
summary(out)

# repeat with ng = 5000, burnin = 500, then plot data:
pl  <- list(trueValues = f$trueValues, SMALLPLOTS=F)
gjamPlot(out, plotPars = pl)

## discrete abundance with heterogeneous effort 
S   <- 5                             
n   <- 1000
eff <- list( columns = 1:S, values = round(runif(n,.5,5),1) )
f   <- gjamSimData(n, S, typeNames='DA', effort=eff)
ml  <- list(ng = 500, burnin = 50, typeNames = f$typeNames, effort = eff)
out <- gjam(f$formula, f$xdata, f$ydata, modelList = ml)
summary(out)

# repeat with ng = 2000, burnin = 500, then plot data:
pl  <- list(trueValues = f$trueValues)
gjamPlot(out, plotPars = pl)
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

