\name{gjamPoints2Grid}
\alias{gjamPoints2Grid}
\title{Incidence point pattern to grid counts}
\description{
  From point pattern data in (x, y) generates counts on a lattice supplied by the user or specified by lattice size or density.  For analysis in \code{\link{gjam} as counts (known effort) or count composition (unknown effort) data}.
}
\usage{
  gjamPoints2Grid(specs, xy, nxy = NULL, dxy = NULL, 
                  predGrid = NULL, effortOnly = TRUE)
}
\arguments{
  \item{specs}{\code{character vector} of species names or codes. }
  
  \item{xy}{\code{matrix} with rows = \code{length(specs)} and columns for (x, y).}
  
  \item{nxy}{length-2 \code{numeric vector} with numbers of points evenly spaced on (x, y).}
  
  \item{dxy}{length-2 \code{numeric vector} with distances for points evenly spaced on (x, y).}
  
  \item{predGrid}{\code{matrix} with 2 columns for (x, y).}
        
  \item{effortOnly}{\code{logical} to return only points where counts are positive (e.g., effort is unknown).}
}
\value{
  
  \item{\code{gridBySpec}}{\code{matrix} with rows for grid locations, columns for counts by species.}
  
  \item{\code{predGrid}}{\code{matrix} with columns for (x, y) and rows matching \code{gridBySpec}.}
  
}

\details{
  
For incidence data with species names \code{specs} and locations \code{(x, y)} constructs a lattice based a prediction grid \code{predGrid}, at a density of \code{(dx, dy)}, or with numbers of lattice points \code{(nx, ny)}.  If \code{effortOnly = T}, returns only points with non-zero values.

A prediction grid \code{predGrid} would be passed when counts by locations of known effort are required or where multiple groups should be assign to the same lattice points.

The returned \code{gridBySpec} can be analyzed in \code{\link{gjam}} with known effort as count data \code{"DA"} or with unknown effort as count composition data \code{"CC"}.
  
}
\examples{
\dontrun{
## random data
n  <- 100
s  <- sample( letters[1:3], n, replace = TRUE)
xy <- cbind( rnorm(n,0,.2), rnorm(n,10,2) )

nx <- ny <- 5                                    # uniform 5 X 5 lattice
f  <- gjamPoints2Grid(s, xy, nxy = c(nx, ny))
plot(f$predGrid[,1], f$predGrid[,2], cex=.1, xlim=c(-1,1), ylim=c(0,20),
     xlab = 'x', ylab = 'y')
text(f$predGrid[,1], f$predGrid[,2], rowSums(f$gridBySpec))

dx <- .2                                          # uniform density
dy <- 1.5
g  <- gjamPoints2Grid(s, xy, dxy = c(dx, dy))
text(g$predGrid[,1], g$predGrid[,2], rowSums(g$gridBySpec), col='brown')

p  <- cbind( runif(30, -1, 1), runif(30, 0, 20) ) # irregular lattice
h  <- gjamPoints2Grid(s, xy, predGrid = p)
text(h$predGrid[,1], h$predGrid[,2], rowSums(h$gridBySpec), col='blue')
}
}
\author{
James S Clark, \email{jimclark@duke.edu}
}
\references{
Clark, J.S., D. Nemergut, B. Seyednasrollah, P. Turner, and S. Zhang. 2016.  Generalized joint attribute modeling for biodiversity analysis: Median-zero, multivariate, multifarious data.  Ecological Monographs 87, 34-56.
}
\seealso{
\code{\link{gjam}}
A more detailed vignette is can be obtained with:

\code{browseVignettes('gjam')}

}

