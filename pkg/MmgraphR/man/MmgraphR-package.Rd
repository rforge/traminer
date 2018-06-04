\name{MmgraphR-package}
\alias{MmgraphR-package}
\alias{MmgraphR}
\docType{package}
\title{
Graphing for Markov, Hidden Markov, and Mixture Transition Distribution Models
}
\description{
MmgraphR produces parallel coordinate plots which map each element in a probability transition matrix as a line. Each line is weighted by probability, and so the thicker the line, the more probable a sequence of states is and vice-versa. Probability transition matrices reflecting higher order dependencies, often useful in life-course studies, can also be mapped. 

Additional functionalities of MmgraphR include probability filters that can easily be applied to emphasize the most (or least) probable state sequences overall, or by initial state. Any specific element(s) of the probability transition matrix can also be highlighted. To render the plot more attractive, various color palettes using the Hue-Chroma-Luminance color scheme have been built-in and can be easily selected by the user. MmgraphR is available through R-forge, and although can be used on its own, it is also compatible with the output of R packages \command{\link[msm]{msm}}, \command{\link[HiddenMarkov]{HiddenMarkov}}, \command{\link[HMM]{HMM}}, \command{\link[depmixS4]{depmixS4}}, \command{\link[seqHMM]{seqHMM}}, as well as \code{march}.
}
\details{
\tabular{ll}{
Package: \tab MmgraphR\cr
Version: \tab 0.2\cr
Date: \tab 2018-06-04\cr
Depends: \tab R (>= 3.1.1), TraMineR, colorspace\cr
License: \tab GPL (>= 2)\cr
LazyData: \tab true\cr
URL: \tab https://r-forge.r-project.org/R/?group_id=743\cr
Encoding: \tab latin1\cr
}
}
\author{
	Pauline (Poulcheria) Adamopoulou [aut, cre, cph],
	Andre Berchtold [ths],
  Gilbert Ritschard [ths],
	Reto Buergin [cbt],
Maintainer: Pauline (Poulcheria) Adamopoulou <padamopo@gmail.com>
}

\references{
	Buergin, R. and G. Ritschard (2014), "A decorated parallel coordinate plot for categorical longitudinal data", \emph{The American Statistician}. Vol. 68(2), pp. 98-103.

	Zeileis, A.; Hornik, K. and P. Murrell (2009), "Escaping RGBland: Selecting Colors for Statistical Graphics", \emph{Computational Statistics & Data Analysis}. Vol. 53, pp. 3259-3270.
}

\examples{

## Plotting a probability transition matrix
###########################################

trMat<-matrix( c (0.1, 0.05, 0.05, 0.80,
                  0.06, 0.02, 0.03, 0.89,
                  0.03, 0.01, 0.01, 0.95,
                  0, 0, 0, 1), nrow = 4, ncol = 4, byrow = TRUE)

trmatplot(trMat)

#--- Setting a seed so that the graphic can be replicated

trmatplot(trMat, seed = 2)

#--- Defining a second order probability transition matrix as an array

trMatArray <- ( array( c (
									0.30, 0.70,
                  0.65, 0.35, 
									0.05, 0.95,
                  0.99, 0.01), dim = c( 1, 2 , 4)) )

trmatplot (trMatArray, seed = 4, morder = 2, lwd = 0.5)

# cspal: ready-to-use color palettes using colorspace
#####################################################

#--- Color palette "dynamic"

trmatplot(trMat, seed = 2, cspal = "dynamic")

#--- Color palette "harmonic"

trmatplot(trMat, seed = 2, cspal = "harmonic")

## pfilter: Filtering out most (or least) probable sequences
############################################################

#--- The most probable sequence given a state

trmatplot(trMat, seed = 2, pfilter = "smax")

#--- The ten least probable sequences

trmatplot(trMat, seed = 2, pfilter = "tmin", num = 10 )
}
\seealso{
	\command{\link{trmatplot}},  
	\command{\link{trmatplot.default}},  
	\command{\link{trmatplot.depmix.fitted}}, 
  \command{\link{trmatplot.hmm}},
	\command{\link{trmatplot.array}},
	\command{\link{march.Dcmm.trmatplot}}
}

\keyword{ package }
