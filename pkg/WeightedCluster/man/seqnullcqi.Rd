\name{seqnullcqi}
\alias{seqnullcqi}
\alias{plot.seqnullcqi}
\alias{print.seqnullcqi}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	 Sequence Analysis Typologies Validation Using Parametric Bootstrap 
}
\description{
\code{seqnullcqi} implements the methodology proposed by Studer (2021) for the validation of sequence analysis typologies using  parametric bootstraps. The method works by comparing the cluster quality of an observed typology with the quality obtained by clustering similar but nonclustered data. Several models to test the different structuring aspects of the sequences important in life-course research, namely, sequencing, timing, and duration (see function \code{\link{seqnull}}). This strategy allows identifying the key structural aspects captured by the observed typology. Plot and print methods of the \code{seqnullcqi} results are also provide. 
}
\usage{
seqnullcqi(seqdata, clustrange, R, model=c("combined", "duration", "sequencing", "stateindep", "Markov", "userpos"),
           seqdist.args=list(), kmedoid = FALSE, hclust.method="ward.D", ...)
\method{plot}{seqnullcqi}(x, stat, type = c("line", "density", "boxplot", "seqdplot"), 
quant = 0.95, norm = TRUE, legendpos = "topright", alpha = 0.2, ...)
\method{print}{seqnullcqi}(x, stat="all", digits=2, ...) 
} 
%- maybe also 'usage' for other objects documented here.
\arguments{
	\item{seqdata}{State sequence object of class \code{stslist}. The sequence data to use. Use \code{\link{seqdef}} to create such an object.}
	\item{clustrange}{The clustering of the data to be validated as an object of class \code{clustrange}.	See \code{\link{as.clustrange}} or \code{\link{wcKMedRange}} to create such an object.}
	\item{model}{String. The model used to generate the similar but nonclustered data. It can be one of \code{"combined"}, \code{"duration"}, \code{"sequencing"}, \code{"stateindep"}, \code{"Markov"} or \code{"userpos"}. See \code{\link{seqnull}} for more information.}
	\item{R}{The number of bootstraps.}
	\item{seqdist.args}{List of arguments passed to \code{\link{seqdist}} for computing the distances.}
	\item{kmedoid}{Logical. If \code{TRUE}, the PAM algorithm is used to cluster the data using  \code{\link{wcKMedRange}}. If \code{FALSE}, \code{hclust} is used.}
	\item{hclust.method}{String. Hierarchical method to use with \code{\link{hclust}}.}
	\item{x}{A \code{\link{seqnullcqi}} object to be plotted or printed.}
	\item{stat}{Character. The statistic to plot or "all" for all statistics. See \code{\link{wcClusterQuality}} for a list of possible values.}
	\item{type}{Character. The type of graphic to be plotted. If \code{type="line"} (default), a transparent line representing the cluster quality index for each bootstrap is plotted using a separate line. If \code{type="density"}, the density of the maximum cluster quality index values among the different number of groups is plotted as well as the original cluster quality values. If \code{type="beanplot"}, beanplot of the distribution of the cluster quality index values for each number of groups is plotted separately. If \code{type="seqdplot"}, a state distribution sequence plot of the sequences generated with the null model is plotted (see \code{\link{seqdplot}}).}
	\item{quant}{Numeric. Quantile to use for the confidence intervals.}
	\item{norm}{Logical. If \code{TRUE}, cluster quality indices are standardized using the mean and standard deviation of the null distribution.}
	\item{legendpos}{Character. legend position, see \code{\link{legend}}.}
	\item{alpha}{Transparency parameter for the lines to be drawn (only for \code{type="line"}).}
	\item{digits}{Transparency parameter for the lines to be drawn (only for \code{type="line"}).}
	\item{\dots}{Additionnal parameters passed to \code{\link{seqnull} (for \code{seqnullcqi}) or \code{\link{plot}} or \code{\link{print}}}.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
