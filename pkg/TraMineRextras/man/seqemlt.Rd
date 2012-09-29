\name{seqemlt}
\alias{seqemlt}
\title{Euclidean metric for longitudinal timelines}
\description{
Computes a Euclidean distance of sequences.
transforms sequences such as distance beween sequences is equivalent to Euclidean distance between transformed sequences. 
The transformed sequences may be used as inputs of any Euclidean algorithm (clustering algorithms, ...).
The distance is built considering the transitions between states at any step.
A weight of steps allows to balance short term/long term transitions.
The background - the duality between the sequence distances and the evolve of objects proximities - is analysed.  
}
\usage{
seqemlt(seqdata, a = 1, b = 1, weighted = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{seqdata}{ a state sequence object defined with the seqdef function.
%%     ~~Describe \code{seqdata} here~~
}
  \item{a}{optional argument for weighting steps in order the balance term/long term transitions. The weighting function is 1/(a*s+b) where s is the step of transition.
%%     ~~Describe \code{a} here~~
}
  \item{b}{see argument a.
%%     ~~Describe \code{b} here~~
}
  \item{weighted}{optional numerical vector containing weights, which may be used by some functions to compute weighted statistics (rates of transitions).
%%     ~~Describe \code{weighted} here~~
}
}
\details{The distance emlt between sequences is equivalent to Euclidean distance between the transformed sequences coord. Using coord as the data input of any clustering algorithm using a euclidean metric is equivalent of clustering with the emlt metric. 
A situation is defined as a state idexed with time, a sequence a timelines of states. The distance between situations is defined from the transitions between situations. The emlt distance between sequencies takes into account the proximity between situations. Transitions are considered at any steps with a weighting balance between long/short terms. A situation may have no occurency when the reffering object is not present during all the duration. The distance between any situation and a situation with no occurency is NA, and has no influence for the distance between sequencies.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  An object of class \code{emlt} with the following componant
  \item{coord}{transformed sequences. Euclidean metric emlt between sequences is equivalent to Euclidean distance between coord. coord is the input of any clustering algorithms using a Euclidean metric}
  \item{states }{list of states}
  \item{situations }{list of situations}
  \item{sit.freq }{frequence of situations}
  \item{sit.transrate }{rate of transitions from a situation to any situation of its own future : vector of transition towards future}
  \item{sit.profil }{profil of situations. The profil is a normalized vector issued from the rate of transition including a balance of short/long term with the weight of time 1/a*s+b, where s is the step of transition }
  \item{sit.cor }{Correlation between situations. Two situations are high correlated when their profiles are similar (ie their transitions towards future are similar). }
}
\references{- Rousset Patrick,  Giret Jean-françois,Classifying Qualitative Time Series with SOM: The Typology of Career Paths in France Lecture Notes in computer science, vol 4507, 2007, Springer Berlin / Heidelberg - 
- Rousset Patrick, Giret Jean-françois, Yvette Grelet (2012) Typologies De Parcours et Dynamique Longitudinale, Bulletin de méthodologie sociologique, issue 114, april 2012. - 
- Rousset Patrick, Giret Jean-françois (2008) A longitudinal Analysis of  Labour Market Data with SOM, Encyclopedia of Artificial Intelligence, Edition Information Science Reference -
%% ~put references to the literature/web site here ~
}
\author{Patrick Rousset, Senior researcher at Cereq, rousset@cereq.fr
}

\seealso{plot.emlt}
\examples{
library(TraMineR)

data(mvad)
mvad.seq <- seqdef(mvad[, 17:86])
alphabet(mvad.seq)
labels(mvad.seq)
## Computing distance
mvad.emlt <- seqemlt(mvad.seq)

## typology1 with kmeans in 8 clusters
km <- kmeans(mvad.emlt$coord, 8)

##Plotting typology1 by clusters
seqdplot(mvad.seq, group=km$cluster)

## typology2 : with ward criterion in 8 clusters for large data : a two step kmeans-cluster
km<-kmeans(mvad.emlt$coord,50)
hc<-hclust(dist(km$centers, method="euclidean"), method="ward")
zz<-cutree(hc, k=8)

##Plotting typology2 by clusters

seqdplot(mvad.seq, group=zz[km$cluster])



## Plotting the evolve of correlation betwwen states  
plot(mvad.emlt, from="employment", to="joblessness",type="cor")
plot(mvad.emlt, from=c("employment","HE", "school", "FE"), to="joblessness", delay=0, leg=TRUE)
plot(mvad.emlt, from="joblessness", to="employment", delay=6)
plot(mvad.emlt, type="pca", cex=0.4, compx=1, compy=2)




}
