% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/treesurgeon.r
\name{multRaw}
\alias{multRaw}
\title{Multiclass Raw error}
\usage{
multRaw(prediction, truth)
}
\arguments{
\item{prediction}{a matrix of predicted probabilities.}

\item{truth}{a numeric vector of length equal to the number of classes with each element encoding the true outcome. Alternatively, a one-hot encoded matrix with rows equal to the number of classes and columns equal to the number of outcomes.}
}
\value{
The multiclass Raw error.
}
\description{
Function to calculate the Multiclass Raw error given a matrix of predicted probabilities and a vector or matrix of actual outcomes.
}
\details{
to do!
}
\examples{
## Simulate some tip data
t1 <- rtree(20)
Q <- structure(c(-0.5, 0.4, 0.05, 0.3, -0.5, 0.5, 0.2, 0.1, -0.55), dim = c(3L, 3L), dimnames = list(c("1", "2", "3"), c("1", "2", "3")))
sim_x <- sim.Mk(t1, Q, internal = T)
tip_states <- head(sim_x, Ntip(t1))
node_states <- tail(sim_x, Nnode(t1))
## Estimate ancestral states
fitARD <- fitMk(t1, x = tip_states, model = "ARD")
ancARD <- ancr(fitARD)
 ## Calculate Multiclass Raw error
multRaw(prediction = ancARD$ace, truth = node_states)
}
