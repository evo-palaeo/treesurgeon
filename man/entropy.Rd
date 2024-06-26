% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/treesurgeon.r
\name{entropy}
\alias{entropy}
\title{Shannon Entropy}
\usage{
entropy(x)
}
\arguments{
\item{x}{a numeric vector of probabilities that sum to 1.}
}
\value{
The Shannon Entropy score.
}
\description{
Function to calculate the Shannon Entropy of a vector or probabilities.
}
\details{
to do!
}
\examples{
## Load data
head(vert_data$morph)
tp <- get_tip_priors(vert_data$morph[,2], extra_state = T)
colnames(tp[[1]]) <- c("bone absent", "secondary osteons absent", "secondary osteons present")
fitER <- fitMk(tree = vert_data, x = tp[[1]], model = "ER")
ancER <- ancr(fitER)
## Calculate mean Entropy
mean(apply(ancER$ace, 1, entropy))
}
