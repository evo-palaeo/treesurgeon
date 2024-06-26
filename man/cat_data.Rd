% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/treesurgeon.r
\name{cat_data}
\alias{cat_data}
\title{Combine phylogenetic data}
\usage{
cat_data(..., use.part.info = F, part.names = NULL)
}
\arguments{
\item{...}{two or more objects of class 'list' comprising named vectors of phylogenetic characters.}

\item{use.part.info}{logical. If TRUE, the function resturns a seperate list for each input list.}

\item{part.names}{a character vector of partition names e.g. 'dna', 'standard', 'protein'. If not supplied, the function will attempt to guess the correct partition type.}
}
\value{
A list, or multiple lists, of phylogenetic data of length T, where T equals the sum of unique taxa in the input lists.
}
\description{
Function to combine lists of phylogenetic data (e.g. those output by the ape function read.nexus.data()).
}
\details{
This function can be used to combine seperate lists of phylogenetic characters with non-overlapping taxa. If a taxon does not occur in all input lists, it will be coded as missing '?' for those input list(s). Can be used, for example, to combine molecular and morphological datasets with non-overlapping taxa. Note that taxa must be named identically in each input for the data to be combined succesfully.
}
\examples{
data(Lavoue2016)
combined_data <- cat_data(Lavoue2016$standard, Lavoue2016$dna, use.part.info = F) 
## Visualise 
if(!require("BiocManager", quietly = TRUE)){
	install.packages("BiocManager")
}
BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
df <- t(as.data.frame(combined_data))
df <- tolower(df)
df[df == "?"] <- NA
df[df == "-"] <- NA
df[df == "n"] <- NA
states <- as.character(na.omit(unique(as.character(df))))
cols <- rep("x", length(states))
names(cols) <- states
morph_states <- suppressWarnings(which(is.na(as.numeric(names(cols))) == F))
mol_states <- which(names(cols) \%in\% c("a", "c", "g", "t"))
cols[morph_states] <- hcl.colors(n = length(morph_states), palette = "Hawaii")
cols[mol_states] <- hcl.colors(n = length(mol_states), palette = "zissou")
cols[-c(mol_states, morph_states)] <- hcl.colors(n = length(cols[-c(mol_states, morph_states)]), palette = "Cividis")
Heatmap(df[,1:1500], row_names_side = "left", col = cols, na_col = "white", name = "states")

}
