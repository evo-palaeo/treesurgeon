% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/treesurgeon.r
\name{csv_to_phyDat}
\alias{csv_to_phyDat}
\title{Comma-seperated values to phyDat}
\usage{
csv_to_phyDat(file, ...)
}
\arguments{
\item{file}{the name of the file which the data are to be read from.}

\item{...}{Further arguments to be passed to read.csv().}
}
\value{
An object of class phyDat.
}
\description{
Function to convert a taxon-by-character matrix of categorical data in csv format into a phyDat object.
}
\details{
This function is a conveniant shortcut for importing a taxon-by-character matrix of categorical data saved as a csv file and converting it directly into phyDat objects. Categorical tokens (e.g. 0, 1, 2... -, ?) are interpreted using the treesurgeon function get_contrast().
}
\examples{
## Load data
data(KeatingDonoghue)
## write as csv file
write.csv(KeatingDonoghue, file = "temp.csv")
## Import using the function
phy_dat <- csv_to_phyDat(file = "temp.csv", row.names = 1)
phy_dat
## delete temporary csv
file.remove("temp.csv")

}
