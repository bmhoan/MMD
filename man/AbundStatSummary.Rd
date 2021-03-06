% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stat_report.R
\name{AbundStatSummary}
\alias{AbundStatSummary}
\title{Function for auto-generate statistical summary of any abudance table (OTU,GENE,KO...)}
\usage{
AbundStatSummary(abund.tab, sp.meta)
}
\arguments{
\item{abund.tab}{as.matrix/as.data.frame: otu,ko, rownames=OTU/GENE/KO names, colnames=sample name}

\item{sp.meta}{sample metadata as.matrix (ID,NICKNAME,CONDITION,SUBJECT,CATEGORY)}
}
\description{
Function for auto-generate statistical summary of any abudance table (OTU,GENE,KO...)
}
\examples{
AbundStatSummary(OTUDATA$otu.tab,OTUDATA$sp.meta)
}
