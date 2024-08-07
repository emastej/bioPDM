% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getDirectionsOfMed.R
\name{getDirectionsOfMed}
\alias{getDirectionsOfMed}
\title{getDirectionsOfMed}
\usage{
getDirectionsOfMed(
  data_list = NULL,
  nPDM = 5,
  doJointPDM = FALSE,
  doSparseThresh = FALSE,
  doSparseEN = FALSE,
  doBootPDM = FALSE,
  doBootJPDM = FALSE,
  bootSamp = 1000,
  returnBootSamples = FALSE,
  saveResults = FALSE,
  saveDir = NULL,
  notes = NULL,
  numCores = 1,
  timeout = 5
)
}
\arguments{
\item{data_list}{List containing X, Y, M, M_reduced, tLoadingMatrix,
and exVar elements}

\item{nPDM}{Number of PDMs to calculate}

\item{doJointPDM}{Calculate the JointPDM. Can be True or False}

\item{doSparseThresh}{Calculate sparse PDM using threshold method. True or False}

\item{doSparseEN}{Calculate sparse PDM using elastic net method. True or False}

\item{doBootPDM}{Bootstrap samples to calculate feature weight stats
(95\% CI, mean, standard deviation). Can be True or False}

\item{doBootJPDM}{Bootstrap samples to calculate feature weight stats
(95\% CI, mean, standard deviation) for JointPDM. Can be True or False}

\item{bootSamp}{Number of bootstrap samples}

\item{returnBootSamples}{Return the PDM weights and path coefficients
calculated from bootstrapped samples. Can be True or False}

\item{saveResults}{Save final list and a .txt file with parameters and any
free text notes}

\item{saveDir}{Directory for which the user wants to save results to}

\item{notes}{Parameter to write a free text note that will be saved in .txt
file (if saveResults = TRUE)}

\item{numCores}{Number of cores to use for parallel execution. Default is
a single core}

\item{timeout}{Time allowed for optimization for each set of starting values.
Reducing timeout may reduce time required for calculating PDMs, but reducing
it too much may lead to reduced accuracy of results.}
}
\value{
A list
\itemize{
\item redFeatWeights - matrix (b x k) Weights for each reduced b feature
calculated for each of the k PDMs
\item featWeights - matrix (p x k) Weights mapped back onto original Mediator
Matrix
\item jointWeights - vector (p x 1) Weights for the JointPDM. Each Mediator
feature has one weight (if applicable)
\item pathCoeff - matrix (5 x k) Mediation Path Coefficients for each of the
k PDMs (c, c', a, b, ab)
\item PDMk - vector (p x 1) PDM vector for each k PDM
\item JointPDM - vector (p x 1) JointPDM
\item Boot - List of bootstrap results (if applicable)
\itemize{
\item featWeightStats -  95\% CI, mean, and standard deviation
(for each feature for each PDM) calculated from the PDM weights
generated from bootstrapped samples
\item featWeightSamples -  The PDM weights generated from bootstrapped
samples (from which the featWeightStats were calculated)
\item pathCoeffSamples - Resulting mediation pathway coefficients from
each bootstrapped sample
\item jointWeightStats - 95\% CI, mean, and standard deviation
(for each feature) calculated from JointW weights generated from
bootstrapped results
\item jointWeightSamples - The JointPDM weights generated from bootstrapped
samples (from which the jointWeightStats were calculated)
}
}
}
\description{
Calculates Principal Directions of Mediation
(PDMs) that mediate the effect of an exposure on an outcome.
}
