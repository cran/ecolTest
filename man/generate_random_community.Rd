% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generate_random_community.R
\name{generate_random_community}
\alias{generate_random_community}
\title{Pseudo-random community generator}
\usage{
generate_random_community(
  H_index,
  shannon.base = exp(1),
  sp_num,
  ntotal,
  tol = 1e-04,
  maxiter = 100,
  silent = FALSE
)
}
\arguments{
\item{H_index}{The value of the generated communty Shannon diversity
index H'.}

\item{shannon.base}{A numeric value indicating the logarithm base for the
Shanon indices. Defaults to \emph{e}.}

\item{sp_num}{Numeric. The number of species in the generated community
sample.}

\item{ntotal}{Numeric. The total number of individuals in the generated
community sample.}

\item{tol}{Numeric. The tolerance of the diffrence of the diversity index
of the generated community relative abundance and the H_index argument.
Defaults to 0.0001.}

\item{maxiter}{Numeric. The maximum number of iterations to be performed.
Defaults to 100.}

\item{silent}{Logical. Indicates if convergence success or failure
messages are omitted. Defaults to FALSE.}
}
\value{
A list containing the following components:
\itemize{
\item community: If convergence is successful, a numeric vector with
the abundance of each species, else \emph{NA}
\item n_ite: If convergence is successful, the number of iterations
used to achieve convergence, else \emph{NA}
}
}
\description{
This function generates a random community sample dataset
that has a user-specified Shannon diversity index, species number, and
total abundance.
}
\details{
This function generates random community sample data: a numeric
vector whose elements represent a sample. The generated sample data is
composed of a number of species with a total number of individuals set by
the user, such that the community has a Shannon diversity index
approximately equal to a user-specified argument.
}
\examples{
# Generate a community with diversity index 2, composed of 20 species
# and 200 total individuals sampled
set.seed(26)
result <- generate_random_community(H_index = 2.7, sp_num = 20, ntotal = 200,
 maxiter = 300)
random_community <- result$community

#Compute H index
total <- sum(random_community)
-sum(random_community/total*log(random_community/total))

#Default maxiter argument will not converge
set.seed(26)
generate_random_community(H_index = 2.7, sp_num = 20, ntotal = 200)
}
\author{
David Ramirez Delgado \email{linfocitoth1@gmail.com}

Hugo Salinas \email{hugosal@comunidad.unam.mx}
}
