#' @title Pseudo-random community generator

#' @description This function generates a random community sample dataset 
#' that has a user-specified Shannon diversity index, species number, and 
#' total abundance.


#' @details This function generates random community sample data: a numeric 
#' vector whose elements represent a sample. The generated sample data is 
#' composed of a number of species with a total number of individuals set by 
#' the user, such that the community has a Shannon diversity index 
#' approximately equal to a user-specified argument.

#' @param H_index The value of the generated communty Shannon diversity 
#' index H'.

#' @param shannon.base A numeric value indicating the logarithm base for the
#' Shanon indices. Defaults to \emph{e}.

#' @param sp_num Numeric. The number of species in the generated community 
#' sample.

#' @param ntotal Numeric. The total number of individuals in the generated 
#' community sample.

#' @param tol Numeric. The tolerance of the diffrence of the diversity index 
#' of the generated community relative abundance and the H_index argument. 
#' Defaults to 0.0001.

#' @param maxiter Numeric. The maximum number of iterations to be performed.
#' Defaults to 100.

#' @param silent Logical. Indicates if convergence success or failure
#'  messages are omitted. Defaults to FALSE.

#' @return A list containing the following components:

#' \itemize{
#'   \item community: If convergence is successful, a numeric vector with
#'   the abundance of each species, else \emph{NA}
#'   \item n_ite: If convergence is successful, the number of iterations
#'    used to achieve convergence, else \emph{NA}
#' }

#' @author David Ramirez Delgado \email{linfocitoth1@gmail.com}

#' @author Hugo Salinas \email{hugosal@comunidad.unam.mx}

#' @examples

#' # Generate a community with diversity index 2, composed of 20 species
#' # and 200 total individuals sampled
#' set.seed(26)
#' result <- generate_random_community(H_index = 2.7, sp_num = 20, ntotal = 200,
#'  maxiter = 300)
#' random_community <- result$community
#'
#' #Compute H index
#' total <- sum(random_community)
#' -sum(random_community/total*log(random_community/total))
#'
#' #Default maxiter argument will not converge
#' set.seed(26)
#' generate_random_community(H_index = 2.7, sp_num = 20, ntotal = 200)

#' @importFrom stats runif

#' @export


generate_random_community <- function(H_index, shannon.base=exp(1), sp_num,
                                      ntotal, tol=0.0001,
                                      maxiter=100, silent=FALSE) {
  if (sp_num < 2) {
    stop("Species number must be > 1")
    }
  if (H_index > log(sp_num,base = shannon.base)) {
    stop("Impossible community (H_index>log(N))")
    }
  if (H_index < 0) {
    stop("H_index must be positive")
    }
  if (!requireNamespace("stats", quietly = TRUE)) {
    stop('Package "stats" must be installed')
    }
  H_index <- -H_index
  a <- 1/sp_num
  community <- runif(sp_num,0,a)
  community <- community+((1-sum(community))*a)
  cost <- function(comun) {
          return(abs(H_index-sum(comun*log(comun, shannon.base))))
    }
  choose_donor <- function(com, ntotal, step) {
    for (d in sample(1:length(com))) {
      if ((com[d]-(com[d]*step))*ntotal > 1) {
        return(d)
        }
      }
    d}
  iter <- 0
  while (abs(H_index - sum(community*log(community, shannon.base),
                           na.rm = T)) > tol) {
    step_size <- 0.1
      if (iter > maxiter) {
        if (!(silent)) {
          message("Convergence failed")}
        return(list(community=NA, n_iter=NA))
      }
    donor <- choose_donor(com = community, ntotal = ntotal, step = step_size)
    for (recipient in sample((1:sp_num)[-donor])) {
      neighbor <- community
      neighbor[recipient] <- neighbor[recipient]+(neighbor[donor]*step_size)
      neighbor[donor] <- neighbor[donor]-(neighbor[donor]*step_size)
      if (cost(neighbor) < cost(community)) {
        community <- neighbor
        break
      }
    }
    iter <- iter+1
  }
  if (!(silent)) {
  message(paste("Convergence succesfull, No. iteration = ", iter))}
  return(list(community=round(community*ntotal), n_iter=iter))
  }
