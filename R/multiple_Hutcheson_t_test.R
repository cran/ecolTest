#' @title Multiple Hutcheson t-tests comparisons between communities.

#' @description This function computes the p-values of the Hutcheson t-test to
#' test the significance of the difference between more than two communities
#' Shannon diversity indexes, in a pairwise way.

#' @details This function performs Hutcheson's t-tests for comparing multiple
#' diversity indexes pairwise. This test is based on the Shannon diversity
#' index computed using a logarithm base specified by the user. The alternative
#'  hypothesis is one-sided, chosen automatically according to the sign of
#'  the difference between each pair of communities tested. The resulting
#'  p-values of the test are returned in a matrix. To see full details of
#'  the results of the test comparing two communities it is better
#'  to use \code{Hutcheson_t_test()}.

#' @note Missing values will be replaced with zero.

#' @param x Numeric dataframe or matrix of abundance of species per community 
#' sample. Columns must correspond to the samples and rows to species.

#' @param shannon.base Numeric value indicating the logarithm base for computing
#' the Shannon indexes. Defaults to \emph{exp(1)}.

#' @return A matrix whose entries are the p-values from the test result rounded to
#' five digits. Self-comparison elements (matrix diagonal) are flagged with \emph{NA}.
#' The names of the rows and columns are the names of the communities with
#' their Shannon diversity index.

#' @seealso See \code{\link{Hutcheson_t_test}} for Hutcheson's t-test details.

#' @author David Ramirez Delgado \email{linfocitoth1@gmail.com}.

#' @author Hugo Salinas \email{hugosal@comunidad.unam.mx}.

#' @examples

#' data("polychaeta_abundance")
#' multiple_Hutcheson_t_test(x = polychaeta_abundance,
#'                            shannon.base = 10)

#' @export

multiple_Hutcheson_t_test <- function(x, shannon.base=exp(1)) {
  x <- drop(as.matrix(x))
  if (!is.numeric(x)) {
    stop("x must be numeric")
    }
  if (any(x < 0, na.rm = TRUE)) {
    stop("x must be non-negative")
    }
  if (any(dim(x) < 2)) {
    stop("x must have at least two columns and rows")
    }
  nx <- ny <- 1:ncol(x)
  N <- colSums(x, na.rm = TRUE)
  H <- (N*log(N, shannon.base)-apply(x*log(x,shannon.base), 2,
                                     sum, na.rm = TRUE))/N
  names(ny) <- names(nx) <- paste(colnames(x), "H =", round(H, 2), sep = " ")
  vectorized_Hutcheson <- Vectorize(function(X, Y) {
    Hutcheson_t_test(x[, X], x[, Y],
                shannon.base = shannon.base, alternative = "auto")$p.value})
  p_values <- outer(nx, ny, FUN = vectorized_Hutcheson)
  p_values<- round(p_values, 5)
  for (i in 1:dim(p_values)[1]) {
    p_values[i, i] <- NA
  }
  return(list(p.values = p_values))
  }
