## Transition Matrixes

#' Calculates e^(-Lambda*t), where Lambda is a matrix.
#'
#' @param Lambda a square matrix of rates
#' @param t a time parameter (default 1).
#' @returns A matrix of the same shape as Lambda
eLt <- function (Lambda, t=1.0) {
  r <- eigen(Lambda)
  V <- r$vectors
  lam <- r$values
  V %*% diag(exp(-lam*t)) %*% solve(V)
}

#' Creates a transition probability matrix from the rate and time.
#'
#' @param Lambda a square matrix of transition rates -- the diagnoal is ingored.
#' @param t A time.
#' @returns a square matrix of transition probabilities for a time interval of size t.
poissonMat <- function(Lambda, t=1.0) {
  mat <- t*Lambda %*% eLt(Lambda,t)
  diag(mat) <- 0
  diag(mat) <- 1 - rowSums(mat)
  mat
}

#' Creates a transition rate matrix from a big a list of learning and forgetting rates.
#'
#' The `learn` parameter is the subdiagonal in the high-first orientation, and the `forget`
#' parameter is superdiagonal.  Thus, their lengths should be the rank of the matrix minus one.
#' As a special case, if the length is 1, then it will be replicated to the rank-1.
#'
#' States should be list of labels for the states.  If not supplied, it will be generated from the
#' based on the length of `learn` (or `forget`).
#'
#' @param learn A vector of learning rates.  Should be one less than number of states.
#' @param forget A vector of forgetting rates.  Should be one less than the number of states.
#' @param states A list of names for the states.
#' @param highFirst A logical value (default true) indicating whether the first state is the highest or lowest value in the sequence.
#' @returns A square matrix with the transition rates.
learnForgetRate <- function (learn, forget=0, states, highFirst=TRUE) {
  if (missing(states)) {
    K <- max(length(learn),length(forget))
    if (highFirst)
      states <- paste("Stage",K:0,sep="")
    else
      states <- paste("Stage",0:K,sep="")
  } else {
    if (!is.character(states))
      stop("Expected character vector for states argument.")
  }
  K <- length(states)
  if (!is.numeric(learn) && all(learn >= 0)) {
    stop ("Expected learn to be positive numbers, got ",learn)
  }
  if (!is.numeric(forget) && all(forget >=0)) {
    stop ("Expected forget to be positive numbers, got ",forget)
  }
  if (length(learn)==1L) learn <- rep(learn,K-1)
  if (length(forget)==1L) forget <- rep(forget,K-1)
  if (length(learn) != K-1L || length(forget) != K-1L) {
    stop("Length of learn and forget must be length(states) -1, or 1")
  }
  if (highFirst) {
    super <- forget; sub <- learn
  } else { ## Do I need this rev?
    super <- learn; sub <- forget
  }
  ratem <- matrix(Inf,K,K,dimnames=list(states,states))
  for (i in 1L:K) {
    if (i > 1L) {
      for (j in 1:(i-1L)) {
        ratem[i,j] <- prod(sub[j:(i-1L)])
      }
    }
    if (i < K-1L) {
      for (j in (i+1L):K) {
        ratem[i,j] <- prod(sub[i:(j-1L)])
      }
    }
  }
  ratem
}

