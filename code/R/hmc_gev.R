# This is slighty modified.  I added the "others" input and
#  changed the format of the output
#
# SIMPLE IMPLEMENTATION OF HAMILTONIAN MONTE CARLO.
#
# Radford M. Neal, 2010.
#
# This program appears in Figure 2 of "MCMC using Hamiltonian dynamics",
# to appear in the Handbook of Markov Chain Monte Carlo.
#
# The arguments to the HMC function are as follows:
#
#   U          A function to evaluate minus the log of the density of the
#              distribution to be sampled, plus any constant - ie, the
#              "potential energy".
#
#   grad_U     A function to evaluate the gradient of U.
#
#   epsilon    The stepsize to use for the leapfrog steps.
#
#   L          The number of leapfrog steps to do to propose a new state.
#
#   current_q  The current state (position variables only).
#
# Momentum variables are sampled from independent standard normal
# distributions within this function.  The value return is the vector
# of new position variables (equal to current_q if the endpoint of the
# trajectory was rejected).
#
# This function was written for illustrative purposes.  More elaborate
# implementations of this basic HMC method and various variants of HMC
# are available from my web page, http://www.cs.utoronto.ca/~radford/

# ------------------------------------------------------------------
# This version of the function returns (as attributes of the result)
# information for assessing how well the No U-Turn sampler might work.
# The the distance of each position along the trajectory from the
# start position is returned as "dist", and with whether the
# trajectory end-point was accepted is returned as "accept".
# ------------------------------------------------------------------


gevHMC = function (U, grad_U, current_q, epsilon = 0.01, L = 10,
                data, beta, xi, a, b, alpha, rho, calc, others, this.param)
{
  dist = numeric(L + 1)

  if (this.param == "a_alpha") {
    infinite <- c(FALSE, FALSE)
  } else {
    infinite <- FALSE
  }

  q = current_q
  p = rnorm(length(q), 0, 1)  # independent standard normal variates
  current_p = p

  # Make a half step for momentum at the beginning

  p = p - epsilon * grad_U(q = q, data = data, beta = beta, xi = xi, a = a,
                           b = b, alpha = alpha, rho = rho, calc = calc,
                           others = others) / 2
  # Alternate full steps for position and momentum

  for (i in 1:L)
  {
    # Make a full step for the position

    q = q + epsilon * p
    if (any(q == Inf | q == -Inf)) {

      if (this.param == "a_alpha") {
        infinite = c(FALSE, FALSE)
        if(any(q[1:(length(q) - 1)] == Inf | q[1:(length(q) - 1)] == -Inf)) {
          print("Proposal is Inf or -Inf for log(a), automatically reject")
          infinite[1] <- TRUE
        }
        if (tail(q, 1) == Inf | tail(q, 1) == -Inf) {
          int("Proposal is Inf or -Inf for logit(alpha), automatically reject")
          infinite[2] <- TRUE
        }
        return(list(q = current_q, accept = FALSE, infinite = infinite))
      } else {
        print(paste("Proposal variable is Inf for ", this.param,
                    ", automatically reject", sep = ""))
        return(list(q = current_q, accept = FALSE, infinite = TRUE))
      }
    }
    if (this.param == "a_alpha") {
      if (any(exp(q[1:(length(q) - 1)]) == 0)) {
        print("Proposal is 0 for a")
        infinite[1] <- TRUE
      }
      if (any(exp(q[1:(length(q) - 1)]) == Inf)) {
        print("Proposal is Inf for a")
        infinite[1] <- TRUE
      }
      if (transform$inv.logit(tail(q, 1)) == 0) {
        print("Proposal is 0 for alpha")
        infinite[2] <- TRUE
      }
      if (transform$inv.logit(tail(q, 1)) == Inf) {
        print("Proposal is Inf for alpha")
        infinite[2] <- TRUE
      }
      if (any(infinite)) {
        return(list(q = current_q, accept = FALSE, infinite = infinite))
      }
    }

    dist[i + 1] = sqrt(sum((q - current_q)^2))

    # Make a full step for the momentum, except at end of trajectory

    if (i != L) p = p - epsilon * grad_U(q = q, data = data, beta = beta,
                                         xi = xi, a = a, b = b, alpha = alpha,
                                         rho = rho, calc = calc,
                                         others = others)
    if (any(is.nan(p))) {
      if (this.param == "a_alpha") {
        infinite <- c(FALSE, FALSE)
        if (any(is.nan(p[1:(length(p) - 1)]))) {
          print("Momentum variable is NaN for a, automatically reject")
          infinite[1] <- TRUE
        }
        if (tail(is.nan(p), 1)) {
          print("Momentum variable is NaN for alpha, automatically reject")
          infinite[2] <- TRUE
        }
        return(list(q = current_q, accept = FALSE, infinite = infinite))
      } else {
        print(paste("Momentum variable is NaN for ", this.param,
                    ", automatically reject", sep = ""))
        return(list(q = current_q, accept = FALSE, infinite = TRUE))
      }
    } else if (any(p == Inf | p == -Inf)) {
      if (this.param == "a_alpha") {
        infinite <- c(FALSE, FALSE)
        if(any(p[1:(length(p) - 1)] == Inf | p[1:(length(p) - 1)] == -Inf)) {
          print("Momemtum variable is Inf for a, automatically reject")
          infinite[1] <- TRUE
        }
        if (tail(q, 1) == Inf | tail(q, 1) == -Inf) {
          print("Momemtum variable is Inf for alpha, automatically reject")
          infinite[2] <- TRUE
        }
        return(list(q = current_q, accept = FALSE, infinite = infinite))
      } else {
        print(paste("Momentum variable is Inf for ", this.param,
                    ", automatically reject", sep = ""))
        return(list(q = current_q, accept = FALSE, infinite = TRUE))
      }
    }
  }

  # Make a half step for momentum at the end.

  p = p - epsilon * grad_U(q = q, data = data, beta = beta, xi = xi, a = a,
                           b = b, alpha = alpha, rho = rho, calc = calc,
                           other = others) / 2

  # Negate momentum at end of trajectory to make the proposal symmetric

  p = -p

  # Evaluate potential & kinetic energies at start & end of trajectory

  current_U = U(current_q, data = data, beta = beta, xi = xi, a = a, b = b,
                alpha = alpha, rho = rho, calc = calc, others = others)
  current_K = sum(current_p^2) / 2
  proposed_U = U(q, data = data, beta = beta, xi = xi, a = a, b = b,
                 alpha = alpha, rho = rho, calc = calc, others = others)
  proposed_K = sum(p^2) / 2

  R <- current_U - proposed_U + current_K - proposed_K

  # Accept or reject the state at end of trajectory, returning either
  # the position at the end of the trajectory or the initial position
  if (!is.nan(R)) {
    if (runif(1) < exp(R))
    {
      return(list(q = q, accept = TRUE, infinite = infinite))
    }
    else
    {
      return(list(q = current_q, accept = FALSE, infinite = infinite))
    }
  } else {
    if (this.param == "a_alpha") {
      infinite = c(FALSE, FALSE)
      if (tail(is.nan(p), 1)) {
        infinite[2] <- TRUE
      }
      if (any(is.nan(p[1:(length(p) - 1)]))) {
        infinite[1] <- TRUE
      }
      return(list(q = current_q, accept = FALSE, infinite = infinite))
    } else {
      return(list(q = current_q, accept = FALSE, infinite = infinite))
    }
  }
}
