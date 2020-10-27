# Side Effect Pipeline Part 1
# input setup creating input fields and user-created functions

library(gmp)
library(stringr)
# create inputs
# inputs (api_key, drug_name, email)
api_key <- [Enter your personal openFDA Key here]
drug_name <- [Enter your selected drug]
drug_name <- str_to_upper(drug_name)
user_agent <- [Enter an email address for API use]

# Framework of Other Drugs
checkQuery <- function(x) {
  if(http_error(x)){
    warning("The request failed")
  } else {
    content(x, as = "text")
  }
}

######################## Create the EBGM Function #######################################3
EBGM <- function(A, B, C, D) {
  N <- A
  E <- (C*B)/D
  RR <- N/E
  PRR <- (A/B)/(C/D)
  
  n <- as.bigz(N)
  e <- as.bigz(E)
  
  logP <- E/log(10) - N*log10(E) + log10(factorial(n)) - log10(1 + E/(N + 1) + E^2/((N + 1)*(N + 2)) + E^3/((N + 1)*(N + 2)*(N + 3)))
  
  a1 <- .2
  B1 = .1
  a2 <- as.bigz(2)
  B2 = 4
  P = 1/3
  
  
  f_nabe_1 <- ((1+B1/E)^-N)*((1+E/B1)^-a1)*((gamma(a1+n))/(gamma(a1)*factorial(n)))
  f_nabe_2 <- ((1+B2/E)^-N)*((1+E/B2)^-a2)*((gamma(a2+n))/(gamma(a2)*factorial(n)))
  
  Qn <- P*f_nabe_1/(P*f_nabe_1 + (1-P)*(f_nabe_2))
  
  a2 <- 2
  
  E_log_gamma_given_N <- Qn*(digamma(a1+N)-log(B1+E))+(1-Qn)*(digamma(a2+N)-log(B2+E))
  EBlog <- as.numeric(E_log_gamma_given_N/log(2))
  EBGM <- 2^EBlog
  
  EB05 <- EBGM*exp(-2/sqrt(N+1))
  EB95 <- EBGM*exp(2/sqrt(N+1))
  return(data.frame(N,
                    E,
                    RR,
                    PRR,
                    EBGM,
                    EB05,
                    EB95))
}
