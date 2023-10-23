library(readxl)
library(zoo)

data <- read_excel("/Users/carlthuresson/Documents/Table 1.xls")
yy <- log(data$rgdpq) #GDP in logs
y <- 100 * diff(yy) #GDP growth in logs
T <- length(y) #Number of quarters

#A number of globals that they have in matlab, not sure what they do
mstarU <- 5 #Nr of lags U-shaped recession persists
mstarL <- 5 #Nr of lags L-shaped recession persists
growth_bdate <- 236 #Assumption of where structural break in growth exist = 2006Q1
vol_bdate <- 149 #Assumption of structural break in volatility exist = 1984Q2

prmtr_in <- c(-3.5, -3.9, 0.6, 0.9, 0.9, -1.3, -2.0, -0.4, -0.2, -1.7) #Initial parameter values for function

#LIKELIHOOD FUNCTION

lik_fcn <- function(prmtr_in, T, y, mstarU, mstarL, vol_bdate, growth_bdate){
  no_st <- 8 #Nr of periods to consider
  dmnsion <- 3^no_st #Nr of cases to consider
  
  #Parameters
  p01 <- prmtr_in[1] #Prob. of going from normal to L-shaped recession
  p02 <- prmtr_in[2] #Prob. of going from normal to U-shaped recession
  p11 <- prmtr_in[3] #Prob. of staying in L-shaped recession
  p22 <- prmtr_in[4] #Prob. of staying in U-shaped recession
  mu0 <- prmtr_in[5] #Growth in normal
  mu1 <- prmtr_in[6] #Effect on growth in L-shaped recession
  mu2 <- prmtr_in[7] #Effect on growth in U-shaped recession
  omega1 <- 0 #No bounceback effect for growth during L-shaped recession
  omega2 <- -mu2/mstarU #Bounceback effect on growth for U-shaped recession
  delta <- prmtr_in[8] #Change in trend growth
  sigma2_pre <- prmtr_in[9]^2 #Volatility of growth before structural break
  sigma2_post <- prmtr_in[10]^2 #Volatility of growth after structural break
  
  PI <- matrix(0, nrow =3, ncol =3) #Matrix of probability for going from one state to another
  PI[1, 1] <- 1 - p01 - p02 #Staying in N
  PI[1, 2] <- p01 #N to L
  PI[1, 3] <- p02 #N to U
  PI[2, 1] <- p11 #L to N
  PI[2, 2] <- 1 - p11 #Staying in L
  PI[2, 3] <- 0 #L to U
  PI[3, 1] <- 1 - p22 #U to N
  PI[3, 2] <- 0 #U to L
  PI[3, 3] <- p22 #Staying in U
  
  A <- rbind(diag(3) - PI, rep(1, 3)) #Complementary probabilities for staying in same state
  en <- c(rep(0, 3), 1) #Vector of zeros and one at the end
  prob__t <- solve(t(A) %*% A) %*% (t(A) %*% en) #Calculates the steady-state prob.(=prob. of staying in one state)
  
  pr_trf0 <- as.vector(PI) #Reshapes PI to a vector
  pr_trf1 <- c(pr_trf0, pr_trf0, pr_trf0) #Stacks the vector on top of itself
  pr_trf2 <- c(pr_trf1, pr_trf1, pr_trf1) #Again, same thing
  pr_trf3 <- c(pr_trf2, pr_trf2, pr_trf2)
  pr_trf4 <- c(pr_trf3, pr_trf3, pr_trf3)
  pr_trf5 <- c(pr_trf4, pr_trf4, pr_trf4)
  pr_trf <- c(pr_trf5, pr_trf5, pr_trf5) #Resulting in a long ass vector of steady state probabilities
  
  #Stack steady state probabilities on top of each other
  prob__t0 <- c(prob__t, prob__t, prob__t) 
  prob__t1 <- c(prob__t0, prob__t0, prob__t0)
  prob__t2 <- c(prob__t1, prob__t1, prob__t1)
  prob__t3 <- c(prob__t2, prob__t2, prob__t2)
  prob__t4 <- c(prob__t3, prob__t3, prob__t3)
  prob__t5 <- c(prob__t4, prob__t4, prob__t4)
  prob__t6 <- c(prob__t5, prob__t5, prob__t5)
  
  #Converting the vectors to a column and multiplying with probability vector to get unconditional probabilities
  prob__t00 <- matrix(prob__t0, nrow = length(prob__t0), ncol = 1) * pr_trf0 
  prob__t11 <- matrix(prob__t1, nrow = length(prob__t1), ncol = 1) * pr_trf1
  prob__t22 <- matrix(prob__t2, nrow = length(prob__t2), ncol = 1) * pr_trf2
  prob__t33 <- matrix(prob__t3, nrow = length(prob__t3), ncol = 1) * pr_trf3
  prob__t44 <- matrix(prob__t4, nrow = length(prob__t4), ncol = 1) * pr_trf4
  prob__t55 <- matrix(prob__t5, nrow = length(prob__t5), ncol = 1) * pr_trf5
  prob__ <- matrix(prob__t6, nrow = length(prob__t6), ncol = 1) * pr_trf #Unconditional prob. for 3*8 cases
  
  #Creating selection matrices for the Markov-switching
  sU <- c(0, 0, 1)
  
  sU0 <- kronecker(ones_vector, kronecker(sU, matrix(1, nrow = 1, ncol = 1)))
  sU1 <- kronecker(rep(1, 3^6), kronecker(sU,rep(1, 3^1)))
  sU2 <- kronecker(rep(1, 3^5), kronecker(sU,rep(1, 3^2)))
  sU3 <- kronecker(rep(1, 3^4), kronecker(sU,rep(1, 3^3)))
  sU4 <- kronecker(rep(1, 3^3), kronecker(sU,rep(1, 3^4)))
  sU5 <- kronecker(rep(1, 3^2), kronecker(sU,rep(1, 3^5)))
  sU6 <- kronecker(rep(1, 3^1), kronecker(sU,rep(1, 3^6)))
  sU7 <- kronecker(rep(1, 3^0), kronecker(sU,rep(1, 3^7)))
  
  sUmat <- matrix(c(sU7, sU6, sU5, sU4, sU3, sU2, sU1, sU0), ncol = 8, byrow = FALSE)
  
  sL <- c(0, 1, 0)
  
  sL0 <- kronecker(ones_vector, kronecker(sL, matrix(1, nrow = 1, ncol = 1)))
  sL1 <- kronecker(rep(1, 3^6), kronecker(sL,rep(1, 3^1)))
  sL2 <- kronecker(rep(1, 3^5), kronecker(sL,rep(1, 3^2)))
  sL3 <- kronecker(rep(1, 3^4), kronecker(sL,rep(1, 3^3)))
  sL4 <- kronecker(rep(1, 3^3), kronecker(sL,rep(1, 3^4)))
  sL5 <- kronecker(rep(1, 3^2), kronecker(sL,rep(1, 3^5)))
  sL6 <- kronecker(rep(1, 3^1), kronecker(sL,rep(1, 3^6)))
  sL7 <- kronecker(rep(1, 3^0), kronecker(sL,rep(1, 3^7)))
  
  sLmat <- matrix(c(sL7, sL6, sL5, sL4, sL3, sL2, sL1, sL0), ncol = 8, byrow = FALSE)
  
  val <- 0 #Stored value that will be updated
  d <- 0 #Same here
  dd <- 0 #And same here
  j_iter <- 1 #Number of iterations
  
  while (j_iter <= T) { #As long as number of iterations smaller than T
    if (j_iter > growth_bdate) { 
      d <- 1 #If number of iterations are bigger than the ass. of struc. break in growth, set d to 1
    }
    if (j_iter > vol_bdate) {
      dd <- 1 #If number of iterations are bigger than the ass. of struc. break in volatility, set dd to 1
    }
    
    