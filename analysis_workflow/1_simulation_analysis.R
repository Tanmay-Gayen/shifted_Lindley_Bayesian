##install.packages("PearsonDS")
##install.packages("nleqslv")

library(PearsonDS)
library(nleqslv)
set.seed(23)
n1=50#####change
mu=1.5######### Change
theta=1.5######### Change
count1=0;
count2=0;
a=matrix(0,1,n1);
b=matrix(0,1,n1);
w1=((1+mu)*theta)/(1+theta*(1+mu))

for(i in 1:n1)
{
  if(runif(1)<w1)
  {
    count1=count1+1
    a[count1]=rpearsonIII(1, shape=1, location=mu, scale=1/theta)
  } else
  {
    count2=count2+1
    b[count2]=rpearsonIII(1, shape=2, location=mu, scale=1/theta)  
  }
}
c1=c(a,b)
d1=c1[c1!=0]
d1
min_value <- min(d1)
print(min_value)
mu1=min_value
mu1
message <- sprintf("The M.L.E of mu1 is = %f", min_value)
print(message)

######################################################################
samp=d1
m31=min_value
x=n1*mu1*(1+mu1)-sum(d1)*(1+mu1)
x
y=n1+n1*mu1-sum(d1)+n1*mu1
y
# Define the coefficients
a <- x
b <- y
c <- 2 * n1

# Calculate the discriminant
discriminant <- b^2 - 4 * a * c
discriminant

# Check if there are real roots
if (discriminant < 0) {
  print("No real roots")
} else if (discriminant == 0) {
  # Calculate the single real root
  root <- -b / (2 * a)
  print(paste("Single real root:", root))
} else {
  # Calculate two real roots
  root1 <- (-b + sqrt(discriminant)) / (2 * a)
  root2 <- (-b - sqrt(discriminant)) / (2 * a)
  print(paste("Root 1:", root1))
  print(paste("Root 2:", root2))
  
}
# Select the positive root
if (root1 > 0) {
  positive_root <- root1
} else {
  positive_root <- root2
}
# Print the positive root
print(positive_root)
theta1=positive_root
theta1
##################################
t=1
########################################################################
###true survival function
R_t_MLE <- function(t, theta1, mu1) {
  result <- ((1 + theta1 * (1 + t)) / (1 + theta1 * (1 + mu1))) * exp(-theta1* (t - mu1))
}
result_value_11<- R_t_MLE(t,theta1,mu1)
print(result_value_11)
R_t_MLE=result_value_11
R_t_MLE

#################################

#####We take the value of theta1 is 1.133, since theta1 is always positive
L11=-(2*n1/theta1^2)+((n1*(1+mu1)^2)/(1 + theta1 * (1 + mu1))^2)
L11
L22=((n1 * theta1^2) / (1 + theta1 * (1 + mu1))^2)
L22
##The value of L12=L21
L21=-(n1/(1 + theta1 * (1 + mu1))^2)+n1
L21
L12=-(n1/(1 + theta1 * (1 + mu1))^2)+n1
L12

L111=(4*n1/theta1^3)-((2*n1*(1+mu1)^3)/(1+theta1*(1+mu1))^3)
L111
L222 <- -((2 * n1 * theta1^3) / (1 + theta1 * (1 + mu1))^3)
L222 

L112 <- ((2 * n1 * (1 + mu1) * (1 + theta1 * (1 + mu1))^2) - (2 * n1 * theta1 * (1 + mu1)^2 * (1 + theta1 * (1 + mu1)))) / ((1 + theta1 * (1 + mu1))^4)
L112

# Create the matrix
matrix<- matrix(c(-L11, -L12, -L21, -L22), nrow = 2, byrow = TRUE)

# Print the matrix
print(matrix)

# Find the inverse of the matrix
inverse_matrix <- solve(matrix)

# Print the inverse matrix
print(inverse_matrix)
T11=inverse_matrix[1,1]
T12=inverse_matrix[1,2]
T21=inverse_matrix[2,1]
T22=inverse_matrix[2,2]

#########change

###T11=0.0001797935
##T12=-0.001927182 
##T21=-0.001927182
##T22=-0.002337156
suppressWarnings({a1 <- a2 <- b1 <- b2 <- 2})
eta1=((a1-1)/theta1)-b1
eta2=((a2-1)/mu1)-b2


##The squared error loss function of theta is
theta_SE=theta1+T11*eta1+T21*eta2+0.5*(L111*T11^2+L222*T22*T12+3*L112*T12*T11)
theta_SE
mu_SE <- mu1 + T22 * eta2 + T12 * eta1 + 0.5 * (L111 * T11 * T12 + L222 * T22^2 + L112 * (T11 * T22 + 2 * T12^2))
mu_SE 

#################################################################
#########################################################################################################  
##Reliability Function (The ML and Bayes estimators of R(t)are evaluated at t=1)
R <- function(t, theta1, mu1) {
  result <- ((1 + theta1 * (1 + t)) / (1 + theta1 * (1 + mu1))) * exp(-theta1 * (t - mu1))
}
result_value1 <- R(t,theta1,mu1)
print(result_value1)
R=result_value1
R

#########################################################################################################  
########################################################################################################
R1 <- function(t, theta1, mu1) {
  numerator <- (1 + theta1 * (1 + t)) * (theta1^2 * (1 + mu1) * exp(mu1 * theta1))
  denominator <- exp(t * theta1) * (1 + theta1 * (1 + mu1))^2
  
  result2 <- numerator / denominator
  
  return(result2)
}              
result2 <- R1(t,theta1, mu1)
print(result2) 
R1=result2
R1
#############################################
R11 <- function(t, theta1, mu1) {
  numerator1 <- (1 + theta1 * (1 + t))
  denominator1 <- exp(t * theta1)
  
  numerator2 <- (1 + theta1 * (1 + mu1))^2 * (theta1^2 * exp(theta1 * mu1) + exp(theta1 * mu1) * theta1^3 * (1 + mu1)) - 
    2 * theta1^3 * (1 + mu1) * exp(theta1*mu1) * (1 + theta1 * (1 + mu1))
  denominator2 <- (1 + theta1* (1 + mu1))^4
  
  result3 <- (numerator1 / denominator1) * (numerator2 / denominator2)
  
  return(result3)
}
result3 <-R11(t, theta1,mu1)
print(result3)
R11=result3
R11

#######################################################                
R12<- function(t, theta1, mu1) {
  numerator1 <- (1 + theta1* (1 + t))
  denominator1 <- exp(t * theta1)
  numerator2 <- ((1 + theta1* (1 + mu1))^2) * (2 * theta1* exp(theta1* mu1) * (1 + mu1) + theta1^2 * exp(theta1* mu1) * mu1* (1 + mu1)) - 
    theta1^2 * (1 + mu1) * exp(theta1* mu1) * 2 * ((1 + theta1* (1 + mu1)) * (1 + mu1))
  denominator2 <- (1 + theta1* (1 + mu1))^4
  
  numerator3 <- (theta1^2 * exp(theta1* mu1) * (1 + mu1))
  denominator3 <- (1 + theta1* (1 + mu1))^2
  
  numerator4 <- ((1 + t) * exp(t * theta1) + (1 + theta1 * (1 + t)) * t * exp(t * theta1))
  denominator4 <- (exp(t * theta1))^2
  
  result4 <- ((numerator1 / denominator1) * (numerator2 / denominator2)) + ((numerator3 / denominator3) * (numerator4 / denominator4))
  
  return(result4)
}
result4 <- R12(t,theta1, mu1)
print(result4)
R12=result4
R12
#######################################################################
R2<- function(t, theta1, mu1) {
  numerator <- theta1 * (mu1^2 + 2 * mu1 - t^2 - 2 * t) + theta1^2 * (mu1 + mu1^2 + mu1^2 * t - t - t^2 - mu1 * t^2) * exp(-theta1 * (t - mu1))
  denominator <- (1 + theta1 * (1 + mu1))^2
  
  result5 <- numerator / denominator
  
  return(result5)
}
result5 <-R2(t,theta1,mu1)
print(result5)
R2=result5
R2
###################################################################################
R22<- function(t, theta1, mu1) {
  term1_numerator <- theta1 * (mu1^2 + 2 * mu1 - t^2 - 2 * t) + theta1^2 * (mu1 + mu1^2 + mu1^2 * t - t - t^2 - mu1 * t^2) * exp(-theta1 * (t - mu1)) * (mu1 - t) * (1 + theta1 * (1 + mu1))^2
  term1_denominator <- (1 + theta1 * (1 + mu1))^4
  
  term2_numerator <- (1 + theta1 * (1 + mu1))^2 * exp(-theta1 * (t - mu1)) * theta1 * (-2 * t - 2 * t^2 - 2 * t^2 * mu1 + 2 * mu1 + 2 * mu1^2 + 2 * mu1^2 * t) + (-2 * t - t^2 + 2 * mu1 + mu1^2)
  term2_denominator <- (1 + theta1 * (1 + mu1))^4
  
  term3_numerator <- 2 * (theta1 * (mu1^2 + 2 * mu1 - t^2 - 2 * t) + theta1^2 * (mu1 + mu1^2 + mu1^2 * t - t - t^2 - mu1* t^2)) * exp(-theta1 * (t - mu1)) * (1 + mu1) * (1 + theta1 * (1 + mu1))
  term3_denominator <- (1 + theta1 * (1 + mu1))^4
  
  result6 <- (term1_numerator / term1_denominator) + (term2_numerator / term2_denominator) - (term3_numerator / term3_denominator)
  
  return(result6)
}
result6 <- R22(t,theta1,mu1)
print(result6)
R22=result6
R22
#################################################
v1=R1*T11+R2*T12
v1
v2=R1*T12+R2*T22
v2
####Squared errorloss function 
###############################################################################
R_SE=R+0.5*(R11*T11+R22*T22+2*R12*T12)+eta1*v1+eta2*v2+0.5*(L111*T11*v1+L222*T22*v2+L112*(2*T12*v1+T11*v2))
R_SE  









#### Llnex loss function###################################################################################
nu=1################################ change##########################################################
h1_LE= -nu*exp(-nu*theta1)
h11_LE=(nu)^2*exp(-nu*theta1)

T_LE <- function(nu,theta1, h11_LE, T11, h1_LE,eta1,eta2, L111, L222, L112, T22,T12){
  result12 <- -1/nu * log(exp(-nu*theta1) + 0.5 * h11_LE * T11 + h1_LE * (T11 *eta1 + T21 * eta2) + 
                            0.5 * h1_LE * (L111 * T11^2 + L222 * T22 * T12 + 3 * L112 * T12 * T11))
  return(result12)
}
result12<-T_LE(nu,theta1, h11_LE, T11, h1_LE,eta1,eta2, L111, L222, L112, T22,T12)
print(result12)
theta_LE=result12
theta_LE
##################################
h2_LE= -nu*exp(-nu*mu1)
h22_LE=(nu)^2*exp(-nu*mu1)

m_LE <- function(nu,mu1, h22_LE, T11, h2_LE,eta1,eta2, L111, L222, L112, T22,T12){
  result22 <- -1/nu * log(exp(-nu*mu1) + 0.5 * h22_LE * T22 + h2_LE * (T22 *eta2 + T12* eta1) + 
                            0.5 * h2_LE * (L111 * T11*T12 + L222 * T22^2  + L112*(T11*T22+2*T12^2)))
  return(result22)
}
result22<-m_LE(nu,mu1, h22_LE, T11, h2_LE,eta1,eta2, L111, L222, L112, T22,T12)
print(result22)
mu_LE=result22
mu_LE
############################################
h_1_LE=-nu*exp(-nu*R)*R1
h_2_LE=-nu*exp(-nu*R)*R2
h_11_LE=nu*exp(-nu*R)*(nu*R1^2-R11)
h_22_LE=nu*exp(-nu*R)*(nu*R2^2-R22)
h_12_LE=nu*exp(-nu*R)*(nu*R1*R2-R12)

R_LE <- function(nu,R,h_11_LE, T11,h_22_LE,h_12_LE,eta1,eta2, v1,v2,L111, L222, L112, T22,T12){
  result33<- -1/nu * log(exp(-nu*R) + 
                           0.5 * (h_11_LE *T11 + h_22_LE*T22+2 * h_12_LE *T12) + 
                           v1 * eta1 + v2 * eta2 + 
                           0.5 * (L111 * T11 * v1 + L222 *T22 * v2 + L112 * (2 *T12 * v1 +T11 * v2))) 
  return(result33)
}
result33<-R_LE(nu,R,h_11_LE, T11,h_22_LE,h_12_LE,eta1,eta2, v1,v2,L111, L222, L112, T22,T12)
print(result33)
R_LE=result33
R_LE


##########################################################################################
###true survival function
R_t <- function(t, theta, mu) {
  result <- ((1 + theta * (1 + t)) / (1 + theta * (1 + mu))) * exp(-theta* (t - mu))
}
result_value_11<- R_t(t,theta,mu)
print(result_value_11)
R_t=result_value_11
R_t
#############################################################################################












###########################################################################################################
##General entropy loss function
w=1#######################################change###############################################
h1_GE= -w*theta1^(-w-1)
h11_GE=w*(w+1)*theta1^(-w-2)

T_GE <- function(w,theta1, h11_GE, T11, h1_GE,eta1,eta2, L111, L222, L112, T22,T12){
  result_12 <- (theta1^-w+ 0.5 * h11_GE * T11 + h1_GE * (T11 *eta1 + T21 * eta2) + 
                  0.5 * h1_GE * (L111 * T11^2 + L222 * T22 * T12 + 3 * L112 * T12 * T11))^(-1/w)
  return(result_12)
}
result_12<-T_GE(w,theta1, h11_GE, T11, h1_GE,eta1,eta2, L111, L222, L112, T22,T12)
print(result_12)
theta_GE=result_12
theta_GE
####################################
h2_GE=-w*mu1^(-w-1)
h22_GE=w*(w+1)*mu1^(-w-2)

m_GE <- function(w,mu1, h22_GE, T11, h2_GE,eta1,eta2, L111, L222, L112, T22,T12){
  result_22 <-(mu1^-w + 0.5 * h22_GE * T22 + h2_GE * (T22 *eta2 + T12* eta1) + 
                 0.5 * h2_GE * (L111 * T11*T12 + L222 * T22^2  + L112*(T11*T22+2*T12^2)))^(-1/w)
  return(result_22)
}
result_22<-m_GE(w,mu1, h22_GE, T11, h2_GE,eta1,eta2, L111, L222, L112, T22,T12)
print(result_22)
mu_GE=result_22
mu_GE
###############################
h_1_GE=-w*R^(-w-1)*R1
h_2_GE=-w*R^(-w-1)*R2
h_11_GE=w*(w+1)*R^(-w-2)*R1^2-w*R^(-w-1)*R11
h_22_GE=w*(w+1)*R^(-w-2)*R2^2-w*R^(-w-1)*R22
h_12_GE=w*(w+1)*R^(-w-2)*R1*R2-w*R^(-w-1)*R12

R_GE <- function(w,R,h_11_GE, T11,h_22_LE,h_12_LE,eta1,eta2, v1,v2,L111, L222, L112, T22,T12){
  result_333<- (R^-w+ 0.5 * (h_11_GE *T11 + h_22_GE*T22+2 * h_12_GE *T12) + 
                  v1 * eta1 + v2 * eta2 + 
                  0.5 * (L111 * T11 * v1 + L222 *T22 * v2 + L112 * (2 *T12 * v1 +T11 * v2)))^(-1/w) 
  return(result_333)
}
result_333<-R_GE(w,R,h_11_LE, T11,h_22_LE,h_12_LE,eta1,eta2, v1,v2,L111, L222, L112, T22,T12)
print(result_333)
R_GE=result_333
R_GE
############################################################################################
n1
#############################################################################################
####For MLE
theta1
mu1
R_t_MLE
#####For squared error loss function 
theta_SE
mu_SE
R_SE
############For Linex loss function 
nu
theta_LE
mu_LE
R_LE
#####FOR General entropy Loss function
w
theta_GE
mu_GE
R_GE



