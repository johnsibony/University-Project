set.seed(26) #to repreduce our result

# Load the library

library(strucchange) #library CUSUM and Chow test
library(lmtest) #library for homoscedasticity and autocorrelation test

# Import the dataset

data<read.table(file="SBData.csv",sep=",",row.names=1,header=TRUE) #load all of the data
attach(data) #direct integration of the variable
sub1_data=data[is.element(data$Year, seq(1967,1988)), ] #keep only the data matching those of the paper
sub2_data=data[is.element(data$Year, seq(1989,2013)), ] #keep the remaining data (not in the paper)

# Run the OLS's

ols_all <- lm(SandP~Dum, data=data) #OLS for all of the data
residuals_all=residuals(ols_all)
ols_sub1 <- lm(SandP~Dum, data=sub1_data) #OLS for the data in the paper
residuals1=residuals(ols_sub1)
ols_sub2 <- lm(SandP~Dum, data=sub2_data) #OLS for the remaining data
residuals2=residuals(ols_sub2)

# First OLS (data matching the paper) 

summary(ols_sub1) #main information on the OLS
# Looking for homoscedasticity
par(mfrow = c(2, 2)) #matrix plot 
plot(fitted(ols_sub1),residuals1^2, xlab="Ŷ", ylab="residuals^2", pch=19, main = "Test for Homoskedasticity") #plot Ŷ against residuals squared
plot(residuals1^2, xlab="time", ylab="residuals", pch=19, main = "Test for Homoskedasticity") #plot residuals against time
bptest(ols_sub1, studentize=FALSE) #Breusch-Pagan test
bptest(ols_sub1, studentize=TRUE) #Breusch-Pagan robust test
# Looking for autocorrelation
plot(residuals1[2:length(residuals1)], residuals1[1:length(residuals1)-1], xlab="e(t-1)", ylab="e(t)", pch=19) #plot the correlation of the residuals of order 1
title("Test for Autocorrelation")
dwtest(ols_sub1) #Durbin Watson test
bgtest(ols_sub1, order=5, type="Chisq") #Breusch-Godfrey test
# Looking for Normality of the residuals
hist(residuals1, xlab="residuals", main="Test for Normality") #histogramme of the residuals
shapiro.test(residuals1) # Normality test
# Looking for specification
resettest(ols_sub1) #RESET test
# Looking for parameter stability
sctest(SandP[is.element(data$Year, seq(1967,1988))]~Dum[is.element(data$Year, seq(1967,1988))], type="Rec-CUSUM") #CUSUM test

# Second OLS (data out of the paper) 

summary(ols_sub2) #main information on the OLS
# Looking for homoscedasticity
par(mfrow = c(2, 2)) #matrix plot 
plot(fitted(ols_sub2),residuals2^2, xlab="Ŷ", ylab="residuals^2", pch=19, main = "Test for Homoskedasticity") #plot Ŷ against residuals squared
plot(residuals2^2, xlab="time", ylab="residuals", pch = 19, main = "Test for Homoskedasticity") #plot residuals against time
bptest(ols_sub2, studentize=FALSE) #Breusch-Pagan test
bptest(ols_sub2, studentize=TRUE) #Breusch-Pagan robust test
# Looking for autocorrelation
plot(residuals2[2:length(residuals2)], residuals2[1:length(residuals2)-1], xlab="e(t-1)", ylab="e(t)", pch=19) #plot the correlation of the residuals of order 1
title("Test for Autocorrelation")
dwtest(ols_sub2) #Durbin Watson test
bgtest(ols_sub2, order=5, type="Chisq") #Breusch-Godfrey test
# Looking for Normality of the residuals
hist(residuals2, xlab="residuals", main="Test for Normality") #histogramme of the residuals
shapiro.test(residuals2) #Normality test
# Looking for specification
resettest(ols_sub2) #RESET test
# Looking for parameter stability
sctest(SandP[is.element(data$Year, seq(1989,2013))]~Dum[is.element(data$Year, seq(1989,2013))], type="Rec-CUSUM") #CUSUM test

# General OLS (all of the data) 

summary(ols_all) #main information on the OLS
# Looking for homoscedasticity
par(mfrow = c(2, 2)) #matrix plot 
plot(fitted(ols_all),residuals_all^2, xlab="Ŷ", ylab="residuals^2", pch=19, main = "Test for Homoskedasticity") #plot Ŷ against residuals squared
plot(residuals_all^2, xlab="time", ylab="residuals", pch=19, main = "Test for Homoskedasticity") #plot residuals against time
bptest(ols_all, studentize=FALSE) #Breusch-Pagan test
bptest(ols_all, studentize=TRUE) #Breusch-Pagan robust test
# Looking for autocorrelation
plot(residuals_all[2:length(residuals_all)], residuals_all[1:length(residuals_all)-1], xlab="e(t-1)", ylab="e(t)", pch=19) #plot the correlation of the residuals of order 1
title("Test for Autocorrelation")
dwtest(ols_all) #Durbin Watson test
bgtest(ols_all, order=5, type="Chisq") #Breusch-Godfrey test
# Looking for Normality of the residuals
hist(residuals_all, xlab="residuals", main="Test for Normality") #histogramme of the residuals
shapiro.test(residuals_all) # Normality test
# Looking for specification
resettest(ols_all) #RESET test
# Looking for parameter stability
sctest(SandP~Dum, type="Rec-CUSUM") #CUSUM test
sctest(SandP~Dum, type="Chow", point=22) #Chow test at point 25

# General inspection of our dataset splitted into 2 period  

plot(sub1_data$Dum, sub1_data$SandP, xlab="Dum", ylab="SandP", col="green", pch=19) #first period
points(sub2_data$Dum, sub2_data$SandP, col="red", pch = 19) #second period
legend("center", legend=c("1967-1988", "1989-2013"), col=c("green", "red"), pch=19)

# Possible representation of our dataset if X was a continuous variable

Y=0.15*c(0:21)+rnorm(22, mean=0, sd=0.8)
Y=c(Y, 3 + 0*c(22:46)+rnorm(25, mean=0, sd=1.5))

my_ols_all=lm(Y~c(0:46))
my_ols_sub1=lm(Y[1:22]~c(0:21))
my_ols_sub2=lm(Y[23:47]~c(22:46))

plot(Y, xlab="X", pch=19)
lines(c(0:46), fitted(my_ols_all), col="green", lty=1, lwd=4)
lines(c(0:21), fitted(my_ols_sub1), col="blue", lty=1, lwd=4)
lines(c(22:46), fitted(my_ols_sub2), col="red", lty=1, lwd=4)
legend("topleft", legend=c("1967-1988", "1989-2013", "1967-2013"), col=c("blue", "red", "green"), lty=1)

bptest(my_ols_all, studentize=TRUE)
bptest(my_ols_all, studentize=FALSE)
sctest(Y~c(0:46), type="Chow", point=22)
summary(my_ols_sub1)$adj.r.squared
summary(my_ols_sub2)$adj.r.squared
summary(my_ols_all)$adj.r.squared
sctest(Y~c(0:46), type="Chow", point=22)
