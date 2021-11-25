rm(list=ls())

#Some packages
library('drc')


#DSS function for viability data
DSS_Viability <- function (rdata, t = 10, DSS.type = 1, Max, Min, log = FALSE, concn.scale = 1e-9, fit = "NPI")
{
  if (t <= 0) {
    stop("Percent inhibition must be greater than Zero. Please, see manual for details.")
  }
  DSS <- rep(0, nrow(rdata))
  for (row in 1:nrow(rdata)) {
    temp <- as.matrix(rdata[row, ])
    if (isTRUE(ncol(rdata) == 5)) {
      c = as.numeric(as.character(str_replace(temp[1], ", ", ".")))
      b = as.numeric(as.character(str_replace(temp[2], ", ", ".")))
      d = 0
      a = as.numeric(as.character(str_replace(temp[3], ", ", ".")))
      Min.Conc <- as.numeric(as.character(str_replace(temp[4], ", ", ".")))
      Max.Conc <- as.numeric(as.character(str_replace(temp[5], ", ", ".")))
      if(fit == "NPI"){
        max <- 100
      }else{
        max <- Max[row]
      }
      y <- max*t/100
    }
    if (isTRUE(ncol(rdata) == 6)) {
      c = as.numeric(as.character(str_replace(temp[1], ", ", ".")))
      b = as.numeric(as.character(str_replace(temp[2], ", ", ".")))
      d = as.numeric(as.character(str_replace(temp[3], ", ", ".")))
      a = as.numeric(as.character(str_replace(temp[4], ", ", ".")))
      Min.Conc_raw <- as.numeric(as.character(str_replace(temp[5], ", ", ".")))
      Max.Conc_raw <- as.numeric(as.character(str_replace(temp[6], ", ", ".")))
      
      if(fit == "NPI"){
        max <- 100
        min <- 0
      }else{
        max <- Max[row]
        min <- Min[row]
      }
      y = (max-min)*t/100 + min
    }
    if (log == FALSE) {
      Min.Conc <- log(Min.Conc_raw * concn.scale)
      Max.Conc <- log(Max.Conc_raw * concn.scale)
      c <- log(c * concn.scale)
    } else {
      Min.Conc <- Min.Conc_raw
      Max.Conc <- Max.Conc_raw
    }
    if (is.na(c) || is.na(b) || is.na(a)) {
      DSS[row] <- as.numeric(NA)
      #} else if (isTRUE(c >= Max.Conc)) { DSS[row] <- 0 }
    }else if (b < 0) {#R equation of LL.4: b>0 means decreasing, but we are working in inhibition after normalisation
      DSS[row] <- 0
    } else {
      if (a > y) {
        if (y > d) {
          x <- (log((a-y)/(y-d))/b)+c
          #for y = d + ((a-d)/(1+exp(b*(x-c))), where x is log(x) and c is log(c), which is what drm uses
          
          if (isTRUE(x > Max.Conc)) {
            x <- Max.Conc
          } else if (isTRUE(x < Min.Conc)) {
            x <- Min.Conc
          }
        } else {
          x <- Max.Conc
        }
        
        int_y <- (((d-a)*log(exp(b*(x-c)) + 1))/b + (a*x)) - (((d-a)*log(exp(b*(Min.Conc-c)) + 1))/b + (a*Min.Conc))
        if(int_y < 0) {
          int_y <- 0
        }
        
        if(a > max){
          x_max <- (log((a-max)/(max-d))/b)+c
          if (isTRUE(x_max > Max.Conc)) {
            x_max <- Max.Conc
          } else if (isTRUE(x_max < Min.Conc)) {
            x_max <- Min.Conc
          }
          extra_area <- (((d-a)*log(exp(b*(x_max-c)) + 1))/b + (a*x_max)) - (((d-a)*log(exp(b*(Min.Conc-c)) + 1))/b + (a*Min.Conc)) - (x_max - Min.Conc)*100
          int_y <- int_y - extra_area
        }
        threshold_area <- (x - Min.Conc) * y
        total_area <- (x - Min.Conc) * (max - y)
        if (DSS.type == 1) {
          norm_area <- ((int_y - threshold_area)/(total_area)) * 100
        }
        if (DSS.type == 2) {
          norm_area <- ((int_y - threshold_area)/(total_area)) * 100 * (log10(max)/log10(a))
        }
        if (DSS.type == 3) {
          norm_area <- ((int_y - threshold_area)/(total_area)) * 100 * ((x - Min.Conc)/(Max.Conc - Min.Conc))
        }
        DSS[row] <- round(100-norm_area, digits = 8)
      }else {
        DSS[row] <- 0
      }
    }
  }
  return(matrix(as.numeric(DSS)))
}





########
# Simulate viability data from a logistic model
#number of replicates
N <- 5 
#Fix concentrations (log10)
n1 <- 10
x1 <- 10^seq(-5,3,length = n1)
n2 <- 6
x2 <- 10^seq(-2,3,length = n2)

#Sample viability by using as mean a logistic model with dose effects
b1 <- 0.75
b2 <- 0.5
b3 <- -0.5
y_mean <- (1 + exp(b1 * matrix(log10(x1), n1, n2) + b2 * matrix(log10(x2), n1, n2, byrow = TRUE)))^(-1)
#plot of means
matplot(log10(x1), y_mean, type = "l", lty = 1, main = "Drug 1", ylab = "Viability", xlab = "Concentration (log10)")      
matplot(log10(x2), t(y_mean), type = "l", lty = 1, main = "Drug 2", ylab = "Viability", xlab = "Concentration (log10)")      

#Add noise simulating the N replicates
set.seed(123)
y <- array(NA, dim = c(N, n1, n2))
for(i in 1:N){
  y[i,,] <- y_mean + matrix(rnorm(n1*n2, 0, sqrt(0.001)), n1, n2)
}
y <- y * 100
#Notice that we have a couple of points outside the range (0,1)
range(y)

#Compute DSS for each monotherapy (by fixing the dose fo the second drug)
DSS1 <- rep(NA, n1)
for(i1 in 1:n1){
  fit_now <- drm(c(t(y[,i1,])) ~ rep(x2, N), fct = LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "IC50")), control = drmc(errorm = FALSE))
  DSS1[i1] <- DSS_Viability(rdata = matrix(c(fit_now$coefficients[c(4,1,2,3)], min(fit_now$data[,1]), max(fit_now$data[,1])), nrow = 1), Min = fit_now$coefficients[2], Max = fit_now$coefficients[3], concn.scale=1e-9, DSS.type = 3, fit = "NPI")      
}
DSS2 <- rep(NA, n2)
for(i2 in 1:n2){
  fit_now <- drm(c(t(y[,,i2])) ~ rep(x1, N), fct = LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "IC50")), control = drmc(errorm = FALSE))
  DSS2[i2] <- DSS_Viability(rdata = matrix(c(fit_now$coefficients[c(4,1,2,3)], min(fit_now$data[,1]), max(fit_now$data[,1])), nrow = 1), Min = fit_now$coefficients[2], Max = fit_now$coefficients[3], concn.scale=1e-9, DSS.type = 3, fit = "NPI")      
}
#We can see and increment in the efficacy, especially for the second drug


##########Evaluate inteaction########
#Now compute the volume under the surface for this combination for each replicate
library("geometry")

getVolume=function(df){
  #find triangular tesselation of (x,y) grid
  res=delaunayn(as.matrix(df[,-3]),full=TRUE,options="Qz")
  #calulates sum of truncated prism volumes
  sum(mapply(function(triPoints,A) A/3*sum(df[triPoints,"z"]), split.data.frame(res$tri,seq_along(res$areas)), res$areas))
}

VUS <- rep(NA, N)
for(i in 1:N){
  df1 <- data.frame( x = c(matrix(x1, n1, n2)),
                     y = c(matrix(x2, n1, n2, byrow = TRUE)),
                     z = c(y[i,,]))
  
  VUS[i] <- getVolume(df1)/diff(range(df1$x))/diff(range(df1$y))
}

plot(VUS, xlab = "Replicate", main = "Volume Under the Surface", ylab = "")



