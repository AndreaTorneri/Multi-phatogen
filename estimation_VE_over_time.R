##################################################################
### Fitting a Gamma distribution to vaccine effectiveness data ###
##################################################################

library(minpack.lm)


####### Influenza ####### 

# Reference: https://academic.oup.com/ofid/article/12/Supplement_1/ofae631.1510/7988530

VE_influenza <- data.frame(
  x = c(21, 45, 75, 105),
  y = c(0.59, 0.52, 0.39, 0.34)
)

VE_waning_influenza <- nlsLM(y ~ a * dgamma(x = x, 
                               shape = s1,
                               scale = sc1),
                             start = list(a = 10, s1 = 1.5, sc1 = 1.5),
                             data = VE_influenza)

parameters <- VE_waning_influenza$m$getPars()
print(parameters)

plot(seq(0,120,1),
     parameters[1]*dgamma(x=seq(0,120,1),
                          shape = parameters[2],
                          scale = parameters[3]),
     ylim = c(0,1))

points(VE_influenza, col="red")


VE_over_time_influenza <- function(t){
  return(parameters[1]*dgamma(t,shape = prms[2],scale = prms[3]))
}


t <- 20
94.748176 * dgamma(t, shape = 1.038451, scale = 129.382569)





