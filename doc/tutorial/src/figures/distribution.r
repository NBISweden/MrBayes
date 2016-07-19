library(msm)
library(gtools)

## gamma distributions
curve(dgamma(x, 0.1, 0.1), from=0, to=3, col="yellow", xlab="rate (r)", ylab="gamma density")
curve(dgamma(x, 0.5, 0.5), from=0, to=3, col="orange", add=T)
curve(dgamma(x, 1.0, 1.0), from=0, to=3, col="red", add=T)
curve(dgamma(x, 2.0, 2.0), from=0, to=3, col="blue", add=T)
curve(dgamma(x, 10., 10.), from=0, to=3, col="purple", add=T)
legend(0.5, 0.2, "a = b = 0.1",bty="n")
legend(0.5, .33, "a = b = 0.5",bty="n")
legend(0.5, 0.5, "a = b = 1",  bty="n")
legend(0.5, 0.7, "a = b = 2",  bty="n")
legend(0.5, 1.5, "a = b = 5",  bty="n")
 
## clock rate priors
curve(dlnorm(x, -7, 0.6), from=0, to=0.005, col="purple", xlab="clock rate", ylab="prob. density")
curve(dgamma(x, 2, 2000), from=0, to=0.005, col="red", add=T)
curve(dtnorm(x, 0.001, 0.0007, lower=0), from=0, to=0.005, col="blue", add=T)
curve(dexp(x, 1000), from=0, to=0.005, col="orange", add=T)
legend(0.002, 800, "lognormal(-7, 0.6)", lty=1, col="purple", bty="n")
legend(0.002, 750, "gamma(2, 2000)", lty=1, col="red", bty="n")
legend(0.002, 700, "normal(0.001, 0.0007)", lty=1, col="blue", bty="n")
legend(0.002, 650, "exp(1000)", lty=1, col="orange", bty="n")

## tree age priors
curve(dtnorm(x, 390, 60, lower=300), from=300, to=600, ylim=c(0,0.012), col="blue", xlab="age", ylab="prob. density")
curve(dexp(x-300, rate=1/(390-300)), from=300, to=600, col="red", add=T)
curve(dgamma(x-300, 2, 2/(390-300)), from=300, to=600, col="orange", add=T)
curve(dlnorm(x-300, 4.0, 1.0), from=300, to=600, col="purple", add=T)
legend(400, 0.011, "offsetlognormal(300, 4, 1)",lty=1, col="purple", bty="n")
legend(400, 0.010, "offsetgamma(300, 2, 0.022)", lty=1, col="orange", bty="n")
legend(400, 0.009, "offsetexp(300, 0.01111)", lty=1, col="red", bty="n")
legend(400, 0.008, "truncatednormal(300, 390, 60)", lty=1, col="blue", bty="n")

## node calibrations
curve(dgamma(x-100, 4, 0.08), from=100, to=600, col="purple", xlab="age", ylab="prob. density")
curve(dtnorm(x, 175, 25, lower=140), from=100, to=600, col="blue", add=T)
curve(dexp(x-300, rate=1/90), from=100, to=600, col="red", add=T)
legend(240, 0.017, "offsetgamma(100, 4, 0.08)", lty=1, col="purple", bty="n")
legend(260, 0.015, "truncatednormal(140, 175, 25)", lty=1, col="blue", bty="n")
legend(280, 0.013, "offsetexp(300, 0.0111)", lty=1, col="red", bty="n")

## for lognormal
u <- 4
s <- 1
mean <- exp(u + s^2/2)
median <- exp(u)
mode <- exp(u - s^2)
sd <- sqrt((exp(s^2) - 1) * exp(2*u + s^2))

## gamma-dirichlet and exp priors
par(mfrow=c(2,2))
# par(mar=c(4, 4, 2, 0.5))
curve(dgamma(x, 17, 10), from=0, to=5, col="red",
      xlab="tree length", ylab="gamma density", main="i.i.d. exp(10)")
legend(2.2, 0.6, "gamma(17, 10)", bty="n")
curve(dgamma(x, 1., 1.), from=0, to=5, col="red",
      xlab="tree length", ylab="", main="gamma-Dirichlet(1,1,1,1)")
legend(2.2, 0.6, "gamma(1, 1)", bty="n")
curve(dgamma(x, 1, 10), from=0, to=1, col="blue", ylim=c(0,5), 
      xlab="branch length", ylab="prob. density")
legend(0.4, 3, "exp(10)", bty="n")
v <- rgamma(1000000, 1, 1) * rdirichlet(1000000, rep(1,17))[,1]
plot(density(v, adjust=2, from=0), col="blue", xlim=c(0,1), ylim=c(0,5),
     xlab="branch length", ylab="", main="")
