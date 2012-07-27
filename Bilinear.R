### File for using the bilinear model/fit for Aux/IAA degradation (AFB+YFP-Aux/IAA system)
#
#
# model: 
#	x_dot = k1 - k2*x
#	y_dot = k3 - k4*y - k5*x*y

# Parameters from Shelly for IAA1 in TIR1 rep 1
parameters <- c(k1 = 0.0168,
		k2 = 0.0426,
		k3 = 0.31519,
		k4 = 0.0405,
		k5 = 0.0137)

# state variables. y is k3/k4. x is k1/k2. I don't know why x is so low, it should be closer to y.
state <- c(x = 0.0,
	   y = 7.78)

# better:
#state <- function(parameters) c(x=as.numeric(parameters[1]/parameters[2])/2,y=as.numeric(parameters[3]/parameters[4]))

# times and resolution
times <- seq(0,200,by=.1)

k3s <-  c(0.300, 0.310, 0.315, 0.320, 0.350)

#Bilinear <- function(t, parameters) {
Bilinear <- function(t, state, parameters) {
	with(as.list(c(state,parameters)), {
		# rate of change
		dx <- k1 - k2*x
		dy <- k3 - k4*y - k5*x*y

		# return rate of change
		list(c(dx,dy))
	})  #end with(as.list ...
}

# Make a bunch of cool stuff (varying k3)

bilinear_demo <- function() {
	library(deSolve)
#	k3s <-  c(0.300, 0.310, 0.315, 0.320, 5)
	k3s <-  c(0.300, 0.310, 0.315, 0.320, 0.350)
	parameters <- c(k1 = 0.0168,
			k2 = 0.0426,
			k3 = 0.31519,
			k4 = 0.0405,
			k5 = 0.0137)
	state <- c(x = 0.0,
		   y = 7.78)
	times <- seq(0,200,by=1)

	out_all <- lapply(k3s,
			function(x) {
				(ode(y=c(state[1],y=as.numeric(x/parameters[4])),
				     times=times,
				     func=Bilinear,
				     parms=c(parameters[1:2],
				     k3=as.numeric(x),
				     parameters[4:5]),
				     method="ode45"))
			})
	out_new <- data.frame(out_all[[1]])
	out_new$k3 <- k3s[1]

	for (i in 2:5) {
		current <- data.frame(out_all[[i]])
		current$k3 <- k3s[i]
		out_new <- rbind(out_new,current)
	}

	qplot(data=out_new,x=time,y=y,color=factor(k3),geom="line") + theme_clean() + scale_color_discrete(name="k3")
	return(out_new)
}

bilinear_plot <- function(x) {
        qplot(data=x,x=time,y=y,color=factor(k3),geom="line") + theme_clean() + scale_color_discrete(name="k3") +
                geom_line(data=x,aes(x=time,y=x),colour=I("black"))
}

