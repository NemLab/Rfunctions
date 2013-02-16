#######################
### General Scripts ###
#######################


### read.Rdata
### Allows .Rdata objects to be read into single variables
### By default, 'load' loads a bunch of unwanted stuff, which this avoids

read.Rdata <- function(file1){
    get(load(file1))
}

### ODfunction
### third-order polynomial fit for OD660 vs yeast concentration
### data from UC Boulder info (same as in CSHL Methods in Yeast Genetics)
### Use format(y,sci=T) to understand the values better

ODfunction <- function(x) {
	y <-  -29893 + 14131618*x -6973847*x^2 + 11618382*x^3
	return(y)
}

### additiondilution
### returns how much culture/whatever you need to add to get
### the desired concentration with variable final volumes.
additiondilution <- function(c1,c2,v2) {
	x <- c2*v2/(c1-c2)
	return(x)
}

### addmediadilution
### returns how much culture/whatever you need to add to get
### the desired concentration with variable final volumes.
addmediadilution <- function(c1,c2,v1) {
    x <- (c1*v1-c2*v1)/c2
    return(x)
}

theme_clean <- function(base_size = 12, base_family = "") {
  library(grid) # sometimes required for 'unit' function
  theme(
    line =               element_line(colour = "black", size = 0.5, linetype = 1,
                            lineend = "butt"),
    rect =               element_rect(fill = "white", colour = "black", size = 0.5, linetype = 1),
    text =               element_text(family = base_family, face = "plain",
                            colour = "black", size = base_size,
                            hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 0.9),
    axis.text =          element_text(size = rel(0.8), colour = "black"),
    strip.text =         element_text(size = rel(0.8)),

    axis.line =          element_line(size = 0.5),
    axis.text.x =        element_text(vjust = 1),
    axis.text.y =        element_text(hjust = 1),
    axis.ticks =         element_line(colour = "black"),
    axis.title.x =       element_text(),
    axis.title.y =       element_text(angle = 90),
    axis.ticks.length =  unit(0.15, "cm"),
    axis.ticks.margin =  unit(0.1, "cm"),

    legend.background =  element_rect(colour = NA),
    legend.margin =      unit(0.2, "cm"),
    legend.key =         element_rect(fill = NA, colour = NA),
    legend.key.size =    unit(1.2, "lines"),
    legend.key.height =  NULL,
    legend.key.width =   NULL,
    legend.text =        element_text(size = rel(0.8)),
    legend.text.align =  NULL,
    legend.title =       element_text(size = rel(0.8), face = "bold", hjust = 0),
    legend.title.align = NULL,
    legend.position =    "right",
    legend.direction =   NULL,
    legend.justification = "center",
    legend.box =         NULL,

    panel.background =   element_blank(),
    panel.border =       element_blank(),
    panel.grid.major =   element_blank(),
    panel.grid.minor =   element_blank(),
    panel.margin =       unit(0.25, "lines"),

    strip.background =   element_rect(fill = NA, colour = NA),
    strip.text.x =       element_text(),
    strip.text.y =       element_text(angle = -90),

    plot.background =    element_rect(colour = NA),
    plot.title =         element_text(size = rel(1.2)),
    plot.margin =        unit(c(1, 1, 0.5, 0.5), "lines"),

    complete = TRUE
  )
}

theme_clean_original <- function(base_size = 12) {
 library(grid)
 structure(list(
    axis.line =         theme_segment(size=.5),
    axis.text.x =       theme_text(size = base_size * 0.8 , lineheight = 0.9, vjust = 1),
    axis.text.y =       theme_text(size = base_size * 0.8, lineheight = 0.9, hjust = 1),
    axis.ticks =        theme_segment(colour = "black", size = 0.4),
    axis.title.x =      theme_text(size = base_size, vjust = 0),
    axis.title.y =      theme_text(size = base_size, angle = 90, vjust = 0.5),
    axis.ticks.length = unit(0.3, "lines"),
    axis.ticks.margin = unit(0.5, "lines"),

    legend.background = theme_blank(),
    legend.key =        theme_blank(),
#    legend.key =        theme_rect(colour = "grey80"),
#    legend.key.size =   unit(6.2, "lines"),
    legend.key.size =   unit(1.2, "lines"),
    legend.text =       theme_text(size = base_size * 0.8),
    legend.title =      theme_text(size = base_size * 0.8, face = "bold", hjust = 0),
    legend.position =   "right",

    panel.background =  theme_blank(),
#    panel.background =  theme_rect(fill = "white", colour = NA),
#    panel.background =  theme_blank(),
    panel.border =      theme_blank(),
    panel.grid.major =  theme_blank(),
    panel.grid.minor =  theme_blank(),
#    panel.margin =      unit(0.3, "lines"),
    panel.margin =      unit(0.4, "lines"),

#    strip.background =  theme_rect(fill = "grey80",colour="grey80", size = 0.25),
    strip.background =  theme_rect(fill = 'NA',colour='NA', size = 0.25),
    strip.label =       function(variable, value) value,
    strip.text.x =      theme_text(size = base_size * 0.8),
    strip.text.y =      theme_text(size = base_size * 0.8, angle = -90),

    plot.background =   theme_rect(colour = NA),
    plot.title =        theme_text(size = base_size * 1.2),
    plot.margin =       unit(c(1, 1, 0.5, 0.5), "lines")
  ), class = "options")
}

###########################
### Third-Party Scripts ###
###########################

### arrange:
### Taken from
### http://gettinggeneticsdone.blogspot.com/2010/03/arrange-multiple-ggplot2-plots-in-same.html
### Plots multiple ggplot2 plots together.  ncol/nrow determines arrangement.
### Call as, e.g. arrange(plot1,plot2,ncol=2)

vp.layout <- function(x, y) viewport(layout.pos.row=x, layout.pos.col=y)
arrange <- function(..., nrow=NULL, ncol=NULL, as.table=FALSE) {
 dots <- list(...)
 n <- length(dots)
 if(is.null(nrow) & is.null(ncol)) { nrow = floor(n/2) ; ncol = ceiling(n/nrow)}
 if(is.null(nrow)) { nrow = ceiling(n/ncol)}
 if(is.null(ncol)) { ncol = ceiling(n/nrow)}
        ## NOTE see n2mfrow in grDevices for possible alternative
grid.newpage()
pushViewport(viewport(layout=grid.layout(nrow,ncol) ) )
 ii.p <- 1
 for(ii.row in seq(1, nrow)){
 ii.table.row <- ii.row
 if(as.table) {ii.table.row <- nrow - ii.table.row + 1}
  for(ii.col in seq(1, ncol)){
   ii.table <- ii.p
   if(ii.p > n) break
   print(dots[[ii.table]], vp=vp.layout(ii.table.row, ii.col))
   ii.p <- ii.p + 1
  }
 }
}


#########################
###  Modeling Scripts ###
#########################

### iaaregress:
### Does a logistic regression on a properly-formatted table: FL1.A first, Time second.
### Note: I haven't updated the function descriptions for the 'log' model
### in log model, x and e are replaced by log(x) and log(e)

iaaregress <- function(table,param=4) {

	if(param==4) {
	    # 4-parameter model
	    # f(x) = c + \frac{d-c}{(1+\exp(b(x - e)))}
	    regress0 <- drm(data=table, mean~time, fct = LL.4()) #log model
#	    regress0 <- drm(data=table, mean~time, fct = L.4()) # non log model
#  		regress0 <- drm(data=table, mean~time, fct = L.4(fixed=c(NA,min(table[,"mean"]),max(table[,"mean"]),NA)))
    } else {
		if(param==3) {
			# 3-parameter model
			# f(x) = c + \frac{d-c}{(1+\exp(b*x)}
			regress0 <- drm(data=table, mean~time, fct = L.3())
		}
        if(param==5) {
            # 5-parameter model
            # f(x) = c + \frac{d-c}{(1+\exp(b(x - e)))^f}
            # Also sometimes called the Boltzmann model
            regress0 <- drm(table, mean~time, fct = L.5())
#            regress0 <- drm(table, mean~time, fct = L.5(fixed=c(NA,min(table[,"mean"]),max(table[,"mean"]),NA,NA)))
        }
    }

    return(regress0)
}

### qplot.logistic:
### Makes plotting time series + fit + thalf a little easier
### Assumes that you have 'time' and 'mean' columns

qplot.logistic <- function(timeseriesdata,minval=0) {

	# Generate fit object
	fitobject <- iaaregress(timeseriesdata)

	# Calculate thalf
	if (minval!=0) {
		thalf <- predict.thalf(fitobject,minval=minval)
	} else {
		thalf <- predict.thalf(fitobject)
	}

	# Calculate text position

	textpos <- range(timeseriesdata[,"mean"])[1] + 0.75*diff(range(timeseriesdata[,"mean"]))

	# Plot
	a <- qplot(data=timeseriesdata,time,mean) +
	geom_line(data=predict.logistic(fitobject),color="blue") +
	geom_vline(xintercept=thalf,color="red") +
	geom_text(aes(x=thalf+0.1*max(timeseriesdata[,"time"]),y=textpos),label=paste("t½ =",signif(thalf,digits=3)))

	return(a)
}

### predict.logistic:
### Generates a data table based on a logistic fit object
### Right now, only works w/ variables "mean" and "time", 
### but could easily be changed to dynamically name

predict.logistic <- function(fitobject) {

	# Generate table for time vs predicted mean
	xrange <- range(fitobject$data[,1])
	predicted <- data.frame(x=seq(0,max(fitobject$data[,1]),length.out=100),y=NA)
	predicted[,2] <- predict(fitobject,newdata=predicted)

	colnames(predicted) <- colnames(fitobject$data)[1:2]

	return(predicted)
}

### thalf:
### Estimates t1/2 using a logistic fit model (drm)
### Currently only built for a 4-parameter model
### You can specify the 'min value' (e.g. from steady state data)
### by setting the minval option.

predict.thalf <- function(fitobject,minval=0) {
	coefficients <- fitobject$coefficients
	b <- coefficients[[1]]
	c <- coefficients[[2]]
	d <- coefficients[[3]]
	e <- coefficients[[4]]

	maxval <- fitobject$data[1,2] # First fluorescence data point
#	maxval <- c+(d-c)/(1+exp(b*-e)) # Only works when data is very logistic-ish
	if (minval == 0) {
		minval <- max(c(min(fitobject$data[,2]),c))
	}
	halfmax <- mean(c(maxval,minval))

	thalf <- exp(log( ( d - c )/( halfmax - c ) - 1 )/b + log(e)) # log model
#	thalf <- log( ( d - c )/( halfmax - c ) - 1 )/b + e # non-log model

	return(thalf)
}


### multiplegillespie:
### Uses GillespieSSA to simulate multiple realizations of a system
### Specify the number as n=.  Defaults to 10.  Returns data frame.
multiplegillespie <- function(x0,a,nu,parms,tf,n=10) {
	#Fencepost and loop
	run1 <- ssa(x0=x0,a=a,nu=nu,parms=parms,tf=tf,method="OTL")$data
	colnames(run1)[1] <- "time"
	run1 <- data.frame(run="run1",run1,row.names=NULL)
	totalruns <- run1

	for (i in 2:n) {
		nextrun <- ssa(x0=x0,a=a,nu=nu,parms=parms,tf=tf)$data
		colnames(nextrun)[1] <- "time"
		currentname <- paste("run",i,sep="")
		nextrun <- data.frame(run=currentname,nextrun,row.names=NULL)
		totalruns <- rbind(totalruns,nextrun)
	}

	return(totalruns)

}


### biggillespie:
### Generates random X values in a normal distribution
### with mean = X and sd = X/3 (similar to data)
biggillespie <- function(x0,a,nu,parms,tf,n=10) {
	Xvector <- rnorm(n,mean=x0["X"],sd=x0["X"]/3)

    #Fencepost and loop
	x0["X"] <- Xvector[1]
    run1 <- ssa(x0=x0,a=a,nu=nu,parms=parms,tf=tf,method="OTL")$data
    colnames(run1)[1] <- "time"
    run1 <- data.frame(run="run1",run1,row.names=NULL)
    totalruns <- run1


    for (i in 2:n) {
		x0["X"] <- Xvector[i]
        nextrun <- ssa(x0=x0,a=a,nu=nu,parms=parms,tf=tf)$data
        colnames(nextrun)[1] <- "time"
        currentname <- paste("run",i,sep="")
        nextrun <- data.frame(run=currentname,nextrun,row.names=NULL)
        totalruns <- rbind(totalruns,nextrun)
    }

    return(totalruns)

}


### runavg:
### Takes a data frame like that generated by biggillespie
### Pools the values together and splits into 100 evenly-spaced time intervals
### Calculates mean + sd at each interval and returns a data frame
### with time, mean, and sd of X.  Good for plotting.
runavg <- function(totalruns) {
	# Get timeframe
	maxtime <- max(totalruns[,"time"])
	timepoints <- seq(0,maxtime,maxtime/100)
	runningavg <- data.frame(time=timepoints[0:100],mean=matrix(100),sd=matrix(100))
	for (i in 1:100) {
		currentint <- subset(totalruns,time>=timepoints[i] & time < timepoints[i+1])
		intmean <- mean(currentint[,"X"])
		intsd <- sd(currentint[,"X"])
		runningavg[i,"mean"] <- intmean
		runningavg[i,"sd"] <- intsd
	}

	return(runningavg)
}

writeplot <- function(plotobject,name="plot",outdir="~/") {
	# write ggplot2 object to .Rdata, svg, and pdf
	save(plotobject,file=paste(outdir,name,".Rdata",sep=""))
	ggsave(plotobject,file=paste(outdir,name,".svg",sep=""))
	ggsave(plotobject,file=paste(outdir,name,".pdf",sep=""))
}

#w303_genotype
# I don't understand how, but the input is fuzzy for some reason, e.g. w303_genotype(l="test") returns leu2::test
# usage:    mat: sets mating type, accepting either "a" or "alpha"
#           markers: the markers seen below can also be set. With the exception of the term "wt", anything defined for a given
#                    locus is assumed to be inserted there, e.g. leu2="test" makes a genotype of leu2::test
#                    If the term "wt" or "WT" is used, it's assumed to be wild-type (e.g. the ade2+ locus is by default in W303)

w303_genotype <- function(
                    mat=NA,
                    leu2=NA,
                    trp1=NA,
                    can1=NA,
                    ura3=NA,
                    ade2=NA,
                    his3=NA) {

    loci <- c(leu2,trp1,can1,ura3,ade2,his3)
    loci[5] <- "wt"

    markers <- c("leu2","trp1","can1","ura3","ade2","his3")

    # leu2 marker
    leu2_d <- "leu2-3,112"
    # trp1 mraker
    trp1_d <- "trp1-1"
    # can1 marker
    can1_d <- "can1-100"
    # ura3 marker
    ura3_d <- "ura3-1"
    # ade2-1 marker
    ade2_d <- "ade2-1"
    # his3 marker
    his3_d <- "his3-11,15"

    loci_default <- c(leu2_d,trp1_d,can1_d,ura3_d,ade2_d,his3_d)

    for (i in 1:6) {
        if (is.na(loci[i])) {
            loci[i] <- loci_default[i]
        } else if (loci[i]=="wt"|loci[i]=="WT") {
            loci[i] <- toupper(markers[i])
        } else {
            loci[i] <- paste(markers[i],"::",loci[i],sep="")
        }
    }

    # Mating type can be either a or alpha
    if (is.na(mat)) {
        message("defaulting to MATa. To switch, set mat='alpha'")
        mat <- "MATa"
    } else if (mat=="a") {
        mat <- "MATa"
    } else if (mat=="alpha") {
        mat <- "MATα"
    } else {
        stop("mat must be set to 'a' or 'alpha'")
    }

    loci_string <- paste(loci,collapse=" ")
    genotype <- paste(mat,loci_string)

    return(genotype)
}

# yeast genotype maker.
# e.g. MATa leu2-3,112 trp1-1 can1-100 ura3-1 ade2-1 his3-11,15. delta sign?

yeast_genotype <- function() {
    print("Δ")
}

# color schemes

# degron color scheme - order is for IAA1, IAA6, IAA28, then IAA28.T2V
# Only works for four (or less?) things
degron_colors <- function() {scale_colour_manual(values=c("#4257A4","#6d426d","#269d9e","#73c59a"))}

# AFB2 + TIR1 + mTIR1 vs all IAAs plot
iaa_grid_colors <- function() {scale_colour_manual(values=c("#999999","#6699cc","#003399"))}
