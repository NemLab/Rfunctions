
### Usage notes for these scripts:
# Almost every function or object depends on the flowCore (and often flowViz) package(s) from biocondcutor.
# To install these, run these commands in R:
# 	source("http://bioconductor.org/biocLite.R")
# 	biocLite("flowCore")
# 	biocLite("flowViz")
#
# The most useful script for an R newbie is 'summary.cyt'.
# It will take a flowSet (see flowCore documentation) and run some QA, gate, and get summary FL1.A (or other channel) data
# as well as experiment information like the time of day and events/µL concentration.

library(flowViz)

#########################
###  Cytometer Gates  ###
#########################

##################################
### Notes on cytometer changes ###
##################################

### Got new cytometer around 2011-02-22-ish (updated gates)
### Switched FL1-A and FL2-A sometime around 02-23.  Was switched back relatively soon afterwards.
### FSC-A was tweaked around 2011-03-8 (early-ish March), changed gates.

### yeastGate
### Defines an SSC.A vs FSC.A gate.  Includes only the yeast population
### from a flowSet

	# As of 2011-10-26, it may be necessary to create a new yeastGate for diploids that truly cuts out dead cells.

	# Used from 2011-02-22 to present, excludes debris/non-yeast
	yeastGate <<- polygonGate(filterId="Yeast",
	    .gate=matrix(c(400000,10000000,10000000,400000, 10000,10000,2300000,60000),
	    ncol=2,nrow=4,dimnames=list(c("1","1","1","1"),c("FSC.A","SSC.A"))))

	# Used until cytometer was switched around 2011-02-22
	oldyeastGate <<- polygonGate(filterId="Yeast",
		.gate=matrix(c(160000,1500000,1500000,160000, 0,0,200000,200000),
	    ncol=2,nrow=4,dimnames=list(c("1","1","1","1"),c("FSC.A","SSC.A"))))


	# Diploid Singlet gate, used 2011-07-09 to present
	# Diploids are slightly larger and have better separation between singlets/doublets
	dipsingletGate <<- polygonGate(filterId="DipSingletGate",
		.gate=matrix(c(
			#x values
			7.5e5,13e5,18e5,15e5,6e5,
			#y values
			9e5,16e5,26e5,30e5,15e5),
			ncol=2,nrow=5,dimnames=list(rep(NA,5),c("FSC.A","FSC.H"))
			      )
			     )

	# Diploid Doublet gate, used 2011-08-09 to present
	# Diploids are slightly larger and have better separation between singlets/doublets
	dipdoubletGate <<- polygonGate(filterId="DipDoubletGate",
		.gate=matrix(c(
			#x values
			10e5,17e5,23e5,22e5,20e5,8e5,
			#y values
			8e5,12.5e5,17e5,20e5,22e5,8.5e5),
			ncol=2,nrow=6,dimnames=list(rep(NA,6),c("FSC.A","FSC.H"))
			      )
			     )


	# Diploid Doublet gate, used 2011-08-09 to present
	# Diploids are slightly larger and have better separation between singlets/doublets
	hapdoubletGate <<- polygonGate(filterId="HaploidDoubletGate",
		.gate=matrix(c(
			#x values
			6.5e5,1.15e6,1.5e6,1.4e6,1.2e6,5e5,
			#y values
			5.75e5,9e5,1.3e6,1.4e6,1.5e6,6.5e5),
			ncol=2,nrow=6,dimnames=list(rep(NA,6),c("FSC.A","FSC.H"))
			      )
			     )

	# Used 2012-02-22 to present
	hapsingletGate <<- polygonGate(filterId="HaploidSingletGate",
		.gate=matrix(c(
			#x values
			5e5,0.8e6,1.15e6,1e6,5e5,
			#y values
			8e5,1.05e6,1.5e6,1.8e6,1e6),
			ncol=2,nrow=5,dimnames=list(rep(NA,5),c("FSC.A","FSC.H"))
			      )
			     )

    ####################
    # Havens et al 2012#
    ####################
	# Used in auxin paper (Havens 2012) to gate all yeast from non-yeast. Also excludes a portion
	# of small-FSC.A, high-SSC.A cells (presumably dead).
	auxinpaper_yeastGate <<- polygonGate(filterId="Yeast",
	    .gate=matrix(c(400000,10000000,10000000,400000, 10000,10000,2300000,60000),
	    ncol=2,nrow=4,dimnames=list(c("1","1","1","1"),c("FSC.A","SSC.A"))))
	# Used in auxin paper to gate for singlets after gating with auxinpaper_yeastGate
	auxinpaper_singletGate <<- polygonGate(filterId="DipSingletGate",
		.gate=matrix(c(
			#x values
			7.5e5,13e5,18e5,15e5,6e5,
			#y values
			9e5,16e5,26e5,30e5,15e5),
			ncol=2,nrow=5,dimnames=list(rep(NA,5),c("FSC.A","FSC.H"))
			      )
			     )
	# Used in auxin paper to gate for doublets after gating with auxinpaper_yeastGate
	auxinpaper_doubletGate <<- polygonGate(filterId="DipDoubletGate",
		.gate=matrix(c(
			#x values
			10e5,17e5,23e5,22e5,20e5,8e5,
			#y values
			8e5,12.5e5,17e5,20e5,22e5,8.5e5),
			ncol=2,nrow=6,dimnames=list(rep(NA,6),c("FSC.A","FSC.H"))
			      )
			     )

	# haploid singlet gate, used from 2011-02-22 to 2012-02-22
	singletGate2 <<- polygonGate(filterId="Singlets",
		.gate=matrix(c(400000,3000000,3000000,160000, 620000,3500000,6000000,500000),
		ncol=2,nrow=4,dimnames=list(c("1","1","1","1"),c("FSC.A","FSC.H"))))



	# Used until cytometer was switched around 2011-02-22
	oldsingletGate <<- polygonGate(filterId="Singlets",
	    .gate=matrix(c(160000,1500000,1400000,160000, 680000,5700000,6000000,750000),
	    ncol=2,nrow=4,dimnames=list(c("1","1","1","1"),c("FSC.A","FSC.H"))))


	# Used up to 2011-02-10
	extraoldsingletGate <<- polygonGate(filterId="Singlets",
	    .gate=matrix(c(160000,800000,800000,160000, 680000,3150000,3500000,750000),
	    ncol=2,nrow=4,dimnames=list(c("1","1","1","1"),c("FSC.A","FSC.H"))))

## experimental gates

    # mChGate - designed around data from cytometer experiment Nick did on 2012-9-19 (NB036)
    # subsets cells with high FL3.A (mCherry) vs FL1.A (EYFP) output.
    # Note that it was designed around cells that have both EYFP and mCherry and
    # may not work for mCherry-only cells
    mChGate <<- polygonGate(filterId='mCherryGate',
                            matrix(c(5e2,2e3,4e4,7e4,8e3,5e2, 2e3,2e3,1.8e4,9e4,8e4,4e3),
                                ncol=2,
                                nrow=6,
                                dimname=list(c(1,1,1,1,1,1),c("FL1.A","FL3.A")))
                            )

    # Same motivatio as the mChGate - separate EYFP+mCh cells from EYFP cells. In this case,
    # subsets for EYFP-only cells. Again, may only work for this particular experiment type
    EYFPGate <<- polygonGate(filterId='EYFPGate',
                             matrix(c(0,4e3,3e4,3e4,3e3,0, 0,0,4e3,9e3,2e3,2e3),
                                ncol=2,
                                nrow=6,
                                dimname=list(c(1,1,1,1,1,1),c("FL1.A","FL3.A"))
                                )
                            )

	# Diploid Singlet gate, first created 2011-07-09
	# Diploids are slightly larger and have better separation between singlets/doublets
	dipsingletGate3 <<- polygonGate(filterId="DipSingletGate2",
		.gate=matrix(c(
			#x values
			5e5,13e5,18e5,15e5,6e5,
			#y values
			7e5,16e5,26e5,30e5,15e5),
			ncol=2,nrow=5,dimnames=list(rep(NA,5),c("FSC.A","FSC.H"))
			      )
			     )

	# for diploids, highlights dead cells
	deadGate <<- polygonGate(filterId="deadGate",
		.gate=matrix(c(7.5e5,13e5,3e6,15e5,6e5,
				1e5,2e5,7.5e5,4e5,2e5),
				ncol=2,nrow=5,dimnames=list(rep(NA,5),c("FSC.A","SSC.A")))
			)

	# for diploids, takes small-sized cells, maybe could use for excluding dead cells when combined with clustering
	deadexcludeGate <<- polygonGate(filterId="deadGate",
                .gate=matrix(c(7.5e5,4e6,3e6,1.5e6,6e5,
				0,0,7.5e5,8e5,2e5),
                                ncol=2,nrow=5,dimnames=list(rep(NA,5),c("FSC.A","SSC.A")))
                        )

	# for diploids, is a little more generous than 'dipsingletGate' Combine with flowClust to clean it up
	dipsingletGate2 <<- polygonGate(filterId="DipSingletGate",
		.gate=matrix(c(
			7.5e5,14e5,30e5,15e5,6e5,
			9e5,16e5,35e5,40e5,15e5),
			ncol=2,nrow=5,dimnames=list(rep(NA,5),c("FSC.A","FSC.H"))))

	# excludes the big-cell subpopulation, it might be sick/weird/less comparable
	# use !excludeBig to exclude them.
	excludeBig <<- polygonGate(filterId="excludeBig",
		.gate=matrix(c(
                        3e6,2e7,2e7,2e7,3e6,
                        1e6,16e5,5e6,1e7,1e7),
                        ncol=2,nrow=5,dimnames=list(rep(NA,5),c("FSC.A","FSC.H"))))

    # Rob's gate for E. coli
	ecoliGate <<- polygonGate(filterId="E.coli",
	    .gate=matrix(c(2e4,8e4,8e4,2e4,  # x points
                       1,1,6e3,6e3),     # y points
	    ncol=2,nrow=4,dimnames=list(c("1","1","1","1"),c("FSC.H","SSC.H"))))


###########################
###  Cytometer Scripts  ###
###########################

### polygate:
### Make a gate easier?

polygate <- function(x,y,filterID="newGate",channels=c("FSC.A","FSC.H")) {
    if( length(x) != length(y) | !is.numeric(x) | !is.numeric(y)) {
        stop("x coordinate vector must be same length as y coordinate vector")
    }

    gate <- polygonGate(filterId=filterID,
                        .gate=matrix(c(x,y),
                        ncol=2,nrow=length(x),dimnames=list(rep(NA,5),channels)))
    return(gate)
}

### ploidy:
### Tries to guess the ploidy of a given flowframe
### Uses FSC.A/FSC.H ratio.
### Diploids are typically 5um x 6um ellipsoids while
### haploids are typically 4um x 4um spheroids
### As a result, diploids are longer and you get a larger 'area/volume' FSC.A
### 'Width' might also be useful.

ploidy <- function(flowframe) {
	# Find FSC.A/FSC.H.  This is close to 1 for diploids and close to .8 for haploids
	# Test this assumption!!!!!
	fsca <- summary(flowframe)[4,1]
	fsch <- summary(flowframe)[4,7]
	quotient <- fsca/fsch
	if(quotient>0.92) {
		return(c("Diploid",quotient))
	} else{
		return(c("Haploid",quotient))
	}
}

### qa.gating:
### Very simple script to check whether a flowSet or flowFrame
### contains empty values, in which case normalization may fail (divide by zero)

qa.gating <- function(x,threshold=100) {
	# Defaults to event count threshold of 100

	print("Running QA...")
	x.class <- class(x)[1]
	if(x.class=="flowFrame") {
		counts <- length(exprs(x[,1]))
	} else if(x.class=="flowSet") {
		counts <- fsApply(x,length,use.exprs=T)
	} else {
		print("Input must be a flowSet or flowFrame")
	}

	# Find all counts less than 100 (now it's a boolean vector)
	counts.boolean <- counts<threshold
	counts.failed.position <- grep(TRUE,counts.boolean) #positions of those that failed

	# Did we get any failed counts?
	# If so, return the position
	counts.failed <- length(counts.failed.position)!=0
	if(counts.failed) {
		print("QA resulted in 1 or more warnings.")
		return(counts.failed.position)
	} else{
		print("QA succeeded")
		return(FALSE)
	}
}


### fl1transform:
### Normalizes FL1.A values in a flowset/flowframe to FSC.A values
### Should control for (at least some) non-linearity in the values
### Also makes FL1.A proportional to fluorescence/cell volume
### Used to multiply by the mean FSC.A value, but this is probably
### statistically questionable. Now multiplies by a constant (10000)
### simply to keep the values in the integer range.
###
### If you specify transform="log", it will simply do a log transform
### to FL1.A instead.

fl1transform <- function(x,transform=F) {
	# Default scaling is 10^4 and is solely for making results human-readable

	# Handle both flowFrames and flowSets
	x.class <- class(x)[1]
	if (x.class=="flowFrame") {
		#stop("This is a flowFrame, not a flowSet")
		return(transform(x,FL1.A=FL1.A/FSC.A*10^4))
	}

	# Protect the input from the modifications
	x <- x[seq(along=x)]

	# Do QA.  Reject all frames with no cells
	qa.result <- qa.gating(x,threshold=1)

	# Remove all 0-valued fluorescence results.
	# These are very likely to be artifacts and screw up some transforms
	print("Removing 0-valued fluorescence outliers")
	x <- Subset(x,rectangleGate(FL1.A=c(0.001,Inf)))

	# Transformation setup
	trans <- function(x) {
		if (transform == "fscanorm") {
			x <- transform(x,FL1.A=FL1.A/FSC.A*10^4)
		} else if (transform == "log") {
			x <- transform(x,FL1.A=log(FL1.A))
		} else if (transform == F) {
			x <- x # do nothing.  Is this necessary?
		} else {
			stop("No legitimate transform set.  Use transform=\"log\" or transform=\"fscanorm\".")
		}
	}

	if (!qa.result) {
		x <- trans(x)
#		x <- transform(x,FL1.A=FL1.A/FSC.A*10^4)
	} else {
		# For loop is inefficient, switch this out for fsApply while maintaining all cells
		for (i in qa.result) {
			x[[i]] <- trans(x)
#			x[[i]] <- transform(x[[i]],FL1.A=FL1.A/FSC.A*10^4)
		}
#		x <- fsApply(x,transform,FL1.A=FL1.A/FSC.A*10^4)
		cat(paste(
			"### Too few cells at this gating level for frame(s) \n### ",
			paste(qa.result,collapse=", "),
			".\n### These frames were not normalized.\n\n",sep=""))
	}

	return(x)
}


### flsummary:
### Get summary statistics for fluorescence, other data

flsummary <- function(flowset,channel="FL1.A",moments=F,split=F,transform=F) {
	# Number of cells (experiments) in the flowSet
	n_experiment <- length(flowset)

	# Initialize empty matrices/data frames to increase efficiency
	warnings <- c()

	if (moments == T) {
		library(moments)
	}

	# Get time of each frame in minutes of the day
	btime_raw <- fsApply(flowset,function(x)as.numeric(unlist(strsplit(keyword(x)$`$BTIM`,split=":"))))
	btime <-  apply(btime_raw,1,function(x)x[1]*60+x[2]+x[3]/60+x[4]/6000)
	time <- btime-min(btime)

	# Acquisition time - how long it took to take the sample, in seconds
	atime <- fsApply(flowset,function(x)as.numeric(keyword(x)$`#ACQUISITIONTIMEMILLI`)/1000)

	events <- fsApply(flowset,function(x)length(x[,1]),use.exprs=T)
	uL <- fsApply(flowset,function(x)as.integer(keyword(x)$`$VOL`)/1000)
	conc <- events/uL

	for (i in 1:n_experiment) {
		if (events[i] < 100) {
			warnings <- c(warnings,i)
		}
	}

	fl_mean <- fsApply(flowset,function(x)mean(x[,channel]),use.exprs=T)
	fl_median <- fsApply(flowset,function(x)median(x[,channel]),use.exprs=T)
	fl_sd <- fsApply(flowset,function(x)sd(x[,channel]),use.exprs=T)
	fl <- data.frame(fl_mean,fl_median,fl_sd)
	colnames(fl) <- paste(channel,c("mean","median","sd"),sep="")

	# Do we want mean fl values for data split into 4 evenly sized chunks?
	if (split==T) {
		split_table <- fsApply(flowset,splitFrame)
		split_table <- data.frame(matrix(unlist(split_table),ncol=4,byrow=T))
		colnames(split_table) <- paste("split",1:4,sep="")
		fl <- cbind(fl,split_table)
	}

	# Do we want the first few moments?
	if (moments == T) {
		require(moments)
		fl_var <- data.frame(fsApply(flowset,function(x)var(x[,channel]),use.exprs=T))
		fl_skew <- data.frame(fsApply(flowset,function(x)skewness(x[,channel]),use.exprs=T))
		fl_kurt <- data.frame(fsApply(flowset,function(x)kurtosis(x[,channel]),use.exprs=T))
		fl_moments <- data.frame(fl_var,fl_skew,fl_kurt)
		colnames(fl_moments) <- paste(channel,c("var","skew","kurt"),sep="")
		fl <- cbind(fl,fl_moments)
	}

	file <- fsApply(flowset,function(x)strsplit(keyword(x)$GUID,".fcs")[[1]])
	colnames(file) <- "file"

	if (length(warnings) != 0) {
		warnings <- paste(warnings,collapse=", ")
		print(paste("Warning: frame(s)",warnings,"had less than 100 events in this gate."))
	}

	# Insert empty strain and colony columns
	strain=matrix(nrow=n_experiment)
	treatment=matrix(nrow=n_experiment)

	# Put it all together
	flsummary <- cbind(time,btime,atime,events,conc,fl,file,strain,treatment)

	# Make rows filename keys
	rownames(flsummary) <- file

	# Rename the 'mean', 'median', and 'sd' columns to reflect transformations done or channel used.
	# 'FL1.A' = no transformation, 'FL1_FSC' = "fsacanorm", 'log' = "log"
	flsummary <- renameflcols(flsummary,channel=channel,transform=transform)

	return(flsummary)
}


# renaming function for flsummary data, keeping it separate for ease of use
# probably slow due to for loop
renameflcols <- function(x,channel="FL1.A",transform=F) {
	cols <- c("mean","median","sd")
	if (transform!=F) {
		if (transform=="fscanorm") {
			tname <- "FL1_FSC"
		} else if (transform=="log") {
			tname <- "log"
		} else {
			stop("invalid transform")
		}
	} else {
		return(x)
	}
	for (i in cols) {
		colnames(x)[grep(i,paste(channel,colnames(x),sep=""))] <- paste(tname,i,sep="")
	}
	return(x)
}

### summary.cyt:
### Gates a sample to all yeast, then singlet, then doublets
### Does the work of singletsummary.cyt,doubletsummary.cyt,yeastsummary.cyt
### Also calculates singlet to doublet ratio
### Returns a list of data frames, e.g. output$singlets, output$doublets, etc.
summary.cyt <- function(
		flowset,
		transform=F,
		channel="FL1.A",
		ploidy=F,
		moments=F,
		split=F,
		only=F) {

	# Number of experiments
	n_experiments <- length(flowset)

	# If using channel="FSC.A", don't use fscanorm
	if (channel=="FSC.A"&transform=="fscanorm") {
		print("Channel FSC.A selected with no transform= setting set.")
		print("Defaulting to no transform (set transform=\"log\" for log transform)")
		transform=F
	}

	# Transform FL1.A
	if (transform != F) {
		print(paste("Transforming FL1.A using",
			     transform,
		    	"transform..."
		   	)
		)
		flowset <- fl1transform(flowset,transform=transform)

	}


	# Gate the samples
	yeast <- Subset(flowset,yeastGate)
    # Remove flowFrames that have no yeast - bug in bioconductor <= 2.13
    # TODO: replace with an empty / NA row in the final dataframe
    empties <- c()
    for (i in 1:length(yeast)) {
        if (length(exprs(yeast[[i]])) == 0) {
            empties <- c(empties, i)
        }
    }
    nonempties <- 1:length(yeast)
    yeast_nonempty <- yeast[seq(along=yeast)]
    if (length(empties) != 0) {
        print(paste("The following indices had no events in the yeast gate and have been excluded from singlets/doublets:",
              paste(empties, collapse=", ")))
        for (i in empties) {
            nonempties <- nonempties[which(nonempties != i)]
        }
        yeast_nonempty <- yeast_nonempty[nonempties]
    }
	if (ploidy=="haploid") {
		print("Gating with haploid gates...")
		singlets <- Subset(yeast_nonempty,hapsingletGate)
		doublets <- Subset(yeast_nonempty,hapdoubletGate)
	} else if (ploidy=="diploid") {
		print("Gating with diploid gates...")
		singlets <- Subset(yeast_nonempty,dipsingletGate)
		doublets <- Subset(yeast_nonempty,dipdoubletGate)
	} else {
		stop('Error: You must define ploidy="haploid" or ploidy="diploid"')
	}

	if (only==F) {
	    # Normalize and summarize each subset
	    print("Summarizing all yeast events...")
	    yeastsum <- flsummary(yeast,channel=channel,moments=moments,split=split,transform=transform)

	    print("Summarizing doublets events...")
	    doubletsum <- flsummary(doublets,channel=channel,moments=moments,split=split,transform=transform)

	    print("Summarizing singlets events...")
	    singletsum <- flsummary(singlets,channel=channel,moments=moments,split=split,transform=transform)
    } else {
		if (only=="singlets") {
    	    print("Summarizing singlets events...")
    	    singletsum <- flsummary(singlets,channel=channel,moments=moments,split=split,transform=transform)
			return(singletsum)
		} else if (only=="doublets") {
    	    print("Summarizing doublets events...")
    	    doubletsum <- flsummary(doublets,channel=channel,moments=moments,split=split,transform=transform)
			return(doubletsum)
		} else if (only=="yeast") {
    	    print("Summarizing all yeast events...")
    	    yeastsum <- flsummary(yeast,channel=channel,moments=moments,split=split,transform=transform)
			return(yeastsum)
		} else {
			print("'only' must be 'singlets','doublets', or 'yeast'")
			stop()
		}
	}

	summary_list <- list(yeast=yeastsum,singlets=singletsum,doublets=doubletsum)
	return(summary_list)
}

### yeastIntSplit:
### Splits a flowSet or flowFrame by FSC.A
### Splits into n equally-spaced intervals
### n defaults to 3 intervals
### Returns a LIST

yeastIntSplit <- function(flowset,nsplit=3){
	maxval <- 800000

	# Normalize FL1.A by FSC.A (temporary)
	flowset <- transform(flowset,`FL1.A` = `FL1.A`/`FSC.A`)

	returnedvalues <- matrix(nrow=nsplit,ncol=6)

	for (i in 1:nsplit) {
		localmin <- floor(1+(i-1)*maxval/nsplit)
		localmax <- ceiling(i*maxval/nsplit)
		tempgate  <- rectangleGate("FSC.A"=c(localmin,localmax))
		tempSet <- Subset(flowset,tempgate)
		currentFL1.As <- exprs(tempSet[,"FL1.A"])

		returnedvalues[i,1] <- i
		if (mean(currentFL1.As) != "NaN") {
			returnedvalues[i,2] <- mean(currentFL1.As)
		}
		returnedvalues[i,3] <- sd(currentFL1.As)
		returnedvalues[i,4] <- localmin
		returnedvalues[i,5] <- localmax
		returnedvalues[i,6] <- length(currentFL1.As)
	}

	finaltable <- data.frame(section=returnedvalues[,1],mean=returnedvalues[,2],sd=returnedvalues[,3],
							 min=returnedvalues[,4],max=returnedvalues[,5],n=returnedvalues[,6])

	return (finaltable)
}

### yeastSampSplit
### Splits a flowFrame by its FSC.A values
### Organizes into n equally-sized populations
### n defaults to 3
### Returns a LIST of tables (kinda weird)

yeastThreePopSplit <- function(flowframe){
	flowFrametable <- flowFrame2table(flowframe)
	flowFrametable <- flowFrametable[order(flowFrametable$FSC.A),]

	firstpop <- flowFrametable[1:(floor(length(flowFrametable[,1])/3)),]
	secondpop <- flowFrametable[(ceiling(length(flowFrametable[,1])/3)):(floor(length(flowFrametable[,1])*2/3)),]
	thirdpop <- flowFrametable[(ceiling(length(flowFrametable[,1])*2/3)):(length(flowFrametable[,1])),]

	poplistmeans <- c(mean(firstpop$FL1.A),mean(secondpop$FL1.A),mean(thirdpop$FL1.A))
	poplistsds <- c(sd(firstpop$FL1.A),sd(secondpop$FL1.A),sd(thirdpop$FL1.A))
	poptable <- data.frame(size=c("Low","Mid","High"),mean=poplistmeans,sd=poplistsds)
#	poplist <- c(firstpop,secondpop,thirdpop)
	return(poptable)
}


### flowFrame2Table:
### Generates a full table (data frame) of a flowFrame's data values
### with appropriate labels
flowFrame2Table <- function(flowframe) {
	flowframetable <- data.frame(exprs(flowframe),ncol=10)
	colnames(flowframetable) <- colnames(flowframe)
	return(flowframetable)
}

# Thalf scripts
# Most of these are junk and only work with properly-formatted data frames
# which contain values of 'time', 'mean', 'strain', and 'treatment'.
# Also assumes that data frame is sorted by time on some level

# For some reason this works even when the treatment column is full of NA values
thalfall <- function(x,minval=F) {
	#HACK
	strain_treatment <- paste(x$strain,x$treatment,sep=",,")
	all_levels <- levels(as.factor(strain_treatment))
	x$strain_treatment <- strain_treatment

	# NOW I KNOW THE DATA FRAME SIZE
	tablelen <- length(all_levels)
	thalftable <- data.frame(matrix(ncol=5,nrow=tablelen))
	colnames(thalftable) <- c("thalf","min","max","strain","treatment")

	for (i in seq(along=all_levels)) {
		current_subset <- subset(x,strain_treatment==all_levels[i])

	        # Generate fit object
	        current_fit <- iaaregress(current_subset)
		c <- current_fit$coefficients[[2]]

		if (minval == F) {
	                current_thalf <- predict.thalf(current_fit)
	                minval_used <- max(c(min(current_fit$data$mean),c))
		} else {
	                current_thalf <- predict.thalf(current_fit,minval=minval)
	                minval_used <- minval
		}

	        maxval <- current_fit$data[1,"mean"] # First fluorescence data point
		thalftable[i,1] <- current_thalf
		thalftable[i,2] <- minval_used
		thalftable[i,3] <- maxval
		thalftable[i,4:5] <- strsplit(all_levels[i],split=",,")[[1]]
	}

	return(thalftable)
}

# Expects a data frame with columns of mean and time values
getthalf <- function(x,minval=0) {
	regression <- iaaregress(x)
	thalf <- predict.thalf(regression,minval=minval)
	return(thalf)
}

### Should make qplot.logistic into a geom, (geom_logistic?)
### This will allow grouping more easily

### qplot.logistic:
### Makes plotting time series + fit + thalf a little easier
### Assumes that you have 'time' and 'mean' columns

qplot.logistic <- function(timeseriesdata,minval=F) {

	# Generate fit object
	fitobject <- iaaregress(timeseriesdata)

	# Calculate thalf
	if (minval==F) {
		thalf <- predict.thalf(fitobject)
	} else {
		thalf <- predict.thalf(fitobject,minval=minval)
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
#	f <- coefficients[[5]]

	maxval <- fitobject$data[1,"mean"] # First fluorescence data point
#	maxval <- c+(d-c)/(1+exp(b*-e)) # Only works when data is very logistic-ish
	if (minval == 0) {
		minval <- max(c(min(fitobject$data$mean),c))
	}
	halfmax <- mean(c(maxval,minval))

	thalf <- exp(log( ( d - c )/( halfmax - c ) - 1 )/b + log(e)) # log model
	thalf <- log( ( d - c )/( halfmax - c ) - 1 )/b + e # non-log model
#	thalf.5 <- e + (1/b) * (log ( 1 - (d - c) / (y - c) ))^1/f

	return(thalf)
}

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
            regress0 <- drm(data=table, mean~time, fct = L.5())
#            regress0 <- drm(table, mean~time, fct = L.5(fixed=c(NA,min(table[,"mean"]),max(table[,"mean"]),NA,NA)))
        }
    }

    return(regress0)
}


### Cytometer data time series normalization
### Assumes:
### 1) a data column with 'strain type' called "strain"
### 2) a data column with 'treatment type' called "treatment"
### 3) a data column with the values to normalize called "mean"
### 4) a data column for relative time passed called, "time"
### 5) a "strain" type called "W303", wherefrom the mean will be subtracted
### 6) That you want to normalize to the very first time point's "mean" value

cyt.normalize <- function(timeseries) {
	# Figure out the metadata
	strains <- levels(factor(timeseries[,"strain"]))
	print(paste("Strains:",paste(strains,collapse=", "),sep=" "))
	treatments <- levels(factor(timeseries[,"treatment"]))
	if (typeof(timeseries[,"treatment"])=="double") {
		treatments <- as.double(treatments)
	}
	print(paste("Treatments:",paste(treatments,collapse=", "),sep=" "))

	# Subtract off mean of all W303 'mean' values
	w303 <- mean(subset(timeseries,strain=="W303")$mean)
	timeseries$mean <- timeseries$mean-w303

	# Remove W303 values because the differences will be hugely exaggerated
	# if normalized to the first point
	timeseries <- subset(timeseries,strain!="W303")

	# Divide by first time value for each strain+treatment subset
	rownames(timeseries) <- 1:length(timeseries[,1])

	for (i in strains) {
		for (j in treatments) {
			current_subset <- subset(timeseries,strain==i&treatment==j)
			current_rows <- rownames(current_subset)
			first_row <- timeseries[current_rows[1],]
			suppressWarnings(timeseries[current_rows,"mean"] <- timeseries[current_rows,"mean"]/first_row$mean)
		}
	}

	return(timeseries)
}

# flowSplit::
# Take a flowSet, split into N evenly-sized pieces M (not random, but from start to finish)
# In all M, calculate mean fluorescence for all N pieces
# Return data frame to summarize this

# n= is currently useless and all this script does is split into 4

splitSet <- function(flowset, n=4) {
        # number of wells
        x <- length(group2)
        x.ind <- seq(along=group2)

        fltable <- data.frame(matrix(nrow=x,ncol=n))
        colnames(fltable) <- paste("mean",1:n,sep="")

        for (i in x.ind) {
                # choose a flowframe
                flowset[[i]]
                # get the raw fluorescence data
                fl_raw <- exprs(flowset[[i]][,3])
                fl_length <- length(fl_raw)
                # piece size
                piece_size <- floor(fl_length/4)

                fltable[i,] <- c(
                        mean(fl_raw[1:piece_size])
                        ,mean(fl_raw[(piece_size+1):(2*piece_size)])
                        ,mean(fl_raw[(2*piece_size+1):(3*piece_size)])
                        ,mean(fl_raw[(3*piece_size+1):fl_length])
                        )
        }

        return(fltable)
}

splitFrame <- function(flowframe,n=4) {
                # get the raw fluorescence data
                fl_raw <- exprs(flowframe[,3])
                fl_length <- length(fl_raw)
                # piece size
                piece_size <- floor(fl_length/n)


                fl_vec <- c(
                        mean1=mean(fl_raw[1:piece_size])
                        ,mean2=mean(fl_raw[(piece_size+1):(2*piece_size)])
                        ,mean3=mean(fl_raw[(2*piece_size+1):(3*piece_size)])
                        ,mean4=mean(fl_raw[(3*piece_size+1):fl_length])
                        )
#	        names(fl_vec) <- paste("mean",1:n,sep="")
		fltable <- data.frame(t(fl_vec))
		return (fltable)
}

# Takes in a flowSet and checks for empty flowFrames (flowCore's Subset fails if there's no events in a flowFrame)
# Also returns a flowSet, with such frames removed if necessary. Posts a note about those which were removed.

qa.flowSet <- function(flowset_in) {
    ### TODO:
    ### Make sure all columns have . substituted for -


	# There is very likely a better way to do this than a for loop
	# Get the flowFrame positions that have zero events in them, for excluding
	pos <- c()
	for (i in 1:length(flowset_in)) {
		if (length(exprs(flowset_in[[i]][,1]))==0) {
			pos <- c(pos,i)
		}
	}
	if (length(pos) > 0) {
		# remove empty frames.
		flowset_out <- flowset_in[(1:length(flowset_in))[-pos]]
		print(paste("The following frames had no events and were removed: ",paste(pos,",",sep=""),".",sep=""))
		print(paste("The flowSet is now ",length(flowset_out)," frames long."))
		return(flowset_out)
	} else {
		# unchanged flowSet
		return(flowset_in)
	}
}

#outputs data frame formatted with several parameters used for modeling.
#assumes diploids and only returns singlet data
#requires flowset and strain vector

# on 2012-1-3, replaced 'split' FL1/FSC values by background-subtracted, changing the order of columns. All 'modelformat' data needs to be reprocessed.
modelingFormat <- function(flowset,strain_vec,baseline="noYFP",normalize=F,ploidy="diploid") {
	# Make sure strain vector includes correct baseline value
	if ( sum(as.numeric(strain_vec==baseline)) == 0 ) {
		stop("No baseline strain found in strain_vec (default is noYFP)")
	}

	# Generate data frames from which to take data
    # Raw FL1.A table
	raw <- summary.cyt(flowset,transform=F,only="singlets",split=T,ploidy=ploidy)
    # Raw FSC.A table
#	fsc <- summary.cyt(flowset,channel="FSC.A",only="singlets")
    # Normalized data table
#	fl1_fsc <- summary.cyt(flowset,transform="fscanorm",only="singlets",split=T)

	out <- raw

	out$grp <- NA
	out$strain <- strain_vec
	out$afb <- NA
	out$rep <- NA


	FL_idx <- which(colnames(out)=="FL1.Amean")
	FL_idx <- c(FL_idx,which(colnames(out)=="FL1.Amedian"))
	FL_idx <- c(FL_idx,which(colnames(out)=="FL1.Asd"))

	colnames(out)[FL_idx] <- c("FL1.A","median","sd")
#	out <- cbind(out,FSC.A=fsc$FSC.Amean)
#	out <- cbind(out,FL1_FSC=fl1_fsc$FL1_FSCmean)

#	out$FL1_FSC_norm <- out$FL1_FSC-mean(subset(out,strain==baseline)$FL1_FSC)

	# Changed 2012-1-3
	colnames(out)[which(colnames(out)=="split1")] <- "FL1.A_bs_1"
	colnames(out)[which(colnames(out)=="split2")] <- "FL1.A_bs_2"
	colnames(out)[which(colnames(out)=="split3")] <- "FL1.A_bs_3"
	colnames(out)[which(colnames(out)=="split4")] <- "FL1.A_bs_4"

#	out <- cbind(out,FL1.A_1=fl1_fsc[,"split1"])
#	out <- cbind(out,FL1.A_2=fl1_fsc[,"split2"])
#	out <- cbind(out,FL1.A_3=fl1_fsc[,"split3"])
#	out <- cbind(out,FL1.A_4=fl1_fsc[,"split4"])

#	out <- cbind(out,FL1.A_1=fl1_fsc[,"split1"])
#	out <- cbind(out,FL1.A_2=fl1_fsc[,"split2"])
#	out <- cbind(out,FL1.A_3=fl1_fsc[,"split3"])
#	out <- cbind(out,FL1.A_4=fl1_fsc[,"split4"])


    # background-subtract all FL1.A values
	out$FL1.A_bs <- out$FL1.A-mean(subset(out,strain==baseline)$FL1.A)

	out$FL1.A_bs_1 <- out$FL1.A-mean(subset(out,strain==baseline)$FL1.A_bs_1)
	out$FL1.A_bs_2 <- out$FL1.A-mean(subset(out,strain==baseline)$FL1.A_bs_2)
	out$FL1.A_bs_3 <- out$FL1.A-mean(subset(out,strain==baseline)$FL1.A_bs_3)
	out$FL1.A_bs_4 <- out$FL1.A-mean(subset(out,strain==baseline)$FL1.A_bs_4)

	if(normalize==T) {
		ddply(out,c("strain","treatment"),transform,norm1=FL1_FSC_norm/min(FL1_FSC_norm[1:3]))
	}

	return(out)
}

# Produces a normalized fluorescence column 'normed'
# Expects the 'FL1.A_bs' column to exist (not hard to extend to others/make it user selectable)
# Has two different methods, version 1 and version 2, described in the script
addnorm <- function(frame,factor_in=c("strain","treatment"),method=1,column="FL1.Amean_bs") {
    library(plyr)
	if ( (sum(colnames(frame)==column)) == 0 ) {
		if( (sum(colnames(frame)=="FL1.A_bs")) == 0 ) {
			stop("Could not find the background-subtracted values column. \
				  This script requires that there be a column named \
				  FL1.Amean_bs, FL1.A_bs, or the user-defined column using\
				  column='desired-column'")
		} else {
			column <- "FL1.A_bs"
		}
	}

	if (method==1) {
		# Default normalization method. Takes highest point in dataset grouped by 'factor_in' and sets it to 1,
		# divides all other values by that number. This method is default because it works regardless of
		# whether the data is a time series.
		estimate_0 <- function(x) {
			x[,"normed"]=x[,column]/max(x[,column])
			return(x)
		}
	} else if (method == 2) {
		# Version 2 - takes the mean value of all time points which are less than 0, after grouped by 'factor_in'.
		# Sets this to the value by which all other data points in that group are divided
		# Therefore, no value is actually '1' except by very rare chance
		# Requires a time series with negative time values to work
		estimate_0 <- function(x) {
			normresult <- x[,column]/mean(x[x$time<0,column])
			x <- cbind(x,normed=normresult)
			return(x)
		}
	} else if (method == 3) {
		# Version 3 makes a fit line to all pre-zero time points and infers the y-intercept
		# Requires a time seriesw ith negative time values to work
		estimate_0 <- function(x) {
			prezero_points <- x[x$time<0,]
			prezero_fit <- lm(prezero_points[,column]~prezero_points[,"time"])
			prezero_intercept <- prezero_fit$coefficients[1] # intercept
			normresult <- x[,column]/prezero_intercept
			x <- cbind(x,normed=normresult)
			return(x)
		}
	} else {
		stop("You must define method=1, method=2, or method=3)")
	}

	# Check for negative time values
	if (sum(frame$time<0)==0) {
		if (method==2|method==3) {
			stop("To use methods 2 or 3, the input data frame must have negative time values for each normalized data subset")
		}
	}

	# Run the chosen estimation function and apply it
	frame <- ddply(frame,factor_in,estimate_0)
	return(frame)
}

addbs <- function(frame,column="FL1.Amean",baseline="noYFP") {
        frame[,paste(column,"_bs",sep="")] <- frame[,column]-mean(subset(frame,strain==baseline)[,column])
	return(frame)
}

# Generate a data frame that's useful for plotting overlapping density plots with ggplot.
# At the moment it's very finnicky. It expects a huge data frame as input that's made with the exprs command on a bunch
# bunch of flowFrames, with an extra column called 'exptime' signifying the time of each well's acquisition.
# It shouldn't be hard to make it accept a flowSet and do this automatically.

density_frame <- function(frame,param="FL1.A") {
	# generate the general plot
	frame_dens <- tapply(frame[,param],frame$exptime,function(x)data.frame(x=density(x)$x,y=density(x)$y))

	# coerce the list into a data frame
	frame_temp <- frame_dens[[1]]
	for (i in 2:length(frame_dens)) {
		frame_temp <- rbind(frame_temp,frame_dens[[i]])
	}

	frame_dens <- frame_temp

	# apply the exptime label to it. Very finnicky here as well, NOT at all happy with arbitrary input
	frame_dens$exptime <- rep(levels(factor(frame$exptime)),each=512)

	# normalize it within the groups
	frame_dens <- ddply(frame_dens,"exptime",transform,y_norm=y/max(y))

	return(frame_dens)
}

# Script to reprocess cytometer data following a given pattern. Uses modelFormat for reprocessing.
# Expects as input a directory with one folder 'csvs' full of prior csvs and the matching source data in 'source'
# Outputs to same directory in 'newcsvs' folder, overwriting anything that exists there.
reprocess <- function(directory_in,nick_type=T) {
    # Delete 'newcsvs' if it already exists
    unlink(paste(directory_in,"/newcsvs",sep=""),recursive=T)

    # Get csv list and create fresh 'newcsv' dir
	csvlist <- list.files(paste(directory_in,"/csvs",sep=""))
    dir.create(paste(directory_in,"/newcsvs",sep=""))

    # massive loop to reprocess each piece of data
    for (i in csvlist) {
        csv <- read.csv(paste(directory_in,"/csvs/",i,sep=""))
        csv_files <- paste(csv$file,".fcs",sep="")
        expname <- gsub("_","",substr(i,1,3))

        # find dir for this experiment
        dirs <- list.files(paste(directory_in,"source/",sep=""))
        dir <- grep(expname,dirs)
        message("Reprocessing ",expname,"...")

        # Read in flowSet, trim to all that match csv
        fs <- read.flowSet(path=paste(directory_in,"source/",dirs[dir],sep=""),alter.names=T)

        fs_exp <- sampleNames(fs)
        fs_matches <- which(fs_exp %in% csv_files)

        fs_trimmed <- fs[fs_matches]

        ####################################
        # Currently uses 'modelingFormat' to reprocess data. Can substitute anything in here
        ####################################

        newcsv <- modelingFormat(fs_trimmed,csv$strain)

        # Re-add the metadata: treatment, grp, afb, rep
        for (cols in c("treatment","grp","afb","rep")) {
            newcsv[,cols] <- csv[,cols]
        }

        write.csv(newcsv,file=paste(directory_in,"newcsvs/",i,sep=""))
    }
    message("Finished. Files are in:")
    message(directory_in,"newcsvs/")
}

flowFrame.gettime <- function(flowframe) {
    time_raw <- as.numeric(unlist(strsplit(keyword(flowframe)$`$BTIM`,split=":")))
    time <- time_raw[1]*60+time_raw[2]+time_raw[3]/60+time_raw[4]/6000
    return(time)
}

YeastCytSummary <- function(inpath,ploidy=F,only="singlets",channel="FL1.A") {
    fs <- read.flowSet(path=inpath,alter.names=T)
    if (ploidy=="diploid"|ploidy=="haploid") {
        fs_sum <- summary.cyt(fs,only=only,ploidy=ploidy,channel=channel)
    } else {
        stop('Must define ploidy= as "haploid" or "diploid"')
    }
    filename=paste( "summary-",format(Sys.time(), "%Y-%m-%d--%H-%M-%S"),".csv",sep="")
#    print(paste("~/Desktop",filename,sep=""))
    write.csv(fs_sum,paste("~/Desktop/",filename,sep=""))
    message("File was written to Desktop/",filename)
}


