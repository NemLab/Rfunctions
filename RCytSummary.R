# R scripts by Miles Gander for E. coli cytometry data processing

# Before using:
# Script requires the flowCore toolset for R. To install run these commands within R:
#    source("http://www.bioconductor.org/biocLite.R")
#    biocLite("flowCore")
#
# Usage notes:
# After loading a fresh R session, you need to load this file:
# source("~/Desktop/RCytSummary.R")
# Substitute ~/Desktop/RCytSummary.R for wherever you actually saved the file.
#
# The command to run is called CytSummmary, and it takes two arguments: 
#   1) the input path for your cytometry data, which should contain *only* a bunch of FCS files.
#   2) the output path for the summary data, in csv format (can open with excel or matlab or R or python)
# Example: CytSummary("~/Desktop/CytometryFiles/experiment-2012-02-17/","~/Desktop/")
# This would read all files in Desktop/CytometryFiles/experiment-2012-02-17, gate to the E. coli population, 
# then provide summary information for each well in a csv on the desktop, automatically named by the date and time.
require(flowCore)

CytSummary <- function(filepath, outpath) {
    require(flowCore)
    library(flowCore)

	## Read in flowSet

	fs = read.flowSet(path=filepath, alter.names=T)

	fs_gated=Subset(fs, ecoligate)  ##gating all data to Rob's ecoli gate
	fs_qad <- fs_gated
#    fs_qad <- qa.flowSet(fs_gated)
	gfp=flsummary(fs_qad, channel="FL1.A") ## gives summary of GFP data from each well

	rfp=flsummary(fs_qad, channel="FL3.A") ## gives summary of RFP data from each well

	FL3_index=which(colnames(rfp)=="FL3.Amean")
	FL1_index=which(colnames(gfp)=="FL1.Asd")

	totalsummary=cbind(gfp[,1:FL1_index],rfp[,FL3_index:length(colnames(rfp))])

	filename=paste( "summary-",format(Sys.time(), "%Y-%m-%d--%H-%M-%S"),sep="")	

	## Following allows for it not to matter if the outpath has a / or not at the end of it 
	## As well as adding .CSV to the end of the file path to make it easier to open in Excel

	outpathlength=length(outpath)
	outpathvector=strsplit(split="",outpath)

	last=outpathvector[outpathlength]

	if(last=="/"){
		outpath=outpath
	} else{

		outpath=paste(outpath,"/",sep="")
		}

  	write.csv(totalsummary,file=paste(outpath,filename,".csv",sep=""))


	return(totalsummary)
}

ecoligate<<- rectangleGate(filterId="ecoligate",.gate=list(SSC.H=c(1e2,1e4), FSC.H=c(1e3,1e5)))

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
        tname <- channel
    }
    for (i in cols) {
        colnames(x)[grep(i,colnames(x))] <- paste(tname,i,sep="")
    }
    return(x)
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
