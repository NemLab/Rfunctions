library(Biostrings)


### TODO:
### Write up visualization of mismatches/overlaps/coverage
### Ideal for such things would be the ability to place all of the sequences in reference to the subject

### Gateway primer overhangs, set up to enable N-terminal fusions.  Also adds 3 Adenine residues
### to resemble ideal yeast Kozak sequence.
gateway1 <- "AAAAAGCAGGCTTCAAA"
gateway2 <- "AGAAAGCTGGGTG"


### Alignment Stuff

# Biostrings apparently doesn't let you coerce a bunch of fastas in a given dir
# into a single DNAStringSet.  This function does that.
# Right now it is very sensitive to user input/formatting.  Make it more robust?
# Expects a dir containing all the .fasta and .seq (both in fasta format) files you want

read.DNAStringSet.path <- function(path="./") {
	library('Biostrings')

	# Grab all *.fasta file paths in the dir
	filelist <- paste(path,
		grep("*.fasta$|*.seq$",list.files(path),value=T),sep="")
	x <- read.DNAStringSet(filelist)
	return(x)
}

# Reads clustal alignment file, realigns as pairwise alignment using 'multipletopairwise', then reports alignanalyze
readAlign <- function(alignmentfile) {
    library('Biostrings')
	a <- read.DNAMultipleAlignment(alignmentfile,format="clustal")
	a <- multipletopairwise(a)
	alignanalyze(a)
	return(a)
}

alignread <- function(path=NULL,a) {
	library('Biostrings')
    # This script absolutely requires a sequence named 'reference' in the FASTA format.

	if (!is.null(path)) {
		# If path doesn't include trailing slash, add it
		path_len <- nchar(path)
		if (substr(path,path_len,path_len)!="/") {
			path <- paste(path,"/",sep="")
		}

		x <- trimSequencing(read.DNAStringSet.path(path))
		references.index <- names(x)=="reference"
	

	    # Prompt for input reference if one doesn't exist
	    if ( sum(references.index) == 0 ) {
	        library(tcltk)
	        getfile <- function() {
	            fileName<-tclvalue(tkgetOpenFile(title="Select fasta format reference sequence file"))
	            if (!nchar(fileName))
	                tkmessageBox(message="No file was selected!")
	            else
	                return(fileName)
	        }
	        ref_in <- getfile()
			ref_in_bs <- read.DNAStringSet(ref_in)
			names(ref_in_bs) <- "reference"
	        x <- append(x,ref_in_bs)
		}
	} else {
		# In the case of no input files supplied to the script, prompt
		# for the reference sequence first, then for the sequencing files
		library(tcltk)
		get_reference <- function() {
			fileName <- tclvalue(tkgetOpenFile(title="Select fasta format reference sequence file"))
			if (!nchar(fileName))
				tkmessageBox(message="No file was selected!")
			else
				return(fileName)
		}
		ref_in <- get_reference()
		x <- read.DNAStringSet(ref_in)
		names(x) <- "reference"
		get_seqs <- function() {
			fileNames <- tk_choose.files(caption="Select fasta format sequencing results (.seq, .fasta)")
			if (!nchar(fileNames))
				tkmessageBox(message="No file was selected!")
			else
				return(fileNames)
		}
		seqs_in <- get_seqs()
		x <- append(x,trimSequencing(read.DNAStringSet(seqs_in)))
	}

	# QA: make sure there's more than one sequence
	if (length(x)==1) {
		stop("Only one sequence input found. Need at least two for an alignment!")
	}

    # align and analyze sequences
    y <- alignseqs(x)
    alignanalyze(y)

	final_result <- list(alignment=y,stringset=x)
	return(final_result)
}

alignseqs <- function(x) {
    # Make an alignment (first pass)
    cat("Aligning... ")
	y <- pairwiseAlignment(
		pattern=x[names(x)!="reference"],
		subject=x[names(x)=="reference"],
		type="overlap"
	)
    # Check for bad scores. They may indicate that sequences need to be reverse complemented
    seq_scores <- score(y)
    if ( any (seq_scores < 10) ) {
        for (i in which(seq_scores < 10) ) {
            x[names(x)!="reference"][i] <- reverseComplement(x[names(x)!="reference"][i])
        }
    	y <- pairwiseAlignment(
    		pattern=x[names(x)!="reference"],
    		subject=x[names(x)=="reference"],
    		type="overlap"
    	)
    }

    seq_scores <- score(y)
    # If there's still bad scores, report them

	bad_alignment_warning <- function() {
        message("")
        message("WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING ")
        message("WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING ")
        message(" ")
        message("Warning! The following sequences had too many mismatches and did not align:")
        message(paste(names(x[names(x)!="reference"])[which(seq_scores < 30)],collapse=", "))
        message(" ")
        message("WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING ")
        message("WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING ")
        message("")
	}
    if ( any (seq_scores < 30) ) {
		bad_alignment_warning()
	}

    message("Done")
    return(y)

}

### alignanalyze:
### Takes a pairwise alignment and produces a report of mismatches, insertions, deletions,
### and sequence coverage.
### Currently broken - it is not correctly reporting the locations of
### deletions in Sequencing/pGP8G-mCherry-IAA17/20120320-alignment/

alignanalyze <- function(pairwise_alignment) {
	# make the alignment variable shorter so it's easier to work with
	y <- pairwise_alignment

	# Make a mismatch+deletion+insertion table
	mid <- gen_mid_table(y)

	# Start reporting any mismatches and sequence coverage
	message("")
	message("##############################################")
	message("#     Mismatch/Deletion/Insertion Report     #")
	message("##############################################")

	if (length(mid)==0) {
		message("No mismatches, deletions, or insertions found.")
		message()
	}

	# Generate string vectors for each result that has a discrepancy,
	# then generate a report
	mismatch_display_list <- list()
	for (i in unique(mid$number)) {
		mismatch_name <- unique(mid[mid$number==i,"name"])
		mismatch_locations <- mid$position[mid$number==i]
		mismatch_display_list[[mismatch_name]] <- gen_aligned_segment(y,mismatch_name,mismatch_locations)
		mismatch_message <- paste("     Disrepancies found in sequencing result",mismatch_name,":")
		message("     -----------------------------------------")
		message(mismatch_message)
		message("     -----------------------------------------")
		message()

		for (j in mismatch_locations) {
			midtype <- as.character(mid[mid$position==j,"type"])
			midtype <- paste(toupper(substring(midtype,1,1)),tolower(substring(midtype,2)),sep="",collapse="")
			message("     ",midtype," at position: ",j)
			report_aligned_segment(mismatch_display_list[[mismatch_name]],j,10)
			message()
		}

	}


	# Can find coverage of sequencing w/ 'coverage' method (coverage(y))
	# 	Can't distinguish (yet) between missing sequencing results at the end
	# 	and sequencing results that extend past the pattern (reference) sequence
	# 	But maybe this only applies to multiple alignments coerced into pairwise alignments?
	message("##############################################")
	message("#              Coverage Report               #")
	message("##############################################")

	y_coverage <- coverage(y)
	y_coverage_length <- runLength(y_coverage)
	y_coverage_val <- runValue(y_coverage)

	if ( all(y_coverage_val != 0) ) {
		message("     No gaps in sequencing coverage.")
	} else {
		# Find the gaps
		y_gaps <- grep(0,y_coverage_val)

		# Inefficiently generate the sites
		y_coverage_stops <- c()
		for (i in 1:length(y_coverage_length)) {
			y_coverage_stops <- c(y_coverage_stops,
					       sum(y_coverage_length[1:i])
						      )
		}

		y_coverage_starts <- y_coverage_stops[y_gaps]-y_coverage_length[y_gaps]+1
		y_coverage_result <- paste("          ",apply(cbind(y_coverage_starts,y_coverage_stops[y_gaps]),1,paste,collapse="-"),collapse=", \n")


		message()
		message("     -----------------------------------------")
        message("     Gaps in sequencing at reference locations: ")
		message("     -----------------------------------------")
		message(y_coverage_result)
	}
	message("")
}

### gen_mid_table:
### Generates a table of mismatches, insertions, and deletions
### for a given pairwise alignment (Biostrings object)

gen_mid_table <- function(pairwise_alignment) {
	y <- pairwise_alignment
	
	#################
	# Find the locations of all mismatches, insertions, and deletions
	#################

	# Get the 'start' and 'width' info from each pairwise alignment
	# the 'range' accessor doesn't seem to work, so used direct attributes
	y_ranges_raw <- attributes(subject(y))$range
	y_ranges <- data.frame(start=start(y_ranges_raw),end=end(y_ranges_raw),width=width(y_ranges_raw))
	# Get the name of each sequencing result
	y_ranges$name <- names(unaligned(pattern(y)))
	# Calculate disjoint binning to minimize number of levels needed
	y_ranges$bin <- disjointBins(IRanges(y_ranges$start,y_ranges$end))

	# Add in full coverage as a single data type thing
	y_coverage <- y_ranges
	y_coverage$name <- "coverage"
	y_coverage$bin <- max(y_ranges$bin)+1


	# mismatch insertion deletion (MID) table. Remember to add insertions and deletions...

	mismatches <- c()
	insertions <- c()
	deletions <- c()

	# mismatches
	if (length(row.names((mismatchTable(y))))!=0) {
		mismatches <- data.frame(position=mismatchTable(y)$SubjectStart,number=mismatchTable(y)$PatternId,type="mismatch")
		mismatches$name <- y_ranges$name[mismatches$number]
	}
	# insertions
	if (length(as.data.frame(insertion(y))[,1])) {
		insertions <- as.data.frame(insertion(y))
		insertions$name <- y_ranges$name[insertions$space]
		colnames(insertions) <- c("number","position","stop","width","name")
		insertions$type <- "insertion"
		# fix positions from relative to pattern to relative to subject
		for (i in y_ranges$name) {
			insertions[insertions$name==i,"position"] <- insertions[insertions$name==i,"position"] + y_ranges[y_ranges$name==i,"start"]
		}
	}
	
	# deletions
	if (length(as.data.frame(deletion(y))[,1])) {
		deletions <- as.data.frame(deletion(y))
		deletions$name <- y_ranges$name[deletions$space]
		colnames(deletions) <- c("number","position","stop","width","name")
		deletions$type <- "deletion"
		# fix positions from relative to pattern to relative to subject
		for (i in y_ranges$name) {
			deletions[deletions$name==i,"position"] <- deletions[deletions$name==i,"position"] + y_ranges[y_ranges$name==i,"start"]
		}
	}

	# columns to keep from each result
	keep_cols <- c("number","position","type","name")
	mid <- rbind(mismatches[,keep_cols],insertions[,keep_cols],deletions[,keep_cols])
	if (length(mid)!=0) {
		# Always use same shapes for mismatches/insertions/deletions, and show key
		mid$type <- factor(mid$type,levels=c("mismatch","insertion","deletion"))
		# Apply bins from before so they plot in the right place:
		mid$bin <- NA

		for (i in levels(factor(mid$name))) {
			mid[mid$name==i,"bin"] <- y_ranges[y_ranges$name==i,"bin"]
		}		
	}
	return(mid)
}

# gen_aligned_segment:
# Makes a graphical representation of matches/mismatches for an alignment vs. the reference
gen_aligned_segment <- function(pairwise_alignment_object,seq_id,mismatch_locations) {
	# Convert alignment into gapped strings
	pa <- as.character(aligngappedseq(pairwise_alignment_object)) # convert alignment to strings
	
	pa_seq <- pa[seq_id] # Isolate just the result desired
	pa_ref <- pa['reference'] # Isolate the reference sequence

	# Generate a string of matches with pipes: |
	matches_string <- rep("|",nchar(pa_ref))
	# Mismatches are non-identities
	matches_string[mismatch_locations] <- " "
	# Gaps are non-identities
	matches_string[strsplit(pa_seq,split="")[[1]]=="-"] <- " "
	matches_string[strsplit(pa_ref,split="")[[1]]=="-"] <- " "

	matches_string <- paste(matches_string,collapse="")
	# Return a vector of 3 strings - reference, matches string, and aligned sequence
	string_vec <- c(pa_ref,matches_string,pa_seq)
	names(string_vec) <- c("reference","identities","seq")
	return(string_vec)
}

### report_aligned_segment:
### prints out a short report of an aligned string (requires a segment vector
### from gen_aligned_segment).
report_aligned_segment <- function(segment_vector,location,n_flanking) {
	seq_start <- location-n_flanking
	seq_end <- location+n_flanking

	# handle the case where the sequence is at the very front or end of an alignment
	if (seq_start < 0) {
		seq_start <- 0
	} else if (seq_end > nchar(segment_vector[1])) {
		seq_end <- nchar(segment_vector[1])
	}

	line1 <- paste("          ",paste(rep(" ",n_flanking),collapse=""),"*",sep="")
	line2 <- paste("          ",substr(segment_vector["reference"],seq_start,seq_end),sep="")
	line3 <- paste("          ",substr(segment_vector["identities"],seq_start,seq_end),sep="")
	line4 <- paste("          ",substr(segment_vector["seq"],seq_start,seq_end),sep="")
	report <- paste(line1,line2,line3,line4,sep="\n")
	message(report)
}

alignmenttostring <- function(alignment) {
    library(Biostrings)
    subject_seq <- toString(unaligned(subject(alignment[1])))
    pattern_seq <- c()
    for (i in 1:length(alignment)) {
        pattern_seq <- c(pattern_seq,toString(pattern(alignment[i])))
    }
    pattern_names <- names(unaligned(pattern(alignment)))
    subject_name <- names(unaligned(subject(alignment)))
    stringset <- DNAStringSet(c(subject_seq,pattern_seq))
    names(stringset) <- c(subject_name,pattern_names)
    return(stringset)
}

### write.pairwise:
### Allows one to write a pairwise alignment into a fasta file (finally), allowing the saving of alignments etc.
### Requires the stringset so it can have the reference sequence and the alignment to get the rest.
### Another idea: mismatchTable/indel provide a lot of info about the sequences. Could use them to reconstruct full alignment
### strings.

write.pairwise <- function(alignment,file=NULL) {
    if (is.null(file)) {
        stop('You must define a file path with file="~/path/to/file"')
    } else {
        stringset <- alignmenttostring(alignment)
        write.XStringSet(stringset,file=file)
    }
}

multipletopairwise <- function(x,subject=1) {
	# Defaults to using the first entry as the subject (reference) sequence
	# apparently I can't make it a fixedsubject version?
	# https://stat.ethz.ch/pipermail/bioc-sig-sequencing/2010-August/001400.html
	# but should instead use 'PairwiseAlignedXString"?
	x_stringset <- unmasked(x)
	x_stringset_names <- names(x_stringset)
	x_stringset_index <- seq(along=x_stringset)

	# Remove gaps before redoing alignment, as it messed them up
	x_stringset_degapped <- x_stringset
	for (i in x_stringset_index) {
		x_stringset_degapped[[i]] <- DNAString(gsub("-","",as.character(x_stringset[[i]])))
	}
	names(x_stringset_degapped) <- names(x_stringset)
	x_stringset <- x_stringset_degapped

	# Assume that first entry is the subject
	x_subject <- x_stringset[[1]]

	# Ignore anything with 'Consensus' in the name, it's an artifact of manipulating seaview
	x_stringset_pattern_index <- grep("Consensus",x_stringset_names,invert=T)
	x_pattern <- x_stringset[x_stringset_pattern_index[2:length(x_stringset_pattern_index)]]

	y <- pairwiseAlignment(
		pattern=x_pattern,
		subject=x_subject,
		type="overlap"
	)

	return(y)

}

### Finnzymes calculator (optimizing)
finnzymes <- function(s, conc=0.5*10^-7, salt=.05, method="f") {
	# Currently only calculates Finnzymes Tm, because we don't care about the others
	# Some code is therefore redundant, but could be adapted to allow other methods
	# Finnzymes defaults to modified breslauer parameters. The only modification seems to be ignoring
	# initiating with AT vs GC as well as 'symmetry correction', 
	# all which contribute to dS in Breslauer, K. J.; Frank, R.; Blocker, H., Marky, L. A. Proc. Natl. Acad. Sci. USA 1986, 83, 3746-3750.

	s <- toupper(s)
	# Pass on NA
	if ( is.na(s) ) {
		return(NA)
	}

	# Breslauer's Parameters
	BRdeltaHParams <- data.frame(AA=-9.1)
	BRdeltaHParams["TT"] = -9.1;
	BRdeltaHParams["AT"] = -8.6;
	BRdeltaHParams["TA"] = -6.0;
	BRdeltaHParams["CA"] = -5.8;
	BRdeltaHParams["TG"] = -5.8;
	BRdeltaHParams["GT"] = -6.5;
	BRdeltaHParams["AC"] = -6.5;
	BRdeltaHParams["CT"] = -7.8;
	BRdeltaHParams["AG"] = -7.8;
	BRdeltaHParams["GA"] = -5.6;
	BRdeltaHParams["TC"] = -5.6;
	BRdeltaHParams["CG"] = -11.9;
	BRdeltaHParams["GC"] = -11.1;
	BRdeltaHParams["GG"] = -11.0;
	BRdeltaHParams["CC"] = -11.0;
	BRdeltaHParams["termAT"] = 0;
	BRdeltaHParams["termGC"] = 0;

	BRdeltaSParams <- data.frame(AA=-24.0)
	BRdeltaSParams["TT"] = -24.0;
	BRdeltaSParams["AT"] = -23.9;
	BRdeltaSParams["TA"] = -16.9;
	BRdeltaSParams["CA"] = -12.9;
	BRdeltaSParams["TG"] = -12.9;
	BRdeltaSParams["GT"] = -17.3;
	BRdeltaSParams["AC"] = -17.3;
	BRdeltaSParams["CT"] = -20.8;
	BRdeltaSParams["AG"] = -20.8;
	BRdeltaSParams["GA"] = -13.5;
	BRdeltaSParams["TC"] = -13.5;
	BRdeltaSParams["CG"] = -27.8;
	BRdeltaSParams["GC"] = -26.7;
	BRdeltaSParams["GG"] = -26.6;
	BRdeltaSParams["CC"] = -26.6;
	BRdeltaSParams["termAT"] = 0;
	BRdeltaSParams["termGC"] = 0;

	deltaSParams <- BRdeltaSParams
	deltaHParams <- BRdeltaHParams

	dS <- 0
	dH <- 0
	# Calculate dS and dH based on base neighbor pairs
	s.length <- nchar(s)

	pairs <- names(BRdeltaHParams)

	counts <- rep(0,length(pairs))
	names(counts) <- pairs

	# For loop solution
#	pair_results <- c()
#	for (i in 1:(nchar(s)-1)) {
#		pair_results <- c(pair_results,substr(s,i,(i+1)))
#	}


	# sapply solution
#	pair_fun <- function(x){
#		substr(s,x,(x+1))
#	}

#	pair_results <- sapply(1:(nchar(s)-1),function(x)pair_fun(x))

#	counts_n <- table(pair_results)

#	countfun <- function(x) {
#		counts[x] <- counts_n[x]
#	}

#	counts <- sapply(names(counts_n),function(x)countfun(x))
    overcount <- function(st,p) {
        x <- gregexpr(paste(paste("(?<=",substr(p,1,1),")",sep=""),substr(p,2,2),sep=""),st,perl=T)[[1]]
        return(length(x[x>0]))
    }

    pairs <- names(BRdeltaHParams)
    counts <- sapply(pairs,function(x)overcount(s,x))


	dH <- sum(counts*BRdeltaHParams)
    dS <- sum(counts*BRdeltaSParams)
	# Corrections
	dH <- dH*10^3 #"cause specs are in kcal"
	dS <- dS - 12.4
	dH <- dH - 3400

	# Salt Calculation
	sc <- 16.6 * log(salt) / log(10.0)

	# Calculate and adjust Tm to °C
	tm <- dH / (1.9872 * log(conc/16) + dS) + sc
	tm <- tm - 273.15

	return(tm)
}


# Take a DNA input, classify it, convert to single character entry
qa.atgc <- function(sequence,output.as.char=T) {
	## I/O
	# Is the sequence object a DNAString?
	if (class(sequence)=="DNAString") {
		sequence <- as.character(sequence)
		sequence <- toupper(sequence)
	# How about a string with As, Ts, Gs, and Cs?
	} else if (class(sequence)=="character") {
		if (length(sequence) > 1) {
			sequence <- paste(sequence,collapse="")
		}
		sequence <- toupper(sequence)
		if (length(grep("[^ATGC]",sequence))) {
			stop("Sequence contains non-ATGC characters")
		}
	} else {
		stop("Input sequence must be character string or DNAString object")
	}
	return(sequence)
}

### Primer Design:
### designprimer:
### Designs primers for a given sequence or set of sequences.
### It expects the sequence to be in the form of a single string. 'qa.atgc' in the script handles other formats to an extent.
### Primers are reported as 'forward' and 'reverse' (I should add an option to design only one or the other)
### where 'forward' is in the direction of the string given, 5' -> 3'
### 	Defaults: constrains primers to less than 60 bp, no tail support yet
### 	Adds bases until a Tm within 1°C of the desired Tm has been reached.  
###	If it overshoots, it sticks with the slightly-overshot primer
###	Defaults to a desired Tm of 72 according to the finnzymes Tm calculator
### TODO: Allow changing of Tm calculation method, automated tail addition
### If it can't find a primer in your desired Tm, it picks the next highest within a threshold (tm.errorplus).

# TODO: handle newlines

designprimer <- function(sequence,tm=72,tm.errorplus=3,tm.errorminus=1,sixtymer=T,tail=F,endGC=T,atleast=F,minstart=10) {

	sequence <- qa.atgc(sequence)

	# See if it's too short to begin with
	if( finnzymes(sequence) < tm-tm.errorminus ) {
		stop("Tm of full-length sequence is lower than desired Tm and error parameters allow.")
	}

	## Iterating primers vs. Tm
	# Starts checking for good Tms at 10 bases
	primer.tm <- 0
	bases <- minstart
	sequence <- substr(sequence,1,70) # Trim down max length to increase efficiency
	primers <- matrix(ncol=2)[0,] # Best way to make empty matrix?

	# Storing primer sequence separate from Tm may increase efficiency
	# First, generate all primers up to Set Tm + Allowed Error. Even though it's iterating, it's faster than generating
	# all reasonably-sized oligos and finding their tms
    # This loop is the slowest part of the program

#	while ( primer.tm <= (tm+tm.errorplus) & bases <= nchar(sequence) ) {
#		bases <- bases + 1
#		primer.current <- substring(sequence,1,bases)
#		primer.tm <- finnzymes(primer.current)
#		primers <- rbind(primers,c(primer.current,primer.tm))
#	}

    # Instead, make all possible primers in from length minstart to 70
    primers <- sapply((minstart:nchar(sequence)),function(x)substr(sequence,0,x))
    # Should really bifurcate or something to get fewer tms, but for now find all tms for all oligos
    primers <- data.frame(oligo=primers,tm=sapply(primers,finnzymes))
    rownames(primers) <- NULL

    if ( primers[length(primers[,1]),2] <= (tm-tm.errorminus) ) {
        warning("Could not generate oligo that fell within desired tm range")
    }

	primers$oligo <- as.character(primers$oligo)

	# If endGC=T, find all of those which end in C or G
	if ( endGC==T ) {
		ends <- substring(primers[,1],nchar(primers[,1]),nchar(primers[,1]))
		ends.G.C <- length(grep("G",ends))+length(grep("C",ends))
		# Do they even exist?
		if (ends.G.C > 0) {
			primers <- primers[which(ends == "G" | ends == "C"),]
#			primers <- as.matrix(primers[which(ends == "G" | ends == "C"),])
		} else {
			warning("No primers in this range ending in G or C could ge generated.")
		}
	}

	if ( atleast == T ) {
		# Trim to primers above set Tm, none below
		primers <- primers[primers[,2]>tm,]
	} else {
		# Trim to those above the desired minimum (Set Tm - Allowed Error)
		primers <- primers[primers[,2]>(tm-tm.errorminus),]
	}

	# Find the primer closest to desired Tm
	error.rel <- abs(tm-as.numeric(primers[,2]))
	error.min.position <- which(error.rel==min(error.rel))
	primer.ideal <- primers[error.min.position,]

	# Add tail
	# TODO handle exceptions
	if (tail != F ) {
		primer.ideal[1] <- paste(tail,primer.ideal[1],sep="")
	}

	# Final nasty clip-down for 60-mer
	if (length(primer.ideal)>60 & sixtymer == T) {
		primer.ideal <- substring(primer.ideal,1,60)
		warning("Primer was sloppily truncated to a 60-mer. To disable this, use sixtymer=F.")
	}

    # Make the rowname be the length of the oligo
    rownames(primer.ideal) <- nchar(primer.ideal$oligo)

	return(primer.ideal)

}

designprimer.gene <- function(sequence,tails=F) {
    # convert to uppercase
    sequence <- toupper(sequence)

    if (tails==F) {
        tails <- ''
    } else if (length(tails) != 2) {
        stop('If tail is supplied, must supply 2 entries (even if one is "")')
    } 

    fwd <- designprimer(sequence,tail=tails[1])
    rev <- designprimer(revcomp(sequence),tail=tails[2])
    result <- rbind(fwd,rev)
    return(result)
}

designprimer.gene.gateway <- function(sequence,tail1=gateway1,tail2=gateway2) {
	#convert to uppercase
	primers <- designprimer.gene(sequence,tails=c(tail1,tail2))

    return(primers)
}

#Expects BioStrings data thing, not just a raw character.  I think.
designprimer.list <- function(sequences,tail1='',tail2='') {
	rawnames <- rep(names(sequences),each=2)
	primernames <- paste(rawnames,c("F","R"),sep="_")

    charfun <- function(dnastring){
        x <- as.character(dnastring)[[1]]
        return(x)
    }

    sequences <- sapply(sequences,charfun)
    names(sequences) <- NULL

	print("Designing forward and reverse Gateway cloning primers...")
    primerframe <- designprimer.gene.gateway(sequences[1],tail1=tail1,tail2=tail2)
    if (length(sequences) > 1 ) {
        for (i in 2:length(sequences)){
            cur_primers <- designprimer.gene.gateway(sequences[i],tail1=tail1,tail2=tail2)
            primerframe <- rbind(primerframe,cur_primers)
        }
    }

	#designprimer_result <- sapply(sequences,designprimer.batch.gateway)
    primerframe <- cbind(primernames,primerframe)
    colnames(primerframe) <- c('name','sequence','notes')

	return(primerframe)
}

designprimer.list.gateway <- function(sequences) {
    out <- designprimer.list(sequences,tail1=gateway1,tail2=gateway2)
    return(out)
}

write.XStringSet.path <- function(XStringObject,path="./") {
	# names of the strings
	objnames <- names(XStringObject)
	for (i in seq(along=objnames)) {
		write.XStringSet(XStringObject[i],filepath=paste(path,names(XStringObject)[i],".fasta",sep=""))
	}
}

write.IDT <- function(primertable,file=paste("IDT-primers-",format(Sys.time(), "%Y.%m.%d-%H.%M.%S"),".csv",sep="")) {
	# Defaults to writing out the primer set + date + time, to the second.
	write.table(primertable,
		    file=file,
		    quote=F,
		    row.names=F,
		    col.names=F,
		    sep=";")
}


# This is a first attempt at analyzing codon frequency/usage/etc with R. It's pretty hackish.
# Absolutely requires input to be a StringSet from the Biostrings package.
# You can get a yeast codon usage table (or cai table, in this case)
# using data(caitab) from the seqinr package.
# In that particular table for $sc, cai < .003 seems to be right for a 'bad' codon'
# Returns a sorted vector of the 'bad' codon positions within the CDS
codonAnalysis <- function(x) {
	# For now, automatically generate 'bad cai' codons using semi-arbitrary metric.
	# An attempted approximation of the results of the GenScript rare codon analysis tool results
	data(caitab)
	sctab <- data.frame(codon=rownames(caitab),cai=caitab$sc)
	badcai <- subset(sctab,cai<.003)

	x_codons <- codons(x)
	x_codon_vector <- as.character(x_codons) #trinucleotide vector of CDS

	results <- c()
	for (i in 1:length(badcai[,1])) {
		current_result <- grep(badcai[i,"codon"],x_codon_vector,ignore.case=T)
		results <- c(results,current_result)
	}

	results <- sort(results)
	return(results)

}

# Given a codon position vector, calculate the position of the ones directly next to one another

codonPos <- function(x,threshold=.003) {
	x_length <- length(x)
	x_diffs <- c()
	for (i in 1:x_length-1) {
		x_diffs_current <- x[i+1]-x[i]
		x_diffs <- c(x_diffs,x_diffs_current)
	}

	x_diffs_bool <- x_diffs==1
	x_space_1 <- as.numeric(x_diffs_bool) #A vector showing all positions where bad codons are side-by-side

	#You can get interesting repeat lengths with the 'rle' command
	x_rle <- rle(x_space_1)
	x_rle_frame <- data.frame(lengths=x_rle$lengths,values=x_rle$values)
	# Add one to the length - a length of '1' repeat = 2 codons
	x_rle_frame$lengths <- x_rle_frame$lengths+1

	# Only care about repeat locations
	x_rle_reps <- subset(x_rle_frame,values==1)
	x_rle_reps_levels <- levels(factor(x_rle_reps$lengths))
	x_rle_reps_levels_n <- length(x_rle_reps_levels)

	x_rle_reps_final <- data.frame(reps=rep(NA,x_rle_reps_levels_n),count=rep(NA,x_rle_reps_levels_n))
	for (i in 1:x_rle_reps_levels_n) {
		current_rep_index <- which(x_rle_reps$lengths==x_rle_reps_levels[i])
		current_rep_sum <- sum(x_rle_reps$values[current_rep_index])
		x_rle_reps_final[i,"reps"] <- x_rle_reps_levels[i]
		x_rle_reps_final[i,"count"] <- current_rep_sum
	}

	return(x_rle_reps_final) # Table of bad codon repeat length vs. counts of those repeats
#	return(x_rle_frame) # rle data table
	return(x_space_1) # Atomic vector of bad codon repeat positions
	return(x_diffs) # If you want the difference vector
}

# plotRanges:
# Designed to plot 'IRanges' objects, derived from IRanges manual

plotRanges <- function(x, xlim = x, main = deparse(substitute(x)),
	col = "black", sep = 0.5, ...) {

	xlim = c(start(x)[1],end(x)[length(end(x))])

	height <- 1
	if (is(xlim, "Ranges"))
	xlim <- c(min(start(xlim)), max(end(xlim)))
	bins <- disjointBins(IRanges(start(x), end(x) + 1))
	plot.new()
	plot.window(xlim, c(0, max(bins) * (height + sep)))
	ybottom <- bins * (sep + height) - height
	rect(start(x) - 0.5, ybottom, end(x) + 0.5, ybottom +
	height, col = col, ...)
	title(main)
	axis(1)
}


gen_overlap <- function(sequence_in,oligo_length=120,tm=63,tm.errorplus=5,minstart=20,overhang="5") {
    library("gtools")

	# Script expects a string (ATGC) input
	sequence_in <- toupper(sequence_in)

	# A particular manipulation of 'designprimer()' is used repeatedly to find good regions for annealing:
	findoverlap <- function(seq_in) {
		overlap <- as.character(revcomp(designprimer(revcomp(seq_in),
							tm=tm,
							tm.errorplus=tm.errorplus,
							endGC=F,
							atleast=T,
							minstart=minstart)[1,1]))
		return(overlap)
	}

	# Iterative method - does not plan things out in advance, just moves along the sequence designing overlaps
	# Make the first oligo

	iter=1
	print("Designing oligos... :")
	cat(paste(iter,""))
	current_pos=1
	current_120 <- substr(sequence_in,start=current_pos,stop=oligo_length)
	current_overlap <- findoverlap(current_120)
	current_pos=current_pos+oligo_length-nchar(current_overlap)

	oligoresult <- data.frame(oligo=current_120,overlap=current_overlap)

	# Generate the rest
	while ( nchar(substr(sequence_in,start=current_pos,stop=current_pos+oligo_length-1)) == oligo_length ) {

		iter=iter+1
		cat(paste(iter,""))
		current_120 <- substr(sequence_in,start=current_pos,stop=current_pos+oligo_length-1)
		current_overlap <- findoverlap(current_120)
		current_pos <- current_pos+oligo_length-nchar(current_overlap)
		current_result <- data.frame(oligo=current_120,overlap=current_overlap)
		oligoresult <- rbind(oligoresult,current_result)
	}

	# Generate the last oligo
	iter=iter+1
	cat(paste(iter,""))
	cat("\n")
	current_120 <- substr(sequence_in,start=current_pos,stop=nchar(sequence_in))
	current_result <- data.frame(oligo=current_120,overlap=NA)

	# note: 'overlap' is meaningless for this last oligo but I can't make it 'NA' without breaking 'finnzymes()'. Fix finnzymes().

	oligoresult <- rbind(oligoresult,current_result)

	# I don't know why it's using some other data type, but this makes the oligos back into characters
	oligoresult$oligo <- as.character(oligoresult$oligo)
	oligoresult$overlap <- as.character(oligoresult$overlap)

	# Calculate overlap tm
	oligoresult$tm <- unlist(lapply(oligoresult$overlap,finnzymes))

	# Check to see whether there's an even number
	if ( even(length(oligoresult[,1])) ) {
		print("Even number of oligos generated")
	} else {
		cat("Tms: ")
		cat(oligoresult[,3])
		cat("\n")
		stop("Generated odd number of oligos - change Tm")
	}


	# Reverse complement every other oligo
	n <- length(oligoresult[,1])/2
	pos <- 2*(1:n)

	if ( overhang == "3" ) {
		pos <- pos-1 # make odd numbers
	} else if (overhang != "5" ) {
		stop('overhang parameter must be "5" or "3"')
	}

	for (i in pos) {
		oligoresult[i,"oligo"] <- revcomp(oligoresult[i,"oligo"])
	}


	oligoresult[length(oligoresult[,1]),2] <- NA
	oligoresult[length(oligoresult[,1]),3] <- NA

	return(oligoresult)

}

# converts ATGC string into reverse complement
revcomp <- function(sequence_in) {
	# Convert Biostrings/other characters into standard format
	sequence <- qa.atgc(sequence_in)

	# Reverse sequence
	sequence <- strsplit(sequence,split="")[[1]]
	sequence <- rev(sequence)

	# Substitute
	seq_in <- sequence
	sequence[seq_in=="A"] <- "T"
	sequence[seq_in=="T"] <- "A"
	sequence[seq_in=="G"] <- "C"
	sequence[seq_in=="C"] <- "G"

	sequence <- paste(sequence,collapse="")
	return(sequence)
}

# Finds largest contiguous non-N read in a sequencing result
trimSequencing <- function(seqs) {
    # Read in files in directory
#    seqs <- read.DNAStringSet.path(indir)
    seqs_trimmed <- rep(NA,length(seqs))
    # Tack on an 'N' on each end (hack) so that if sequencing continues to end, script still works
    seqs_fixed <- DNAStringSet(paste("N",as.character(seqs),sep=""))
    seqs_fixed <- DNAStringSet(paste(as.character(seqs_fixed),"N",sep=""))
    names(seqs_fixed) <- names(seqs)
    seqs <- seqs_fixed
    for (i in 1:length(seqs)) {
        split_seq <- strsplit(as.character(seqs[[i]]),split="")[[1]]
        n_locs <- which(split_seq=="N")
        if (length(n_locs)!=0) {
            max_loc <- which(max(diff(n_locs))==diff(n_locs))
            max_start <- n_locs[max_loc]
            max_stop <- n_locs[max_loc+1]
            best_seq <- paste(split_seq[(max_start+1):(max_stop-1)],collapse="")
            seqs_trimmed[i] <- best_seq
        } else {
            seqs_trimmed[i] <- as.character(seqs[i])
        }
    }
    seq_set <- DNAStringSet(seqs_trimmed)
    names(seq_set) <- names(seqs)
    return(seq_set)
}

# Designed to work with a 'list' object with 'stringset' as saved DNA sequences used in an alignment.
# Defaults to "out.fasta" in one's main directory (My Documents for Windows, /Users/YourUser for Mac).
# To specify your own file, add file="your/file.fasta" (e.g.)
write.alignment <- function(alignment,file="~/out.fasta") {
    write.XStringSet(alignment$stringset,format="fasta",file=file)
}

#GibsonVol <- function(conc,bp,volume=5) {
    # volume is in uL, default is 5
    # conc is a vector of concentrations
    # bp is a vector of corresponding base pairs

    # x1 + x2 + ... + x(n-1) + xn = volume
    # (x1*conc1)/bp1 = (x2*conc2)/bp2 = ... = (xn*concn/bpn)
    # solve!
#    if (length(conc)!=length(bp) {
#        stop("concentration vector must be same length as base pair vector")
#    }

#    firstcol <- c(1,
#}


### randGS
### one-off program to make random GS linkers using most common codons.

randgs <- function() {
    # top three glycine codons - the fourth is lowest-usage
    gly_codons <- c("GGA","GGT","GGC")
    # top four serine codons - other two are lower-usage
    ser_codons <- c("AGT","TCA","TCT","TCC")

    desired_seq <- c("GSGSGSGSGS")
    dna_seq <- c()

    # iterate through that desired sequence, pick random codon
    for (i in 1:nchar(desired_seq)) {
        if(substr(desired_seq,i,i)=="G") {
            # random integer, 1-3
            rnum <- ceiling(runif(1,0,3))
            dna_seq <- c(dna_seq,gly_codons[rnum])
        } else {
            # random number
            rnum <- ceiling(runif(1,0,4))
            dna_seq <- c(dna_seq,ser_codons[rnum])
        }
    }

    dna_seq <- paste(dna_seq,collapse="")
    return(dna_seq)
}

# script to hackishly optimize GS linkers using nupack against a similar possible partner
# i.e. we want to Gibson together pieces with flanking GS linkers, but they are so similar that 
# undesired pairing and lose the insert. Want to design less-matchy linkers with the same coding sequence.
# n determines how many times to repeat this
# TODO: make this more generalized. Need to write up preferred codon table (complete) to allow arbitrary input.
# could be useful for any two target sequences that are the ends of Gibsons and very similar.
optgs <- function(n=10) {
    ### user-defined stuff:
    # Target
    target <- "ggttctggatcaggtagtggctcaggatct"
    target <- toupper(target)

    nupackdir <- "/home/nick/nupack/nupack3.0"
    outfolder <- "/home/nick/nupackout" # must exist beforehand (fixme later)

    ### methods/iterator
    # top three glycine codons - the fourth is lowest-usage
    gly_codons <- c("GGA","GGT","GGC")
    # top four serine codons - other two are lower-usage
    ser_codons <- c("AGT","TCA","TCT","TCC")

    ### step 1: generate 10 random 10X gs linkers
    linker <- vector(length=n)
    for (i in 1:n) {
        linker[i] <- randgs()
    }

    ### step 2: Using NUPACK, find the one with the lowest expected concentration versus the target

    # set up NUPACK environment. Works in linux, probably mac, probably not windows.
    # Actually, all commands need to be a one-liner or something (horrible), so need to export the dir before every command
    nupackprefix <- paste("export NUPACKHOME=",nupackdir,sep="")

    # Make new dir in outfolder named after the exact time (prevents overwriting runs)
    outfolder_sub <- paste(outfolder,"/",format(Sys.time(), "%Y-%m-%d--%H-%M-%S"),sep="")
    dir.create(outfolder_sub)

    # Make the dirs for input/workfiles for nupack
    dir.create(paste(outfolder_sub,"/workfiles",sep=""))

    scores <- vector(length=n)

    for (i in 1:n) {
        if (n > 100 & i%%(n%/%25)==0) {
            perc_done <- format(i/n*100,digits=3)
            message(perc_done,"% complete.")
        }

        dir.create(paste(outfolder_sub,"/workfiles","/",i,sep=""))
        # Prepare input file
        inputfile <- paste("2\n",target,"\n",revcomp(linker[i]),"\n2",sep="")

        writeChar(inputfile,paste(outfolder_sub,"/workfiles","/",i,"/input.in",sep=""),eos=NULL)

        # run 'complexes' from NUPACK, generating .ocx file (is ocx vs cx preferable?)
        complexescommand <- paste(nupackprefix," && ","cd ",outfolder_sub,"/workfiles"," && ",nupackdir,
                "/bin/complexes -T 50 -material dna -mfe ",outfolder_sub,"/workfiles/",i,"/input",sep="")
        system(complexescommand,ignore.stdout=T)

        # Define concentrations file
        inputfile_conc <- "0.5e-6\n0.5e-6"
        writeChar(inputfile_conc,paste(outfolder_sub,"/workfiles","/",i,"/input.con",sep=""),eos=NULL)

        # run 'concentrations' from NUPACK, generating .eq file (equilibrium concentrations)
        concentrationscommand <- paste(nupackprefix," && ","cd ",outfolder_sub,"/workfiles"," && ",nupackdir,
            "/bin/concentrations -sort 3 ",outfolder_sub,"/workfiles/",i,"/input",sep="")
        system(concentrationscommand,ignore.stdout=T)

        # hackish (probably) way to get last 5 lines of output file, convert to data frame
        results_i <- as.character(readLines(paste(outfolder_sub,"/workfiles/",i,"/input.eq",sep="")))
        results_i <- results_i[(length(results_i)-4):length(results_i)]
        results_i <- data.frame(matrix(unlist(strsplit(results_i,split="\t")),ncol=5,byrow=T),stringsAsFactors=F)
        colnames(results_i) <- c("c_int","s1","s2","dG","conc")

        # Add up the score - right now this only works for comparing 2 sequences, as it expects 5 possible combinations.
        # The ones we want to minimize are the self-self of the tested oligo and the self-target equilibrium concentrations.
        # These are rows 5 and 4 respectively
        scores[i] <- sum(as.numeric(results_i[4:5,5]))
    }

    # report the winner (lowest score)
    lowest_pos <- which(scores==min(scores))
    lowest <- linker[lowest_pos]
    result <- data.frame(sequence=lowest,score=min(scores))
    rownames(result) <- lowest_pos
    return(result)
}


# Make all possible 10XGS sequences, then run through NUPACK?
# General way to write 'em out: repeat a1 n/length(a1) times. Repeat a2 n/(length(a1)*length(a2)) times. etc etc.
# simpler: repeat a_i by: n/cumprod(length(a_n))[i]
# first example is just for 10X GS linker
# warning: it's very very very big and will take forever to display raw table. 
# Should probably hide it in a class to prevent easy printing.

allgs <- function() {

    # top three glycine codons - the fourth is lowest-usage
    gly <- c("GGA","GGT","GGC")
    # top four serine codons - other two are lower-usage
    ser <- c("AGT","TCA","TCT","TCC")
    input <- list(gly,ser,gly,ser,gly,ser,gly,ser,gly,ser)

    n <- 3^5*4^5
    lengths <- rep(c(length(gly),length(ser)),5)
    lengths_cumu <- cumprod(lengths)
    big_table <- data.frame(matrix(ncol=10,nrow=lengths_cumu[length(lengths)]))
    for (i in 1:length(lengths)) {
        current_vec <- rep(input[[i]],each=n/lengths_cumu[i])
        big_table[,i] <- current_vec
    }
    out_strings <- apply(big_table,1,paste,collapse="")
    return(out_strings)
}

# Takes the list returned by alignseqs and plots an alignment
# TODO:
#	1) Remove IRanges dependency? Is the disjoint binning calculation really that fancy?
#	2) Due to using fill="black" for the coverage rectangle, it doesn't show up in the key. Fix this!
#	3) Add sequencing result length info somehow
alignplot <- function(list_in) {
	library(ggplot2)
	library(IRanges)
	# rename the list's alignment to keep things simple
	x <- list_in$alignment
	# Get the 'start' and 'width' info from each pairwise alignment
	# the 'range' accessor doesn't seem to work, so used direct attributes
	x_ranges_raw <- attributes(subject(x))$range
	x_ranges <- data.frame(start=start(x_ranges_raw),end=end(x_ranges_raw),width=width(x_ranges_raw))
	# Get the name of each sequencing result
	x_ranges$name <- names(unaligned(pattern(x)))
	# Calculate disjoint binning to minimize number of levels needed
	x_ranges$bin <- disjointBins(IRanges(x_ranges$start,x_ranges$end))

	# Add in full coverage as a single data type thing
	x_coverage <- x_ranges
	x_coverage$name <- "coverage"
	x_coverage$bin <- max(x_ranges$bin)+1

	# Get full size of reference sequence so that xlim can be set on plot
	ref_width <- width(list_in$stringset[names(list_in$stringset)=="reference"])

	# mismatch insertion deletion (MID) table. Remember to add insertions and deletions...

	mismatches <- c()
	insertions <- c()
	deletions <- c()

	# mismatches
	if (length(row.names((mismatchTable(x))))!=0) {
		mismatches <- data.frame(position=mismatchTable(x)$SubjectStart,number=mismatchTable(x)$PatternId,type="mismatch")
		mismatches$name <- x_ranges$name[mismatches$number]
	}
	# insertions
	if (length(as.data.frame(insertion(x))[,1])) {
		insertions <- as.data.frame(insertion(x))
		insertions$name <- x_ranges$name[insertions$space]
		colnames(insertions) <- c("number","position","stop","width","name")
		insertions$type <- "insertion"
		# fix positions from relative to pattern to relative to subject
		for (i in x_ranges$name) {
			insertions[insertions$name==i,"position"] <- insertions[insertions$name==i,"position"] + x_ranges[x_ranges$name==i,"start"]
		}
	}
	
	# deletions
	if (length(as.data.frame(deletion(x))[,1])) {
		deletions <- as.data.frame(deletion(x))
		deletions$name <- x_ranges$name[deletions$space]
		colnames(deletions) <- c("number","position","stop","width","name")
		deletions$type <- "deletion"
		# fix positions from relative to pattern to relative to subject
		for (i in x_ranges$name) {
			deletions[deletions$name==i,"position"] <- deletions[deletions$name==i,"position"] + x_ranges[x_ranges$name==i,"start"]
		}
	}

	# columns to keep from each result
	keep_cols <- c("number","position","type","name")
	mid <- rbind(mismatches[,keep_cols],insertions[,keep_cols],deletions[,keep_cols])
	if (length(mid)!=0) {
		# Always use same shapes for mismatches/insertions/deletions, and show key
		mid$type <- factor(mid$type,levels=c("mismatch","insertion","deletion"))
		# Apply bins from before so they plot in the right place:
		mid$bin <- NA

		for (i in levels(factor(mid$name))) {
			mid[mid$name==i,"bin"] <- x_ranges[x_ranges$name==i,"bin"]
		}		
	}
		
	# Fix factor levels to make it plot 'coverage' on top, (ascending order)
	x_ranges$bin <- factor(x_ranges$bin,rev(levels(factor(x_ranges$bin))))
	x_coverage$bin <- factor(x_coverage$bin,rev(levels(factor(x_coverage$bin))))
	
	# Plot
	p <- qplot(data=x_ranges,xmin=start,xmax=end,ymin=0,ymax=10,fill=factor(name),geom="rect") +
	geom_rect(data=x_coverage,fill=I("black"),aes(xmin=start,xmax=end,ymin=0,ymax=10)) +
	facet_grid(bin~.) + 
	theme_bw() + 
	opts(axis.ticks=theme_blank(),axis.text.y=theme_blank(),axis.title.y=theme_blank(),axis.title.x=theme_text(vjust=0,size=12)) + 
	opts(strip.text.y=theme_blank(),strip.background=theme_blank(),title="Sequence Coverage") + 
	scale_fill_discrete("Result") +
	scale_shape_discrete("Disrepancies") +
	xlab("Reference Sequence Position") + 
	xlim(c(0,1.1*ref_width)) + 
	scale_x_continuous(breaks=1e3*(0:100)/2) + 
	coord_cartesian(xlim=c(0,ref_width))

	if( is.list(mid) ) {
		p + geom_point(data=mid,aes(x=position,xmin=NULL,xmax=NULL,ymin=NULL,ymax=NULL,y=5,shape=factor(type)),fill=I("black"))
	} else {
		p
	}
}

# alignwrite - writes out alignment svg+png and fasta files of input + aligned sequences
# TODO: append the date/time to each file to avoid overwriting+add useful info?
alignwrite <- function(alignseqs_alignment,name=NULL,svg=F) {
	library(Biostrings)
	if (is.null(name)) {
		name="last-alignment"
	} 

	message("Writing files to:")

	# reassign alignment for simplicity
	x <- alignseqs_alignment

	# Write out input sequences to a single fasta file
	message("\t~/",name,"-input.fasta")
	write.XStringSet(x$stringset,file=paste("~/",name,"-input.fasta",sep=""))

	# Write out aligned sequences to a single fasta file
	x_aligned <- aligngappedseq(x$alignment)
	write.XStringSet(x_aligned,file=paste("~/",name,"-alignment.fasta",sep=""))
	message("\t~/",name,"-alignment.fasta")

	# Write out summary plot for visualization
	x_plot <- alignplot(x)
	if (svg == T) {
		ggsave(x_plot,height=4,width=12,file=paste("~/",name,".svg",sep=""))
		message("\t~/",name,".svg")
	}
	ggsave(x_plot,height=4,width=12,file=paste("~/",name,".pdf",sep=""))
	message("\t~/",name,".pdf")

	ggsave(x_plot,height=4,width=12,dpi=90,file=paste("~/",name,".png",sep=""))
	message("\t~/",name,".png")
}

# aligngappedseq: takes pairwise alignment (Biostrings), constructs a gapped StringSet (collection to write to fasta)
aligngappedseq <- function(aligned_in) {
	library(Biostrings)

	# Simplify alignment name
	x <- aligned_in
	
	# Get gapped representation of sequencing results. Does not include reference sequence (subject)
	# It also manages to forget the names of the results
	x_aligned <- aligned(x,degap=F)
	x_aligned_names <- names(unaligned(pattern(x)))

	# Get the reference sequence:
	x_subject <- unaligned(subject(x))
	x_subject_name <- names(unaligned(subject(x)))

	# Since we're using separate tools to get subject + aligned sequences, make sure they're the same length
	if (any(width(x_aligned)!=width(x_subject))) {
		stop("The script messed up - aligned sequences are of different length than reference sequence.")
	}

	# Combine the subject and aligned sequences:
	x_subject_char <- as.character(x_subject)
	x_aligned_char <- as.character(x_aligned)
	x_writable <- DNAStringSet(c(x_subject_char,x_aligned_char))
	names(x_writable) <- c(x_subject_name,x_aligned_names)
	return(x_writable)
}

### new_gen_oligos:
### generates oligos from an input sequence, good for assembly of genes
### from an oligo order

new_gen_oligos <- function(sequence_input,tm=65,oligo_size=120,forceeven=F,start_5=T) {
	# sequence_input: must be plain text gene sequence
	# tm: default overlap Tm min is 65°C
	# oligo_size: max size of the oligo - this script makes all oligos this size
	# forceeven: Sets whether there must be an even number of oligos,
	#       which forces the ends of the assembly construct to 
	#		have the same directionality (5' or 3')
	# start_5: Determines which oligos to reverse complement - 
	#		 if even (start_5=T), then first oligo overhang is 5'
	#		 if odd (start_5=F), the first oligo overhang is 3'
	# even=T and start_5=F would make an even number of oligos,
	# the overlapped construct of which would have 3' overhangs.
	
	
	# This ratio gets reused - sequence length / oligo max size
	seq_oligo_ratio <- nchar(sequence_input)/oligo_size

	# Figure out minimum possible number of oligos needed
	oligo_n <- ceiling(seq_oligo_ratio)

	# if 'even' setting is true, make sure that oligo_n is even
	oligo_n_isodd <- oligo_n%%2==1
	if (forceeven==T&oligo_n_isodd) {
		oligo_n <- oligo_n+1
	}
	# Set initial overlap tm to arbitrary low number, forcing random choice
	o_f <- data.frame(overlap_tm=1)

	expansion_loop <- function(sequence_input,oligo_n) {
		# Run a loop that generates oligos of length 'oligo_size' 'fairly' using
		# the finnzymes calculator
		# TODO: write this as a function
		# function starts here

		message("Trying with ",oligo_n," oligos.")

		# determine the starting location of each overlap - discretize seq length / oligo_size * n_vector
		cont_overlap_locations <- 1:(oligo_n-1)*(nchar(sequence_input)/oligo_n)
		disc_overlap_locations <- ceiling(cont_overlap_locations)

		# initial vector of overlaps
		overlaps <- c()
		for (i in disc_overlap_locations) {
			overlaps <- c(overlaps,substr(sequence_input,i,i))
		}
        
        # initial vector of oligos
		oligos <- substr(sequence_input,0,disc_overlap_locations[1])
		if (length(overlaps)>1) {
			for (i in 1:(length(overlaps)-1)) {
				current_oligo <- substr(sequence_input,disc_overlap_locations[i],disc_overlap_locations[i+1])
				oligos <- c(oligos,current_oligo)
			}
		}
		final_oligo <- substr(sequence_input,
							  disc_overlap_locations[length(disc_overlap_locations)],
							  nchar(sequence_input))
		oligos <- c(oligos,final_oligo)

		# data frame to hold start/stop/overlap/oligo information
		o_f <- data.frame( id = 1:oligo_n,
						   start=c(0,disc_overlap_locations),
						   stop=c(disc_overlap_locations,nchar(sequence_input)),
						   oligo=as.character(oligos),
						   overlap=as.character(c(overlaps,NA)),
						   overlap_tm=as.numeric(NA),
						   maxed=F
						 )
		o_f$overlap <- as.character(o_f$overlap)
		o_f$oligo <- as.character(o_f$oligo)
		total_oligos_len <- sum(nchar(o_f$oligo))

		# Prioritizing those with the lowest Tm
		overlap_nonmaxed <- which(rep(T,oligo_n-1))
		o_f$overlap_tm <- sapply(o_f$overlap,function(x)as.numeric(finnzymes(x)))
		#o_f$overlap_tm[overlap_nonmaxed] <- sapply(o_f$overlap[overlap_nonmaxed],function(x)as.numeric(finnzymes(x)))
		oligo_changed <- 1 #arbitrary working intialization
		j <- 0
		while (total_oligos_len < oligo_size*oligo_n) {
			# Find available overlap with the lowest Tm
			smallest_overlap <- min(o_f$overlap_tm[overlap_nonmaxed],na.rm=T)
			min_location <- which(o_f$overlap_tm==smallest_overlap)

			# If two or more Tms are equal (as they will be at first), pick one at random
			if (length(min_location)>1) {
				min_location <- sample(min_location,1)
			}

			# First, see if we've met the oligo size threshold
			if (nchar(o_f$oligo[min_location])==oligo_size &
			nchar(o_f$oligo[min_location+1])==oligo_size) {
				# if both hit the threshold, set 'maxed' column to true
				new_slice <- o_f
			} else if (nchar(o_f$oligo[min_location])==oligo_size) {
				# if the first hits the threshold, increase size of the second
				# (decrement overlap start site))
				# also set 'maxed' variable to true
				o_f$start[min_location+1] <- o_f$start[min_location+1]-1
			
				new_slice <- gen_new_o_f_slice(sequence_input,o_f,min_location,lr="r")
			} else if (nchar(o_f$olig[min_location+1])==oligo_size) {		
				# if the second hits the threshold, increase the size of the first
				# (increment overlap stop site)
				# also set 'maxed' variable to true
				o_f$stop[min_location] <- o_f$stop[min_location]+1
			
				new_slice <- gen_new_o_f_slice(sequence_input,o_f,min_location,lr="l")
			} else {
				# increase the size of the smaller one or arbitrarily if they're the same size
				if (nchar(o_f$oligo[min_location])==nchar(o_f$oligo[min_location+1])) {
					# If they're the same size, increase the first one's size
					o_f$stop[min_location] <- o_f$stop[min_location]+1

					new_slice <- gen_new_o_f_slice(sequence_input,o_f,min_location,lr="l")
				} else if (nchar(o_f$oligo[min_location])<nchar(o_f$oligo[min_location+1])) {
					# If Second is bigger, increase size of first
					o_f$stop[min_location] <- o_f$stop[min_location]+1

					new_slice <- gen_new_o_f_slice(sequence_input,o_f,min_location,lr="l")
				} else {
					# First is bigger, so increase size of second
					o_f$start[min_location+1] <- o_f$start[min_location+1]-1

					new_slice <- gen_new_o_f_slice(sequence_input,o_f,min_location,lr="r")
				}
			}
			oligo_changed <- min_location
			# Put the new slice back into the data frame
			o_f <- new_slice

			o_f$maxed <- nchar(o_f$oligo)==oligo_size

			# Record overlaps that have maxed out the oligo size
			overlap_nonmaxed <- sapply(1:(length(o_f$maxed)-1),
											function(x) {
											o_f$maxed[x]==F|
											o_f$maxed[x+1]==F})

			total_oligos_len <- sum(nchar(o_f$oligo))
			# Find the Tm of only the overlaps that changed
			o_f$overlap_tm[oligo_changed] <- finnzymes(o_f$overlap[oligo_changed])

			#print(o_f)

			j <- j+1
			if (j%%10==0) {
				#message("Oligo Lengths: ",paste(nchar(o_f$oligo),collapse=" "))
			}
		}
		return(o_f)
	}

	# Apply setting of 'even' argument to oligo number increment size
	if (forceeven==T) {
		oligo_increment <- 2
	} else if (forceeven==F) {
		oligo_increment <- 1
	} else {
		stop("Invalid 'even' setting, must be T or F")
	}

	# Loop to make sure the minimum threshold is hit
	j <- 0
	while((any(o_f$overlap_tm < tm,na.rm=T))) {
		if (j > 0) {
			oligo_n <- oligo_n+oligo_increment
		}
		j <- j+1
		o_f <- expansion_loop(sequence_input,oligo_n)
	}

	# Reverse complement every other oligo depending on 'start_5' argument
	if (start_5==T) {
		torev_comp <- which(1:length(o_f$oligo)%%2==0)
	} else if (start_5==F) {
		torev_comp <- which(1:length(o_f$oligo)%%2==1)
	} else {
		stop("Invalid start_5 option (Must be T or F)")
	}

	o_f$oligo[torev_comp] <- sapply(o_f$oligo[torev_comp],function(x)revcomp(x))

	# Pare down the data frame before returning it
#	o_f$id <- NULL
	o_f$start <- NULL
	o_f$stop <- NULL
	o_f$maxed <- NULL

	return(o_f)
}

# function to regenerate pieces of oligo_frame on the fly given new info
gen_new_o_f_slice <- function(sequence_input,o_f_in,overlap_id,lr=F) {

	nstart <- o_f_in$start
	nstop <- o_f_in$stop

	# generate new overlap
	new_overlap <- substr(sequence_input,nstart[overlap_id+1],nstop[overlap_id])
		
	# generate new oligo
	if (lr=="l") {
		id=overlap_id
		new_oligo <- substr(sequence_input,nstart[overlap_id],nstop[overlap_id])
	} else if (lr=="r"){
		id=overlap_id+1
		new_oligo <- substr(sequence_input,nstart[overlap_id+1],nstop[overlap_id+1])
	} else {
		stop("didn't supply 'l' or 'r'")
	}

	o_f_out <- o_f_in

	o_f_out[overlap_id,"overlap"] <- new_overlap
	o_f_out[id,"oligo"] <- new_oligo

	return(o_f_out)
}










