#' @export

simpleImport <- function(bamFile, testRanges, samplename=NULL, 
                         #region options
                         style="region",
                         method="bin", binMethod = "mean", nOfWindows=100, 
                         expand = NULL, flank=NULL, flankUp=NULL, flankDown=NULL, 
                         #import options
                         format="bam", paired=FALSE, removeDup=FALSE, FragmentLength=150, 
                         forceFragment=NULL, minFragmentLength=NULL, maxFragmentLength=NULL,
                         normalize="RPM", downSample=NULL,
                         strand.aware = FALSE,
                         seqlengths=NULL, 
                         #pwm options
                         genome=NULL, cutoff=80){
  ### TO DO: parameter checking here
  check_params(style, c("region", "centered", "flank"))
  if (!is.null(method)){
    check_params(method, c("bin", "spline"))
  }
  format <- tolower(format)
  check_params(format, c("bam", "bigwig", "bw", "pwm", "rlelist", "granges"))
  
  ## If format is bam, read header and get contig information
  if(format == "bam"){
    ## Get all chromosomes in bamFile
    message("Reading Bam header information...",appendLF = FALSE)
    allchrs <- names(scanBamHeader(bamFile)[[1]]$targets)
    lengths <- as.vector(scanBamHeader(bamFile)[[1]]$targets)
    names(lengths) <- allchrs
    message("..Done")
  }
  
  ## For remaining formats import data and find contig information
  # Import bigwig.
  if(format %in% c("bigwig", "bw")){
    message("Importing BigWig...",appendLF = FALSE)
    genomeCov <- import.bw(bamFile, as = "RleList")
    if(is.null(seqlengths)){
      seqlengths(genomeCov) <- unlist(lapply(genomeCov,length))
    }else{
      allchrs <- intersect(names(seqlengths), names(genomeCov))
      if (length(allchrs)==0){ error("No overlapping chromosomes between coverage and supplied seqlengths") }
      genomeCov <- genomeCov[allchrs] #subset to get only desired chrs
      seqlengths(genomeCov)[allchrs] <- seqlengths[allchrs] #set seqlengths
    }
    lengths <- seqlengths(genomeCov)
    allchrs <- names(lengths)
    message("..Done")
  }
  
  # If format is pwm, calculate motifs on forward and reverse strand.
  if(format=="pwm"){
    bamFile <- pwmToCoverage(bamFile,genome,min=cutoff,removeRand=FALSE)
    format <- "rlelist"   
  }
  
  # If format is granges, simply covert to coverage (This would work for GenomicInterval or GenomicAlignments)
  if(format=="granges"){
    genomeCov <- coverage(bamFile)
    format <- "rlelist"   
  }  
  
  # If format is rlelist, import rle and set widths/contigs by seqlengths.
  if(format=="rlelist"){
    message("Importing rlelist",appendLF = FALSE)
    genomeCov <- bamFile
    if(is.null(seqlengths)){
      seqlengths(genomeCov) <- unlist(lapply(genomeCov,length))
    }else{
      allchrs <- intersect(names(seqlengths), names(genomeCov))
      if (length(allchrs)==0){ stop("No overlapping chromosomes between coverage and supplied seqlengths") }
      genomeCov <- genomeCov[allchrs] #subset to get only desired chrs
      seqlengths(genomeCov)[allchrs] <- seqlengths[allchrs] #set seqlengths
    }
    lengths <- seqlengths(genomeCov)
    allchrs <- names(lengths)
    message("..Done")
  }

  original_ranges <- testRanges #save for later!
  testRanges$original_order <- paste0("range_", 1:length(testRanges))

  # Expand ranges as specified by 'style'
  testRanges <- make_ranges(testRanges, style = style, expand = expand, 
                            flank = flank, flankUp = flankUp, flankDown = flankDown)
  
  # Exclude and count regions which when extended are outside contig boundaries.
  message("Filtering regions which extend outside of genome boundaries...",appendLF = FALSE)
  chr_lengths <- lengths[as.vector(seqnames(testRanges))]
  oob_idx <- which(end(testRanges) > chr_lengths | start(testRanges) < 0)
  message("..Done")
  if (length(oob_idx) > 0){
    testRanges <- testRanges[-oob_idx]
    message(paste("Filtered regions with indexes:", paste(oob_idx, collapse= ", ")))
  }
  message("Filtered ", length(oob_idx)," of ",length(original_ranges)," regions")
  
  # if style is flank, flanks are specified in bp, and  regions are not of equal width,
  # split ranges out into region plus left and right flanks
  # otherwise starts and ends of regions will not line up in matrix
  # e.g. 100 bp flank with 800 bp region will have 10/80/10 bins
  # while 100 bp flank with 300 bp region will have 20/60/20 bins
  # other cases are fine without splitting
  
  to_split <- (style == "flanked") & 
    (length(unique(width(testRanges))) > 1) & 
    (mode(flank) == "numeric" | (mode(flankUp) == "numeric" & mode(flankDown) == "numeric"))
                              
  if (to_split){
    ranges_list <- split_ranges(testRanges, flank, flankUp, flankDown)
  } else {
    ranges_list <- list(testRanges)
  }
  
  # if regions are not of equal width and method is not 'bin' or 'spline', 
  # set method to 'bin' with 100 windows
  
  if (length(unique(width(testRanges))) > 1 & is.null(method)){
    warning("Cannot make matrix with regions of varying width: setting method to 'bin'")
    method <- "bin"
  }
 
  if (!is.null(method) && method == "bin") {
    # then filter any regions that are too small
    # Exclude regions which are smaller than the number of windows to be used
    message("Filtering regions which are smaller than windows into region...",appendLF = FALSE)
    
    too_small_idx <- unique(unlist(lapply(ranges_list, function(x){
      which(width(x) < nOfWindows)
    })))
    
    ranges_list <- lapply(ranges_list, function(x){ x[-too_small_idx]})
    
    message("..Done")
    message(paste0("Filtered ", length(too_small_idx), " regions that are smaller than nOfWindows (", 
                   nOfWindows, "bp)"))
    if (length(too_small_idx) > 0){
      message(paste("Filtered regions:", paste(too_small_idx, collapse= ", ")))
    }
  }
  
  # then read in data and create matrices of signal, per contig
  
  ## Create GRanges to be used in scanBamParam while reading in Bamfile regions.
  reducedExtTestRanges <- reduce(testRanges)

  ## Set up scanBanParam for reading in bam file.
  
  ## if format is bam read in bamfile.
  if(format == "bam"){
    
    if(!removeDup){
      Param <- ScanBamParam(which=GRanges(seqnames=seqnames(reducedExtTestRanges[seqnames(reducedExtTestRanges) %in% allchrs]),IRanges(start=start(reducedExtTestRanges[seqnames(reducedExtTestRanges) %in% allchrs]),end=end(reducedExtTestRanges[seqnames(reducedExtTestRanges) %in% allchrs]))))
    }else{
      Param <- ScanBamParam(flag=scanBamFlag(isDuplicate=FALSE),which=GRanges(seqnames=seqnames(reducedExtTestRanges[seqnames(reducedExtTestRanges) %in% allchrs]),IRanges(start=start(reducedExtTestRanges[seqnames(reducedExtTestRanges) %in% allchrs]),end=end(reducedExtTestRanges[seqnames(reducedExtTestRanges) %in% allchrs]))))
    }
    
    message("Reading tags from ",bamFile,appendLF=FALSE)
    #totalReads <- alignmentStats(bamFile)[,"mapped"]
    # QuasR is breaking builds so for now relies on USER input
    totalReads <- 10^6
  ##  if data is single end reads then import reads and reset reads to fragment length.
  ##  Calculate fragment length from cross-coverage if not provided
  
    if(paired==FALSE){
      total <- readGAlignments(bamFile,param=Param)
      message("..Done.\nRead in ",length(total)," reads")
      
      if(is.null(FragmentLength)){
        FragmentLength <- getShifts(total,lengths,shiftWindowStart=1,shiftWindowEnd=400)
      }
      
      message("Extending reads to fragmentlength of ",FragmentLength,appendLF=FALSE)
      temp <- resize(as(total,"GRanges"),FragmentLength,"start")
      message("..done")
    }

    ##  if data is paired end reads then import reads.
    ##  Reset fragment size and/or filter to fragment length range is specified.
  
    if(paired==TRUE){
      
      gaPaired <- readGAlignments(bamFile, param=ScanBamParam(what=c("mpos"),
                                                              flag=scanBamFlag(isProperPair = TRUE,isFirstMateRead = TRUE)))      
      tempPos <- GRanges(seqnames(gaPaired[strand(gaPaired) == "+"]),
                         IRanges(
                           start=start(gaPaired[strand(gaPaired) == "+"]),
                           end=mcols(gaPaired[strand(gaPaired) == "+"])$mpos
                           +qwidth(gaPaired[strand(gaPaired) == "+"])))
      tempNeg <- GRanges(seqnames(gaPaired[strand(gaPaired) == "-"]),
                         IRanges(
                           start=mcols(gaPaired[strand(gaPaired) == "-"])$mpos,                        
                           end=end(gaPaired[strand(gaPaired) == "-"])
                         )) 
      temp <- c(tempPos,tempNeg)                
      #temp <- GRanges(seqnames(tempPaired),IRanges(start(left(tempPaired)),end(right(tempPaired))))
      message("..Done.\nRead in ",length(temp)," reads")
      if(!is.null(forceFragment)){
        message("Forcing fragments to be centred and set to ",forceFragment,"..",appendLF=FALSE)
        temp <- resize(temp,forceFragment,"center")
        message("..done")        
      }
      if(!is.null(minFragmentLength)){
        temp <- temp[width(temp) > minFragmentLength]
      }
      if(!is.null(maxFragmentLength)){
        temp <- temp[width(temp) < maxFragmentLength]
      }   
      message("..done")
    }
  
    ## Downsample single or paired end reads if specified and create coverage rlelist
    message("Calculating coverage..",appendLF=FALSE)
    
    seqlengths(temp)[match(names(lengths),names(seqlengths(temp)))] <- lengths
    if(!is.null(downSample)){
        if(downSample < 1 & downSample > 0){
          temp <- temp[sample(length(temp),round(length(temp))*downSample),]
        }else if(downSample > 1){
          temp <- temp[sample(length(temp),downSample),]
        } 
    }

    genomeCov <- coverage(temp)
    lengths <- seqlengths(genomeCov)
    allchrs <- names(lengths)
    message("..done")
  }
  
  chrs_to_use <- intersect(unique(seqnames(testRanges)), unique(names(genomeCov)))
  message("Filtering data to only common chrs: ", paste(chrs_to_use, collapse = ", "))
  testRanges <- testRanges[seqnames(testRanges) %in% chrs_to_use]
  genomeCov <- genomeCov[chrs_to_use]
  
  message("Extracting coverage data...")
  mat <- extract_data_in_ranges(testRanges, genomeCov, method = method, 
                                binMethod = binMethod, nOfWindows = nOfWindows)
  
  message("Processing matrix and creating ChIPprofile object...")
  if(strand.aware == TRUE){
    mat[which(strand(testRanges)=="-"),] <- mat[which(strand(testRanges)=="-"), ncol(mat):1] #flip rows
  }
  
  colnames(mat) <- 1:ncol(mat)
 
  if (length(oob_idx) > 0){
    filteredRanges <- original_ranges[-oob_idx]
    filteredRanges <- filteredRanges[seqnames(filteredRanges) %in% chrs_to_use]
  } else{
    filteredRanges <- original_ranges[seqnames(original_ranges) %in% chrs_to_use]
  }

  profileSample <- SummarizedExperiment(mat, rowRanges = filteredRanges)
  
  if(is.null(samplename)){
    if(is.character(bamFile)){
      samplename <- bamFile # give name of file if bam or bigwig
    } else if (format %in% c("granges", "rlelist")){
      deparse(substitute(bamFile)) #give name of object if granges or rlelist
    } else {
      samplename <- "Sample"
    }
  }
  
  names(assays(profileSample)) <- samplename
  
  paramList <- list(samplename = samplename,
                    #region options
                    style = style,
                    method = method, binMethod = binMethod, nOfWindows = nOfWindows, 
                    expand = expand, flank = flank, flankUp = flankUp, flankDown = flankDown, 
                    #import options
                    format = format, paired= paired, removeDup = removeDup, 
                    FragmentLength = FragmentLength, forceFragment = forceFragment, 
                    minFragmentLength = minFragmentLength, maxFragmentLength = maxFragmentLength,
                    normalize = normalize, downSample = downSample,
                    seqlengths = seqlengths, 
                    #pwm options
                    genome = genome, cutoff = cutoff)
  
  message("Done!")
  return(new("ChIPprofile",profileSample,params=paramList))
}


#' Make ranges for data import

#' Make ranges according to specified style: use regions as-is, expand all to same width,
#' add a constant flank to each side (in bp), or extend by some % to each side.
#'
#' @import GenomicRanges 

make_ranges <- function(testRanges, style = "region", expand = NULL, flank = NULL, 
                        flankUp = NULL, flankDown = NULL){
  
  # Expand ranges as specified by 'style'
  if (style == "region"){
    message("Using testRanges as regions.")
    testRanges <- testRanges
    
  } else if (style == "centered"){
    if (is.null(expand)){
      warning("Range expansion width not specified; using 100bp.")
      expand <- 100
    }
    message(paste("Expanding testRanges to", expand, "bp to use as regions."))
    testRanges <- resize(testRanges, fix = "center", width = expand)
    
  } else if (style == "flanked"){
    
    if (mode(flank) == "character"){ 
      flank_p <- strsplit(flank, "%")[[1]][1] # anything before the % will be taken
      
      message(paste("Expanding testRanges by ", flank_p, "% to use as regions."))
      flank_widths <- width(testRanges) * as.numeric(flank_p) / 100
      start(testRanges) <- start(testRanges) - flank_widths
      end(testRanges) <- end(testRanges) + flank_widths
      
    } else if (mode(flank) == "numeric" | (mode(flankUp) == "numeric" & mode(flankDown) == "numeric")) {
      #either flank or both flankUp andflankDown must be specified
      
      if (is.null(flankUp) & is.null(flankDown)){
        #by default up and down flanks are equal size
        flankUp <- flank
        flankDown <- flank
      }
      
      message(paste("Expanding testRanges by ", flankUp," bp upstream and ", flankDown, "bp downstream to use as regions." ))
      
      if(any(strand(testRanges)=="*")){ stop("strand '*' should be set to '+'")}
      
      start(testRanges) <- start(testRanges) - ifelse(strand(testRanges) == "+", flankUp, flankDown) #correct flank depending on strand
      end(testRanges) <- end(testRanges) + ifelse(strand(testRanges) == "+", flankDown, flankUp)
      
    } else {
      stop("For style 'flanked', 'flank' should be either a character string e.g. '100%', or, if numeric, the number of base pairs to use to flank the regions.")
    }
    
  } else {
    stop(paste0("style '", style,"' not recognised! Must be one of region, centered or flanked."))
  }
  
  return(testRanges)
}

split_ranges <- function(testRanges, flank = NULL, flankUp = NULL, flankDown = NULL){
  #either flank or both flankUp andflankDown must be specified
  if (is.null(flankUp) & is.null(flankDown)){
    #by default up and down flanks are equal size
    flankUp <- flank
    flankDown <- flank
  }
  
  flank_width_up <- flankUp
  flank_width_down <- flankDown
  
  #split into flanking and original ranges
  
  orig_start <- start(testRanges) + ifelse(strand(testRanges) == "+", 
                                           flank_width_up, flank_width_down)
  orig_end <- end(testRanges) - ifelse(strand(testRanges) == "-", 
                                       flank_width_down, flank_width_up)
  
  up_ranges <- testRanges
  end(up_ranges) <- orig_start - 1
  
  down_ranges <- testRanges
  start(down_ranges) <- orig_end +1
  
  orig_ranges <- testRanges
  start(orig_ranges) <- orig_start
  end(orig_ranges) <- orig_end
   
  return(list(up_ranges, orig_ranges, down_ranges))
}

## run this function per chromosome and then combine matrices?
extract_data_in_ranges <- function(testRanges, genomeCov, method = "bin", binMethod = "mean", nOfWindows = 100){
  if (is.null(method)){
    # subset coverage, make vector, return to rbind
    mat <- do.call("rbind", lapply(genomeCov[testRanges], as.vector))
    
  } else if (method == "spline"){
    
    mat <- do.call("rbind", lapply(genomeCov[testRanges], function(x){
      spline(as.vector(x), n = nOfWindows)$y
    }))
    
  } else if (method == "bin"){
    stop("not yet implemented!")
    #     
    #     #calculate bin size including final bin which may be different
    #     window_size <- width(testRanges) %/% nOfWindows
    #     remainder <- width(testRanges) %% nOfWindows
    #     
    #     if any((remainder > window_size)){
    #       warning(paste0("Final bin size (", window_size + remainder,
    #                      " bp) is more than twice the size of other bins (", window_size,
    #                      " bp); consider using 'spline' method instead"))
    #     }
    #     
    #     binstarts <- seq(from = start(testRange), by = window_size, length.out = nOfWindows)
    #     binends <- c(seq(from = start(testRange) + window_size - 1, by = window_size, 
    #                      length.out = nOfWindows - 1), end(testRange))
    #     
    #     bins <- GRanges(seqnames = seqnames(testRange),
    #                     ranges = IRanges(start = binstarts, end = binends))
    #     
    #     bin_cov <- genomeCov[bins]
    #     vec_list <- lapply(unname(bin_cov), as.vector)
    #     
    #     vec <- vapply(vec_list, binMethod, FUN.VALUE = numeric(1))
    #     #binMethod can be any function that returns a single value
    #     #how to get a more useful error message here??
    #     
  } else { 
    stop(paste0("method '", method, 
                "' not recognised! Must be NULL for base-pair resolution,, 'bin', or 'spline'."))
  }
  
  rownames(mat) <- testRanges$original_order
  return(mat)
}


name_cols <- function(mat, style = "region", method = NULL, expand = NULL, flank = NULL, 
                      flankUp = NULL, flankDown = NULL){
  
  if(style == "centered"){
    col_names <- seq(-(expand/2), expand/2, length.out = ncol(mat))
    col_names <- paste(col_names, "bp")
    
  } else if (style == "region"){
    col_names <- seq(-(ncol(mat)/2), ncol(mat/2), length.out = ncol(mat))
    
    if (is.null(method)){
      col_names <- paste(col_names, "bp")
    } else {
      col_names <- paste(col_names, "bins")
    }
  } else if (style == "flank"){
    
    col_names <- seq(-(ncol(mat)/2), ncol(mat/2), length.out = ncol(mat))
    
  }
  
  
}

