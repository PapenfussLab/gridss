# Script to generate a chart to display GC bias based upon read starts observed
# in windows along the genome.
#
# @author Tim Fennell

# Parse the arguments
args <- commandArgs(trailing=T)
metricsFile  <- args[1]
outputFile   <- args[2]
bamName      <- args[3]

# Figure out where the metrics and the histogram are in the file and parse them out
startFinder <- scan(metricsFile, what="character", sep="\n", quiet=TRUE, blank.lines.skip=FALSE)

firstBlankLine=0

for (i in 1:length(startFinder)) {
  if (startFinder[i] == "") {
    if (firstBlankLine==0) {
      firstBlankLine=i+1
    } else {
      secondBlankLine=i+1
      break
    }
  }
}

histogram <- read.table(metricsFile, header=TRUE, sep="\t", skip=secondBlankLine, comment.char="", quote='', check.names=FALSE)

headers <- sapply(sub(".", "", names(histogram), fixed=TRUE), "[[" , 1)

if (any(duplicated(headers))) {
  print(paste("Not creating insert size PDF as there are duplicated header names:", headers[which(duplicated(headers))]))
} else {
  pdf(pdfFile)

  for (i in 2:ncol(histogram)) {
  	label <- sub(".", "", names(histogram)[i], fixed=TRUE)
  	# pad out missing values
  	df <- data.frame(mapq=0:255, count=rep(0, 256))
  	df$count[histogram$mapq + 1] <- histogram[ ,i]

    plot(x=df$mapq, y=df$count,
         type="s",
    	 col="blue",
         main=paste("Mapq Histogram for", label, "\nin file", bamName),
         xlab="Mapq",
         ylab="Count",
         xlim=range(0, max(histogram$mapq + 1)))
  	#lines(df$mapq, df$count,  type="h", col="red")
  }

  dev.off()
}
