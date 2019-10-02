write("Loading dependencies...", stdout())
suppressMessages(suppressWarnings(library("optparse")))
suppressMessages(suppressWarnings(library("minfi")))
suppressMessages(suppressWarnings(library("conumee")))
suppressMessages(suppressWarnings(library("stringr")))
suppressMessages(suppressWarnings(library("utils")))
suppressMessages(suppressWarnings(library("rtracklayer")))
suppressMessages(suppressWarnings(library("foreach")))
suppressMessages(suppressWarnings(library("doParallel")))
suppressMessages(suppressWarnings(library("parallel")))
suppressMessages(suppressWarnings(library("biomaRt")))

# Parse input arguments
parser = OptionParser()

parser <- add_option(parser, c("--data"), help = "zip or gzip containing DNA methylation microarray data")
parser <- add_option(parser, c("--controlsdata"), help = "controls methylation data")
parser <- add_option(parser, c("--controls"), help = "names of control samples")
parser <- add_option(parser, c("--genesfile"), help = "file with list of genes to highlight")
parser <- add_option(parser, c("--ignorefile"), help = "bed file of regions to exclude")
parser <- add_option(parser, c("--xy"), help = "include XY chromosomes")
args <- parse_args(parser)

preprocess.minfi <- function(rawdata) {
    # Path to folder containing raw experiment .IDAT files. Will read Illumina format SampleSheet.csv to load data if found.
    if (endsWith(rawdata, ".tar.gz")) {
    untar(rawdata, exdir = "rawdata")
    data.paths <- paste("rawdata", untar(rawdata, exdir = "rawdata", list = TRUE), sep = "/")
    } else if (endsWith(rawdata, ".gz")) {
    gunzip(rawdata, exdir = "rawdata")
    data.paths <- paste("rawdata", gunzip(rawdata, exdir = "rawdata", list = TRUE), sep = "/")
    } else if (endsWith(rawdata, ".zip")) {
    unzip(rawdata, exdir = "rawdata")
    data.paths <- paste("rawdata", unzip(rawdata, exdir = "rawdata", list = TRUE), sep = "/")
    } else {
    write("Only tar.gz, .gz, and .zip archive file extensions are supported for sample data input.", stdout())
    stop()
    }

    # Read one directory down if top of archive is a single folder
    if (length(data.paths) == 1) {
    data.folder <- paste(getwd(), "rawdata", strsplit(data.paths[[1]], "/")[[1]], sep = "/")
    } else {
    data.folder <- paste(getwd(), "rawdata", sep = "/")
    }

    write("Loading IDAT files...", stdout())
    # Try Load sample sheet
    if (any(grepl(".csv$", list.files(data.folder)))) {
    targets <- read.metharray.sheet(data.folder)

    # remove unannotated samples
    targets.rmdups <- targets[targets$Basename != "character(0)", ]
    experiment.rgset <- read.metharray.exp(targets = targets.rmdups)

    # Recursively find all .idat files in data.folder if sample sheet not found
    } else {
    experiment.rgset <- read.metharray.exp(base = data.folder, recursive = TRUE)
    }

    # Perform preprocessing
    write("Perform background subtraction and control normalization as implemented by Illumina's GenomeStudio.", stdout())
    methyl.set <- preprocessIllumina(experiment.rgset)

    write("Cleaning up intermediate files...", stdout())
    unlink("rawdata", recursive = TRUE)
    write("Done.", stdout())

    return(methyl.set)
}

queryBiomart <- function(args){
	# Query biomart for gene annotations
	if (!("genesfile" %in% names(args))) {
	  detail_regions <- NULL
	} else {
	  write("Querying biomaRt for gene locations...", stdout())
	  gene.list <- read.table(args$genesfile, header = FALSE)
	  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "grch37.ensembl.org")
	  detail_regions <- getBM(attributes = c("hgnc_symbol", "chromosome_name", "start_position",
		"end_position", "strand"), filter = "hgnc_symbol", values = gene.list, mart = ensembl)

	  valid_regions <- detail_regions$chromosome_name %in% seq(1, 22)
	  detail_regions <- detail_regions[valid_regions, ]

	  # Include 5kb upstream for promoter region
	  for (i in 1:dim(detail_regions)[1]) {
		entry <- detail_regions[i, ]
		if (entry$strand == 1) {
		  entry$start_position <- entry$start_position - 5000
		} else {
		  entry$end_position <- entry$end_position + 5000
		}
	  }

	  # Convert to GRange object
	  detail_regions$strand[detail_regions$strand == 1] <- "+"
	  detail_regions$strand[detail_regions$strand == -1] <- "-"
	  detail_regions$chromosome_name <- paste("chr", detail_regions$chromosome_name,
		sep = "")
	  detail_regions <- makeGRangesFromDataFrame(detail_regions, seqnames = "chromosome_name",
		start.field = "start_position", end.field = "end_position", strand.field = "strand",
		keep.extra.columns = TRUE)
	  detail_regions$name <- detail_regions$hgnc_symbol
	  detail_regions$hgnc_symbol <- NULL
	  detail_regions <- trim(detail_regions)
	}
	return(detail_regions)
}

setGeneric("CNV.fit", function(query, ref, anno, ...) {
    standardGeneric("CNV.fit")
})

# CNV.fit from conumee with controls n=1 bugfix
#' @rdname CNV.fit
setMethod("CNV.fit", signature(query = "CNV.data", ref = "CNV.data", anno = "CNV.anno"), 
    function(query, ref, anno, name = NULL, intercept = TRUE) {
        if (ncol(query@intensity) == 0) 
            stop("query intensities unavailable, run CNV.load")
        if (ncol(ref@intensity) == 0) 
            stop("reference set intensities unavailable, run CNV.load")
        
        if (ncol(query@intensity) != 1) 
            stop("query contains more than one sample.")
        if (ncol(ref@intensity) == 1) 
            warning("reference set contains only a single sample. use more samples for better results.")
        
        p <- names(anno@probes)  # ordered by location
        if (!all(is.element(p, rownames(query@intensity)))) 
            stop("query intensities not given for all probes.")
        if (!all(is.element(p, rownames(ref@intensity)))) 
            stop("reference set intensities not given for all probes.")
        
        object <- new("CNV.analysis")
        object@date <- date()
        object@fit$args <- list(intercept = intercept)
        
        if (!is.null(name)) {
            names(object) <- name
        } else {
            names(object) <- colnames(query@intensity)
        }
        object@anno <- anno
        
        r <- cor(query@intensity[p, ], ref@intensity[p, ])
        if (!is.matrix(r)){ # cor returns numeric if 1x1, matrix otherwise...
            r <- as.matrix(r)
            colnames(r) <- colnames(ref@intensity)
        }
        r <- r[1, ] < 0.99
        if (any(!r)) message("query sample seems to also be in the reference set. not used for fit.")
        if (intercept) {
            ref.fit <- lm(y ~ ., data = data.frame(y = query@intensity[p, 
                1], X = ref@intensity[p, r]))
        } else {
            ref.fit <- lm(y ~ . - 1, data = data.frame(y = query@intensity[p, 
                1], X = ref@intensity[p, r]))
        }
        object@fit$coef <- ref.fit$coefficients
        
        ref.predict <- predict(ref.fit)
        ref.predict[ref.predict < 1] <- 1
        
        object@fit$ratio <- log2(query@intensity[p, 1]/ref.predict[p])
        object@fit$noise <- sqrt(mean((object@fit$ratio[-1] - object@fit$ratio[-length(object@fit$ratio)])^2, 
            na.rm = TRUE))
        
        return(object)
    })

###########
# Script ##
###########

# Parse args

# Include/exclude sex chromosomes
if (tolower(args$xy) == "yes") {
  xy <- TRUE
} else {
  xy <- FALSE
}

# Ignorefile (highly polymorphic regions of genome)
if (!("ignorefile" %in% names(args))) {
  ignore_regions <- NULL
} else if (endsWith(args$ignorefile, ".bed")) {
  write("Loading regions to exclude...", stdout())
  ignore_regions <- args$ignorefile
} else {
  write(paste("Ignore regions value, '", args$ignorefile, "' not recognized.",
    sep = ""), stdout())
  stop()
}

# Genesfile
detail_regions <- queryBiomart(args)


# Load datasets
samples.data <- preprocess.minfi(args$data)

# Identify control samples, first check name, then check for separate control dataset
if (tolower(args$controls) != "none") {
  if (!all(controls.names %in% colnames(samples.data))) {
    write("One or more of the provided names of control samples are not in the dataset.",
      stdout())
    stop()
  }
  
  write("Using subset as controls...", stdout())
  controls.names <- strsplit(args$controls, ",")
  controls.names <- lapply(controls.names, trimws)
  controls.names <- unlist(controls.names)

  arraytype <- str_sub(samples.data@annotation['array'], -4)
  controls.data <- samples.data[controls.names]

} else {
  if (!("controlsdata" %in% names(args))) {
    write("Please specify control sample names or provide a separate controls methylation data file.", stdout())
    stop()
  }

  controls.data <- preprocess.minfi(args$controlsdata)
  controls.names <- colnames(controls.data)

  # Auto-detect array type
  # arraytype = input to conumee create_anno.
  if (samples.data@annotation['array'] == controls.data@annotation['array']) {
    all.data <- combineArrays(samples.data, controls.data, outType = samples.data@annotation['array'])
    arraytype <- str_sub(samples.data@annotation['array'], -4)
  } else {
    all.data <- combineArrays(samples.data, controls.data, outType = "IlluminaHumanMethylation450k")
    arraytype <- "overlap"
  }
}


# ---------- QC plots ----------
pdf("qcPlots.pdf")

# Plot median log2 meth vs unmeth intensities
qc <- getQC(all.data)
plotQC(qc)

# Plot Beta value distributions
phenoData <- pData(all.data)
phenoData$isControl <- colnames(all.data) %in% controls.names
phenoData$isControl <- ifelse(phenoData$isControl == TRUE, "Control", "Query")

densityPlot(all.data, sampGroups = phenoData$isControl)
dev.off()

# ----------- conumee analysis -----------

# Create annotation object. Set to look at probes common to both 850k and 450k.
write("Creating CNV annotation object...", stdout())
anno <- CNV.create_anno(array_type = arraytype, exclude_regions = ignore_regions,
  detail_regions = detail_regions, chrXY = xy, bin_minprobes = 15, bin_minsize = 50000,
  bin_maxsize = 5e+06)
if (arraytype == 'overlap'){
	anno@probes <- subset(anno@probes, names(anno@probes) %in% rownames(samples.data))
	anno@probes <- subset(anno@probes, names(anno@probes) %in% rownames(controls.data))
}

# Create CNV object from methylation data
write("Creating CNV object...", stdout())
samples.data <- CNV.load(samples.data)
controls.data <- CNV.load(controls.data)

# parallelization parameters
all_samples <- colnames(samples.data@intensity)
# numCores <- min(max(1, length(all_samples) - length(controls.names)), detectCores())
# cl <- makeCluster(numCores)
# registerDoParallel(cl)

# CNV analysis/plot function
cnv.analyze.plot <- function(sample, controls.data, sample.data, anno) {
  write(paste(sample, " - fitting", sep = ""), stdout())

  cnv.analysis <- CNV.fit(sample.data[sample], controls.data, anno)
  cnv.analysis <- CNV.bin(cnv.analysis)
  cnv.analysis <- CNV.detail(cnv.analysis)
  cnv.analysis <- CNV.segment(cnv.analysis)

  # Open file to write
  write(paste(sample, " - plotting", sep = ""), stdout())
  plot.filename <- paste(sample, ".cnvPlots.pdf", sep = "")

  pdf(plot.filename, width = 16, height = 8)

  # Overview plots
  CNV.genomeplot(cnv.analysis)
  CNV.detailplot_wrap(cnv.analysis)

  # Chromosome plots
  for (chrom in rownames(anno@genome)) {
    CNV.genomeplot(cnv.analysis, chrom)
  }

  # Detail plots
  genes <- names(cnv.analysis@detail$probes)
  for (gene in genes) {
    CNV.detailplot(cnv.analysis, gene)
  }
  dev.off()

  write(paste(sample, " - writing seg files", sep = ""), stdout())
  # Export seg file
  CNV.write(cnv.analysis, what = "segments", file = paste(sample, ".cnv.seg", sep = ""))

  # Export details file
  CNV.write(cnv.analysis, what = "detail", file = paste(sample, ".detail.cnv.seg", sep = ""))
}

# Analyze in parallel
write(paste("Using ", 1, " cores...", sep = ""), stdout())
for (sample in all_samples) {
    cnv.analyze.plot(sample, controls.data, samples.data, anno)
}

write("Done.", stdout())




