# Plotting biological data with Gviz

# load the package
library(Gviz)

# Let's get some data. Genomic features may either be extracted from public data
# bases like ENSEMBL or UCSC, or they may be generated or curated in-house.

# Gviz validates that the chromosome names provided match with existent
# UCSC Chromosome Names. All of them must start with 'chr'. You may turn off
# the validation with options(ucscChromosomeNames=FALSE)

# A sample set of gene coordinates has been saved in the cpgIslands data set.

data(cpgIslands)
cpgIslands

# note that this data set is a GRanges object
class(cpgIslands)

# now we need to transform the GRanges object to an AnnotationTrack object
atrack <- AnnotationTrack(cpgIslands, name = "CpG")
atrack

# and plot the features available for CpG
plotTracks(atrack)

# Let's add a genomic axis as position reference
gtrack <- GenomeAxisTrack()

plotTracks(list(gtrack, atrack))

# Let's add a chromosome ideogram, which is a simplified visual representation
# of a chromosome
#
# first, we need to identify the genome corresponding to our data
gen <- genome(cpgIslands)

# and store the sequence names of the features
chr <- as.character(unique(seqnames(cpgIslands)))

# now we are ready to annotate the features
itrack <- IdeogramTrack(genome = gen, chromosome = chr)

# Alright! Now we can plot all the elements. Note the red mark indicating the
# position of the chromosome amplified when plotting the features
plotTracks(list(itrack, gtrack, atrack))

# We can add more tracks, such  as the Gene models
data(geneModels)
grtrack <- GeneRegionTrack(geneModels, genome = gen,
                           chromosome = chr, name = "Gene Model")
plotTracks(list(itrack, gtrack, atrack, grtrack))

# If we need to focus on an specific region, we can select a range
plotTracks(list(itrack, gtrack, atrack, grtrack),
           from = 26700000, to = 26750000)

# Plot customization
# We can change the display parameters directly. Note that we are adding gene symbols
# and changing the title background color
grtrack <- GeneRegionTrack(geneModels, genome = gen, chromosome = chr,
                           name = "Gene Model",
                           transcriptAnnotation = "symbol",
                           background.title = "brown")

# Or, we can see the display parameters
head(displayPars(grtrack))

# and change what we want to customize
displayPars(grtrack) <- list(background.panel = "#FFFEDB", col = NULL)

# Now the plot looks much better
plotTracks(list(itrack, gtrack, atrack, grtrack))
