# Using GenomicRanges we can have more information from each sequence, such as
# the chromosome, strand and genomic position

# Load library
library(GenomicRanges)

# Let's see how is the basic structure of a GRanges object


gr <- GRanges(seqnames = c("geneA", "geneB", "geneC"),
              ranges = IRanges(start = c(10, 20, 32),
                               end = c(15,27,42) ),
              strand = c("+", "+", "-")
              )
gr

# Note that we used a special way to define the genomic position (ranges)
IRanges(start = c(10, 20, 22),
        end = c(15,27,32) )

# We can add more information in the sequence name field by using Rle objects

# note that we can specify the range width instead of the end position
IRanges(1:5, width=10:14, names=head(letters, 5))

gr <- GRanges(seqnames = Rle(c("chr2", "chr2", "chr1", "chr3"), c(1, 1, 2, 1)),
              ranges = IRanges(1:5, width=10:14, names=head(letters, 5)),
              strand = Rle(strand(c("-", "+", "*")), c( 1,3, 1)),
              score=1:5, GC=seq(1, 0, length=5) # these are metadata columns
              )
gr

# GenomicRanges has functions to select genomic regions such as introns and exons.
gr <- GRanges(seqnames = "geneM",
              ranges = IRanges(start = c(10, 32),
                               end = c(15,42) ),
              strand = "+",
              exon_id = c(1,2) )
gr

# we can know the pre-processed gene range, including introns and exons
range(gr)

# Using flank() we can retrieve the upstream regulatory region

gr <- GRanges(seqnames = "geneO",
              ranges = IRanges(start = c(120, 180),
                               end = c(150,240) ),
              strand = "+",
              exon_id = c(1,2) )
# we indicate the number of upstream nucleotides to consider as the regulatory region
flank(gr, 100)

