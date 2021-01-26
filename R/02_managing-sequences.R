# Manipulation of genomic sequences

# The package **Biostrings** have functions to read, write and handling genomic sequences.

# Load library
library(Biostrings)

## Import data with readDNAStringSet(), usually in fasta format.
eco <- readDNAStringSet("../data/eco.fasta")
eco

# We can get some attributes of the data

# Get the number of sequences
length(eco)

# Get the number of characters on each sequence.
nchar(eco)

# Get the frequency of specific characters
letterFrequency(eco, "GC")

# Also we can Subset sequences

# Get the first two sequences from the eco dataset
eco[1:2]

# Get the first 10 characters from each sequence in the eco dataset
subseq(eco,start = 1, end = 10)

# Bioostriings als have functions to:
#
# Translate sequences to amin acids
translate(eco$`eco-b0001`)

# Get the reverse sequence
reverse(eco)

# Get the complement sequence
complement(eco)

# Get the reverse complement sequence
reverseComplement(eco)
