# AnnotationHub connect with multiple genomic data providers and give us access to their
# data

# Load library
library(AnnotationHub)

# Connect with the AnnotationHub database
ah = AnnotationHub()
ah

# We can see which species are available
unique(ah$species) %>% head(10)
unique(ah$species) %>% length()

# Also we can search by genome name and provider
genomes <- unique(ah$genome)
genomes[grep("hg",genomes)]

unique(ah$dataprovider)

# We can also see the class of data contained
unique(ah$rdataclass)

# query() helps to retrieve the information from a specie of interest
hg <- query(ah, c("UCSC", "Homo sapiens", "GRanges") )

# We still have a lot of information in hg, let's use display() to facilitate
# exploration of the data
display(hg)

# Once we know the data we want to retrieve, we can select the data with [ ]
hg_genes <- hg[["AH5036"]]

# Note that the result is a GRanges object

################################################################################

# There are some packages that give us access to databases from a specific
# organism, such as **org.Hs.eg.db** that maps Gene identifiers to GenBank
# Accession Numbers for the Human genome.

#  Load the libraries
library(AnnotationDbi)
library(org.Hs.eg.db)

# How does the object look?
org.Hs.eg.db

# Which information is available in this database?
x <- org.Hs.eg.db

# Retrieve the available data sets names
keytypes(x)

# Then, let's see the available names for the desired data set
keys(x, keytype = "GENENAME")

# What kind of information can we retrieve?
columns(x)

# Now let's access to some information
#
# As and example, let's see the information for the "ATRX chromatin remodeler"
AnnotationDbi::select(x,
                      keys = "ATRX chromatin remodeler",
                      keytype = "GENENAME",
                      columns = c("ENTREZID", "ALIAS", "UNIPROT")
                      )
