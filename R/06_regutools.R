# An Introduction to the R/Bioconductor Ecosystem

# Part 05: regutools package demostration


# First, we need to install the regutools package. Let's use our BiocManager
# function to installe it from Bioconductor.

BiocManager::install("regutools")

# Load regutools
library(regutools)

# The first step every time you use regutools is to connect with the database.
# This step might take a minute.
regulondb_conn <- connect_database()

# Now we need to create a special object to use all the tables contained in the
# database

e_coli_regulondb <-
    regulondb(
        database_conn = regulondb_conn,
        organism = "E.coli",
        database_version = "1",
        genome_version = "1"
    )

#


# Integration of regutools with the Bioconductor Ecosystem


# Let's extract some gene data in the regular way:

res <- get_dataset(
    regulondb = e_coli_regulondb,
    dataset = "GENE",
    attributes = c("posleft", "posright", "strand", "name"),
    filters = list("name" = c("araC","crp","lacI"))
)

# How does the result look?
res

# Its class is regulondb-result, but looks very similar to a data frame
class(res)

res$name

# Note that this table have the colums position, strand and name. I think we
# have seen this before...in the GRanges objects.

# The function convert_to_granges() converts a regulondb_result object into
# a GRanges object to facilitate the integration with other Bioconductor
# workflows.

convert_to_granges(res)

# If we evaluate the class of the result, is a GRanges object now.
GR_res <- convert_to_granges(res)
class(GR_res)

# This integration allows to use your results with other packages, for example
# to visualize their genomic position.
grange <- GenomicRanges::GRanges("chr",IRanges::IRanges(5000, 10000))

res <- get_dna_objects(
    regulondb = e_coli_regulondb,
    grange = grange,
    elements = c("gene", "promoter")
)

# How does the result look?
res

# This does not seem like a data frame or a regulondb_result anymore
class(res)

# Because it is a GRanges object that can be used to plot its elements using
# the dataviz package

e_coli_regulondb <-
    regulondb(
        database_conn = regulondb_conn,
        organism = "chr",
        database_version = "1",
        genome_version = "1"
    )

plot_dna_objects(
    regulondb = e_coli_regulondb,
    grange = grange,
    elements = c("gene", "promoter")
)


