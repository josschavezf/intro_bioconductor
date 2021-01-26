# Let's install all packages needed for today's session.

# BiocManager
install.packages("BiocManager")

# Biostrings
BiocManager::install("Biostrings")

# GenomicRanges
BiocManager::install("GenomicRanges")

# AnnotationHub
BiocManager::install("AnnotationHub")

# org.Hs.eg.db
BiocManager::install("org.Hs.eg.db")

# Gviz
BiocManager::install("Gviz")

# SummarizedExperiment
BiocManager::install("SummarizedExperiment")

# airway
BiocManager::install("airway")

# regutols
BiocManager::install("regutools")
