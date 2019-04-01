# Synteny #

## clade.py ##
Usage clade.py [--debug] [--caller_list] --clades CLADES_FILE --clusters GENE_CLUSTER_FILE


## CLADES_FILE ##
A tab-separated file with one header line. Headers must include "genome" and "clade". 

## GENE_CLUSTER_FILE ##
A tab-separated file with one header line. Headers must include "gene_cluster_id", "genome_name", and "gene_callers_id". 
Gene Callers Id rows are expected to be unique integers, and in order based on genomic position. 
