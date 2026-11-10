
# Application specification: CodonTree

This is the application specification for service with identifier CodonTree.

The backend script implementing the application is [App-CodonTree.pl](../service-scripts/App-CodonTree.pl).

The raw JSON file for this specification is [CodonTree.json](CodonTree.json).

This service performs the following task:   Computes a phylogenetic tree based on protein and DNA sequences of PGFams for a set of genomes

It takes the following parameters:

| id | label | type | required | default value |
| -- | ----- | ---- | :------: | ------------ |
| output_path | Output Folder | folder  | :heavy_check_mark: |  |
| output_file | File Basename | wsid  | :heavy_check_mark: |  |
| genome_ids | Main genomes | list  |  | ARRAY(0x560bf4ac3948) |
| genome_groups | Main genomes | list  |  | ARRAY(0x560bf4b56438) |
| optional_genome_ids | Optional genomes (not penalized for missing/duplicated genes) | list  |  | ARRAY(0x560bf4ac3cd8) |
| genome_metadata_fields | Genome Metadata Fields | string  |  |  |
| number_of_genes | Desired number of genes | int  |  | 20 |
| bootstraps | Number of bootstrap replicates | int  |  | 100 |
| sequence_type | Sequences to analyze (both, dna, protein) | enum  |  | both |
| max_genomes_missing | Number of main genomes allowed missing from any PGFam | int  |  | 0 |
| max_allowed_dups | Number of duplications allowed for main genomes in any PGFam | int  |  | 0 |

