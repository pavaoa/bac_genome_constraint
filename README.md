# Constraints on Bacterial Genome Size  

HST.508 Final Project (Fall 2025)  

**Authors:** Aidan Pavao, Terry Cho  



## Project Description  

Bacterial genomes vary nearly an order of magnitude in size. Interestingly, bacteria from similar environments across diverse taxa tend to have similar genome sizes (e.g., human skin commensals ~2.5 Mb). It has been proposed that this convergence reflects metabolic complexity of the environment, while an alternative hypothesis suggests that genome size is constrained by nutrient limitation (e.g., nitrogen conservation).  



The goal of this project is to model how environmental factors constrain bacterial genome size. We will analyze genome size, metabolic coding capacity (via KEGG), and other genomic signatures of nutrient limitation (e.g., amino acid usage, pathway composition) across multiple well-defined environments.



## Data  

Filtered reference genomes were obtained from NCBI:  

[NCBI Reference Genome Dataset](https://www.ncbi.nlm.nih.gov/datasets/genome/?taxon=2&reference_only=true&annotated_only=true&refseq_annotation=true&genbank_annotation=true&typical_only=true&exclude_mags=true&exclude_multi_isolates=true&assembly_level=3:3)  



Criteria:  

- Reference genomes only  

- Complete assembly level  

- Released since 2010  

- Annotated (RefSeq and GenBank)  

- Excluded atypical, MAGs, and multi-isolate assemblies  



**Total genomes:** 5,988  

(`0_data_ncbi_bac_table.tsv`)



## Preliminary Analysis  

Initial exploration includes:  

- Histogram and scatter plots of genome size distribution  

- Assessment of genome size variation across environments  

- Verification that data are appropriate for downstream KEGG and nutrient limitation analyses  



## Planned Analysis  

- KEGG annotation of protein sequences to estimate metabolic capacity  

- Correlation of genome size with environment, metabolic coding capacity, and nutrient limitation signatures  

- Model building (linear regression or machine learning) incorporating phylogeny and environment  



## Contributions  

- **Aidan Pavao:** Primary data analysis, visualization, modeling  

- **Terry Cho:** Data setup, repository organization, documentation  



## Repository Structure  



data/                # Input datasets (NCBI genome tables, KEGG mappings)

scripts/             # Analysis and plotting scripts

results/             # Output figures and summary tables

README.md            # Project overview



## Timeline  

- **Data criteria due:** Nov 7, 2025  

- **Preliminary analysis due:** Nov 17, 2025  

- **Presentation:** Dec 8–10, 2025  

- **Final paper due:** Dec 12, 2025  

