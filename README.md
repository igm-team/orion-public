# Orion: Hunting The Intolerant Regions of the Genome
This repository contains the code used for the Orion methodology.

How to calculate Orion scores:  
>src/Orion.py CalculateOrionScores \  
> --sample-size \<sample_size\> \  
> --coverage-file \<coverage_file\> \  
>	--allele-counts-file \<allele_counts_file\> \  
>	--window-length \<window_length\> \  
>	--output-path \<output_path:final output file\> \  
>	--output-directory \<output_directory:temporary output files directory\> \  
>	--workers \<num_workers>

Orion regions can be calculated using getMCFRegions.R. Run:
>Rscript src/getMCFRegions.R --help

For more information.
