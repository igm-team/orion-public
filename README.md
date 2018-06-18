# Orion: Hunting The Intolerant Regions of the Genome
This repository contains the code used for the Orion methodology (http://bit.ly/2hOI39X).

## Citation
Gussow AB, Copeland BR, Dhindsa RS, Wang Q, Petrovski S, Majoros WH, Allen AS, Goldstein DB. Orion: Detecting regions of the human non-coding genome that are intolerant to variation using population genetics. PLoS ONE. 2017; 12(8): e0181604.

## Preparation To Run Orion
We recommend using [conda](https://conda.io/docs/) to create an environment in which to run Orion. Orion requires Python 2.7.

Orion relies on the [luigi](https://github.com/spotify/luigi) framework, with a couple of in-house modifications. Specifically, we added two new 
parameter types, named InputFileParameter and OutputFileParameter.

After installing the luigi framework, you'll need to add these new parameter types to your luigi installation. Append the contents of src/luigi\_parameter\_extension.py to the existing parameter.py file within the installed luigi package (luigi/parameter.py):

>cat src/luigi\_parameter\_extension.py >> path/to/luigi/parameter.py

Following this, you'll need to edit the \_\_init\_\_.py file in the installed luigi package so that it loades the new parameter types.

In your installation of luigi, in luigi/\_\_init\_\_.py edit line that begins with:
>from luigi.parameter import [...here you'll see a list of luigi parameter types...]

Add the new parameter types (InputFileParameter, OutputFileParameter) at the end of the import list on this line.

## Calculating Orion Scores And Regions
How to calculate Orion scores:  
>src/Orion.py CalculateOrionScores \  
> --sample-size \<sample_size\> \  
> --coverage-file \<coverage_file\> \  
>	--allele-counts-file \<allele_counts_file\> \  
>	--window-length \<window_length\> \  
>	--output-path \<output_path:final output file\> \  
>	--output-directory \<output_directory:temporary output files directory\> \  
>	--workers \<num_workers>

Orion regions can be calculated using getMCFRegions.R. For more information, run:
>Rscript src/getMCFRegions.R --help

## Data Generated Using The Orion Methodology
Throughout the study we generated three datasets:
1. [Orion Scores](https://doi.org/10.6084/m9.figshare.4541632.v1)
2. [Orion regions](https://doi.org/10.6084/m9.figshare.4536101.v1)
3. [Coordinates](https://doi.org/10.6084/m9.figshare.4536095.v1) of defined Orion scores, non-repeat autosomal regions that were covered in our sample
