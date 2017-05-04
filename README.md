[Master thesis] 
# Reproducable methods for network analysis of high-throughpot genomic data
#### A reimplementation of a Regulatory Network for Th17 Cell Specification
[Source for reference method on Cell](http://www.cell.com/cell/abstract/S0092-8674(12)01123-3)

This project is part of my master thesis at the Center for Bionformatics, University of Hamburg, and the University Medical Center (UKE) in Hamburg.

Given the information provided on the Cell page as well as [available data on NCBI GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE40918), I am building a similar version which adapts methods of the described pipeline.

The pipeline follows the 'Computational Methods' described in Ciofani et al. [1] and adaptations in response to examination of source code provided by Aviv Madar. 

[Project Trello board](https://trello.com/b/E6WcAF7I/th17-regulatory-networks)

## Requirements
- Git
- R
- Bash (for setup.sh)
- Disk space for downloaded GEO data (~ 150MB)
- [Cytoscape](https://www.cytoscape.org) (+ [AllegroLayout](http://allegroviva.com/allegrolayout2/) plugin)
- Recommended: LibreOffice available on command-line so xlsx/xls-files can be automatically converted to CSV.

## Getting ready
Clone the git repository into your preferred file directory by running the following command.

`git clone https://github.com/ckeil0689/master-thesis.git`

From the directory master-thesis/ run the setup script.

`./setup.sh`

Follow instructions on screen. It is important to convert mmc4.xlsx and mmc5.xls to CSV-format to get comma separated tables. If LibreOffice is installed and functional on the command-line, the script will take care of it. Otherwise both files have to be manually converted.

## Creating a network
When setup is complete and all NCBI GEO data has been downloaded, change to the /scripts/ directory to run the master script.

`./master.R`

You should now have a file with network interactions (e.g. `kc_single_1.65_cs-cut_<date>.csv`). The columns are `nodeA | interaction | nodeB | confidence_score`

If you o not want to load ChIP-seq and KO-files again and reuse existing matrices (after at least one run), you can use the `noload` flag.

`./master.R noload`

## Visualizing the network
Open Cytoscape. You can load the network file (e.g. `kc_single_1.65_cs-cut_<date>.csv`) as a network table. Z-scores for genes have been extracted from NCBI GEO GSE40918 and are loaded into a separate file `zscores.txt`. Load this as node table to achieve node coloring by differential expression. Import the XML-style file from /scripts/cyt_styles to achieve the same style as used in [1]. 

Install and run the AllegroLayout plugin in Cytoscape. For example, use the Fruchterman-Reingold algorithm with randomized starting positions.


[1] Ciofani, Maria, Aviv Madar, Carolina Galan, MacLean Sellars, Kieran Mace, Florencia Pauli, Ashish Agarwal et al. "A validated regulatory network for Th17 cell specification." Cell 151, no. 2 (2012): 289-303.
