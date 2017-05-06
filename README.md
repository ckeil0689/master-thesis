[Master thesis] 
# Reproducible methods for network analysis of high-throughput genomic data
#### A reimplementation of a Regulatory Network for Th17 Cell Specification ([Source for reference method on Cell](http://www.cell.com/cell/abstract/S0092-8674(12)01123-3))

This project is part of my master thesis at the Center for Bionformatics, University of Hamburg, and the University Medical Center (UKE) in Hamburg.

Given the information provided on the Cell page as well as [available data on NCBI GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE40918), I am building a similar version which adapts methods of the described pipeline.

The pipeline follows the 'Computational Methods' described in Ciofani et al. [1] and adaptations in response to examination of source code provided by Aviv Madar. 

[Project Trello board](https://trello.com/b/E6WcAF7I/th17-regulatory-networks)

---
### Present Status
- Existing data made available by the authors of [1] is used to create a interaction network between transcription factors and target genes.
- Currently, ChIP-seq and DESeq data is used to calculate the interaction network (2 data types).

### Working On
- Use your own custom DESeq data to create an interaction network.
- Use your own custom RNA-seq results to display differential expression as node colors (calculate z-scores instead of importing them from mmc5.xls)

---
## Requirements
- GNU/Linux or macOS
- Git
- R
- Bash (for setup.sh)
- Disk space for downloaded GEO data (~ 150MB)
- [Cytoscape](https://www.cytoscape.org) (+ [AllegroLayout](http://allegroviva.com/allegrolayout2/) plugin)
- Recommended: LibreOffice available on command-line so xlsx/xls-files can be automatically converted to CSV.

## Getting ready
Clone the git repository into your preferred directory by running the following command.

`git clone https://github.com/ckeil0689/master-thesis.git`

From the directory `/master-thesis/` run the setup script.

`./setup.sh`

Follow instructions on screen. It is important to **convert mmc4.xlsx and mmc5.xls to CSV-format** to get comma separated tables. If LibreOffice is installed and functional on the command-line, the script will take care of it. Otherwise both files have to be manually converted.

## Testing Setup
In order to test whether setup was completed as expected, run the main test script. From `/scripts/` change to the `/test/` directory and run 

`./run-tests.R`

If the run completes without errors, you can proceed to the next step. If errors occur make sure you

1) Have an internet connection so all required libraries can be loaded from cran. 
2) Converted mmc4.xlsx and mmc5.xls to comma-separated CSV files.

Otherwise, please screenshot or copy the error output and post it in the Issues subsection.

## Creating a Network
When setup is complete and all NCBI GEO data has been downloaded, change to the `/scripts/` directory to run the master script.

`./master.R`

You should now have a file with network interactions located at `/suppl/data/analysis/cyt/` (e.g. `kc_single_1.65_cs-cut_<date>.csv`). The columns are `nodeA | interaction | nodeB | confidence_score`

If you do not want to load ChIP-seq and KO-files again and reuse existing matrices (after at least one run), you can use the `noload` flag.

`./master.R noload`

## Visualizing the Network
Open Cytoscape. You can load the network file (e.g. `kc_single_1.65_cs-cut_<date>.csv`) as a network table. Z-scores for genes have been extracted from NCBI GEO GSE40918 and are loaded into a separate file `zscores.txt` (also located at `/suppl/data/analysis/cyt/`). Load this as node table to achieve node coloring by differential expression. Import the XML-style file from /scripts/cyt_styles to achieve the same style as used in [1]. 

Install and run the AllegroLayout plugin in Cytoscape. For example, use the Fruchterman-Reingold algorithm with randomized starting positions.


*[1] Ciofani, Maria, Aviv Madar, Carolina Galan, MacLean Sellars, Kieran Mace, Florencia Pauli, Ashish Agarwal et al. "A validated regulatory network for Th17 cell specification." Cell 151, no. 2 (2012): 289-303.*
