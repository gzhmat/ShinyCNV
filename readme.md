ShinyCNV
================

A Shiny/R application to view and annotate copy number variations
-----------------------------------------------------------------

ShinyCNV is developed by wrapping up the graphic and data-table processing functions in R package, and the interactive features are implemented by Shiny package. Users can visually check normalized SNP data (either from Illumina or Affymetrix platform) together with reported CNVs from any CNV detection tools, and semi-atomically edit and update the CNVs. Detailed steps are listed below and a video tutorial is available at YouTube.

Installation
------------

1.  Install [R environment](https://www.r-project.org/)
2.  Install [RStudio](https://www.rstudio.com/).
3.  Install ShinyCNV
    -   Download ShinyCNV released version (stable) or testing version (from master branch).
        The full version includes 11 samples' SNP data (~350Mb), and the lite version includes 2 paired samples' SNP data (~60Mb).
    -   Unzip the ShinyCNV package and open either "ui.R" or "server.R" in RStudio
    -   Click "Run App" button at the top-right corner, the availability of relied R packages will be automatically checked and installed.
    -   After successful installation, you should be able to see this without error:

    ![first view](./readme_files/fig/1.PNG "fig:")

Input files
-----------

1.  Normalized SNP data
    -   Same input file for [OncoSNP](https://sites.google.com/site/oncosnp/user-guide/input-files). The following 5 columns are required:
        -   SNP Name
        -   Chr
        -   Position
        -   Log R Ratio
        -   B Allele Freq
    -   Detailed steps on preparing [**Illumina**](http://penncnv.openbioinformatics.org/en/latest/user-guide/input/) and [**Affymetrix**](http://penncnv.openbioinformatics.org/en/latest/user-guide/affy/) data are available from PennCNV's website.
2.  Reported CNV regions
    -   The reported CNVs could be called by software, as long as the following columns are provided (use the exact same column names):
        -   index `/CNV id asigned by users, which is used to compare the updated CNVs with the original ones`
        -   caseID `/case ID; SNP data file is caseID.{suffix}`
        -   controlID `/control ID; SNP data file is control.{suffix}; for unpaired cases, use other control`
        -   chr `/chromosome`
        -   start `/CNV start position`
        -   end `/CNV end position`
        -   CN `/copy number`
        -   normRate `/normal sample contamination rate; set 0 if unknown`
        -   gender `/Female or Male; set 'Unknown' if unknown`
        -   paired `/Yes or No`

Annotate CNV table
------------------

Now we are ready to go! Click the "Browse" button to import the CNV region file:

![CNV table panel](./readme_files/fig/2.PNG)

Within this panel, you could:

-   select a CNV segment by clicking
-   read in SNP data based on the imported case ID and control ID
-   use "Prev/Next CNV" button to navigate the selected CNV
-   add/delete CNV
-   mark selected CNV as germline/covered/false/true by clicking the buttons
-   set chr/copy number(CN)/normRate(NR); "set chr" is only for newly added CNVs
-   clear the imported SNP data, which usually take the longest time, so be cautious on doing this

Check CNVs and update their breakpoints
---------------------------------------

The most useful function of this app is to manually check each CNV and adjust inaccurate breakpoints, which is in the LRR/BAF panel as below:

![LRR/BAF panel](./readme_files/fig/3.PNG)

-   This figure will be shown once CNV "A004" is selected. Each dot represents a SNP probe, with X axis along chromosome coordination and Y axis showing normalized LRR and BAF.
-   Two red vertical lines indicate the reported breakpoints of selected CNV, which is obviously correct according to the figure.
-   The RefSeq genes are shown as bars at the bottom panel, and the COSMIC cancer genes are marked in red. Detailed gene information could be checked by clicking. E.g., the red bar is clicked and the information of gene *ETV6* is shown at the bottom.
-   To change the breakpoints of selected CNV, users could zoom in LRR/BAF graph by mouse-swipe and then pick the correct position by clicking. Chromosome position of the SNP nearest to the click will be shown in box "Pos:" and the start/end positions could be updated by "Set Start/End" buttons.
-   In case of marking whole chromosome gain/loss, users need the start and end positions of that chromosome. To do this, the dropdown list "Chr:" is very handy. Usually the start of each chromosome is 0, but for chromosome 13, 14, 15, 21 and 22, the P arm is not assessable and thus would start from the first SNP probe according to the most recent Infinium Omni2.5Exome-8 array.

CNV spectrum for imported samples
---------------------------------

After loading SNP data , LRR across genome will be shown in the spectrum panel below:

![Genome spectrum panel](./readme_files/fig/4.PNG)

-   Blue means LRR is below 0 (copy number loss) while red means above 0 (copy number gain).
-   Case IDs are on the left side and chromosomes are on the top.
-   Genders are marked at the right side: pink for female; skyblue for male.
-   Based on LRR intensity on X and Y chromosomes, gender information could checked.

Users could navigate to specific chromosome, gene, region through input box "Chr/Gene":

![Region spectrum panel](./readme_files/fig/5.PNG)

-   "Chr/Gene" accepts:
    -   chromosomes: "1-22, X, Y"
    -   gene symbol: e.g. "TP53"
    -   region: "12:11329569-13134790"
-   The input type is automatically detected, and if it cannot fit into any of the 3 types, whole genome spectrum will be shown.
-   Within the figure, mouse swipe zoom in is supported.
-   Genes from RefSeq are shown as bars at bottom, and the COSMIC cancer genes are highlighted in red.
-   This figure is useful for checking key focal lesions like *IKZF1*, *CDKN2A/B*, *PAX5* etc.

###### ----------------------------END----------------------------
