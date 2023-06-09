
In this folder you will find all files necessary to re-generate our manuscript and supplementary materials 

1. sim-manuscript 

  A folder containing all files necessary to re-generate our manuscript in the style of the journal: Statistics in Medicine. 
  
  - sim-manuscript.Rmd : R markdown file to generate the manuscript
  - sim-manuscript.tex : latex file produced when sim-manuscript.Rmd is knitted
  - WileyNJD-AMA.bst   : styling file generated by the R package "rticles"
  - WileyNJD-v2.cls    : styling file generated by the R package "rticles"
  - bibfile.bib        : bibliography 


2. supplementary_mat.Rmd and supplementary_mat.pdf 

  The R markdown file generating the supplementary materials for this project and the corresponding pdf, respectfully. 


3. figures 

   A folder containing all figures included in both the manuscript and the supplementary materials, stored as .png images. 
   

4. tables 

  A folder containing the necessary files to generate all tables included in both the manuscript and supplementary materials. 
  
  - all_tables.Rmd : includes the code used to generate all tables.
  - all_table.tex  : includes the LaTeX file resulting from knitting the file all_tables.Rmd 
  
  To produce the tables seen in the manuscript and supplementary materials, the LaTeX code for each table in the file all_tables.tex
  was copy and pasted into the appropriate file (sim-manuscript.Rmd or supplemtary_mat.Rmd) and table captions were added once
  pasted into the appropriate document. 
  

