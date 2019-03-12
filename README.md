# DRC_neutralization
DRC script for Leslie's team 

# Instructions 
There are two options: 

1. A .Rmd R Markdown file that can be opened with RStudio and executed interactivelly 
from within Rstudio. 

2. A .R script that can be used from the terminal and will produce 2 pdf reports, one for 
each half of the plate (containing what is currently in the output folder of this repository). 
This script accepts a single .csv file as an input and can be used in the following way :

i) Make sure R is installed and in your path. To install R on Mac go here: 
https://cran.r-project.org/bin/macosx/

ii) From terminal, assuming that you are in DRC_neutralization/scripts : 
$ Rscript --vanilla process_csv_to_pdf.R <full path to input csv> <pdf_label> <full path to output directory>
  
Output pdf will be named   <pdf_label>_left_plate_reg and <pdf_label>_right_plate_reg
