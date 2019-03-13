#### Calculation of ED50, RSD values and plotting of corresponding dose-response curves
### Instructions 

## There are three options: 

**1**. A .Rmd R Markdown file that can be opened with RStudio and executed interactivelly 
from within Rstudio. 
*(drc_script.Rmd)*

**2**. A .R script that can be used from the terminal and will produce 2 pdf reports, one for 
each half of the plate (containing what is currently in the output folder of this repository).
*(process_csv_to_pdf.R)*
This script accepts a single .csv file as an input and can be used in the following way :

i) Make sure R is installed and in your path. To install R on Mac go [here](https://cran.r-project.org/bin/macosx/): 
ii) From terminal, assuming that you are in DRC_neutralization/scripts : 

'Rscript --vanilla process_csv_to_pdf.R "full path to input csv" "pdf_label" "full path to output directory"'
  
Output pdf will be named: "pdf_label"_left_plate_reg and "pdf_label"_right_plate_reg


**3**. A .R script that can be used from the terminal and will produce 2 pdf reports, one for 
each half of the plate. The script requires two parameters, the input folder and the output folder. 
It generates two pdf reports (one for each half of the plate) in the output folder for ALL .csv files
in the input folder (NOTE: if you have unrelated csvs in that folder the script will not work as it 
currently performs the same operation on ALL .csv files in the input folder). 
*(process_csv_to_pdf_batch.R)*

In addition, the script generates a master table of all assay attributes (ED50, SDR, etc) in the output folder that can be
used to query across all experiments. 

Usage (assuming you are in the scripts folder):

'Rscript --vanilla process_csv_to_pdf_batch.R "input_folder" "output_folder"'
