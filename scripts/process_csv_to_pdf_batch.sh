# Run from the directory containing the input files
dir='~/Documents/Labwork/Projects/Biohub/DRC_neutralization/input/'

for i in *.csv 
 do
	echo "plate file is:" $i
	echo "filename: " "$dir$i"
	/usr/local/bin/Rscript --vanilla ~/Documents/Labwork/Projects/Biohub/DRC_neutralization/scripts/process_csv_to_pdf.R $dir$i &i ~/Documents/Labwork/Projects/Biohub/DRC_neutralization/output/ 

done
 
