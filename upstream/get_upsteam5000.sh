#!/bin/bash
#Get upstream5000.fa.  Clean this up later for public consumption.
DIRECTORY=/mnt/kann2/ds2/data2/SHARE/transfac_data/October_2011/
scp pon2@kann2.bioinf.umbc.edu:$DIRECTORY/upstream5000.fa .

# split them up into chunks of 300 refseqs and rename appropriately
split -l 30300 upstream5000.fa
for file in `ls x*`;do mv $file $file.fa;done
