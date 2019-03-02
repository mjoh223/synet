#!/usr/bin/env python
import glob, subprocess
#purpose: combine all fasta files and create a single blastdb for the concatenated fasta file
#input: list of all fna files in fna_files/
#output: blastdb of all fna files

#open up a new file called cat.fna
f = open("/Users/matt/OneDrive/UCSF/JBD_Lab_Rotation/synet_data/fna_files/cat.fna", "a+")

#append all of the fasta files from fna_files/ to cat.fna
for file in glob.glob("/Users/matt/OneDrive/UCSF/JBD_Lab_Rotation/synet_data/fna_files/*"):
    f.write(open(file).read()) #dirty way to create master fna file probably will need to use seq objects instead, also it assumes cat is empty

#call makeblastdb using the cat.fa
subprocess.call("/anaconda3/envs/synet/bin/makeblastdb -in " +
        "/Users/matt/OneDrive/UCSF/JBD_Lab_Rotation/synet_data/fna_files/cat.fna " + 
        "-out /Users/matt/OneDrive/UCSF/JBD_Lab_Rotation/synet_data/targets/synet_blastdb " + 
        "-dbtype prot ", shell=True)
