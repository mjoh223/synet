from functions.synet_module import create_gff_db
import glob
import os
#Function: find all gff files and create a SQL database for each one
#Input: list of gff files in gff_files/
#output: SQL database for each gff file found
gff_files = glob.glob("/Users/matt/OneDrive/UCSF/JBD_Lab_Rotation/synet_data/gff_files/*.gff")
for gff_file in gff_files:
    db_path = os.path.join("/Users/matt/OneDrive/UCSF/JBD_Lab_Rotation/synet_data/gff_db/", os.path.basename(gff_file[:-4]+".SQL"))
    #Function input: path of .gff file input and path to output sql database
    create_gff_db(gff_file, db_path)
    print("Created database: {}".format(db_path))
