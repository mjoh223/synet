import gffutils
import re
from pyfaidx import Fasta
import glob
from pathlib2 import Path
from Bio.Blast.Applications import NcbiblastpCommandline
import os
import numpy as np

def create_gff_db(gff_file, db_dir):
    #Purpose: create a gff SQL database
    #Input: gff file and location of database to be stored
    #Return: Database created in the location indicated
    gffutils.create_db(gff_file, db_dir, id_spec={"gene":"Dbxref"}, merge_strategy="create_unique") #use the gff file attributes "Name" as the table lookup IDs, If there are duplicates then create uniques by adding _1.

def gff_db_query(db_name, target_gene, window):
    #Purpose: Get information on the target gene's locus
    #Input: Gff SQL database to be looked through, target gene, and size of locus window in nucleotides
    #Return: gffutils Features of the locus members
    db = gffutils.FeatureDB(db_name)
    gene = db[target_gene]
    start, stop = gene.start, gene.end
    seqid = gene[0] #NC_002516.2
    #probably do a check here to see if CDS is in GFF as a featuretype
    neighborFeatures = db.region(region=(seqid, start-window, stop+window), featuretype="gene", completely_within=False)
    neighbor_list = []
    for feat in neighborFeatures:
        neighbor_list.append(feat["Dbxref"][0])
    return neighbor_list #neighborFeatures #[example.attributes["Name"], example.start, example.end]

def recognitionBuild(fasta_file, db_xref_list):
    #Purpose: keep a list of genes of interest and a pointer to what's been used in the BLAST step
    #Input: Genes to be added to recognition list and pointer start
    #Return: Genes to be BlASTed and the pointer end
    genes = Fasta(fasta_file, read_long_names=True, key_function = lambda x: str(re.findall("db_xref=(\w*:\w*)", x)[0])) #clean headers to only be the db_dxref geneid
    for key in db_xref_list:
        print(genes[key].name)
        print(genes[key])

if __name__ == "__main__":
    #blast_db = "/Users/matt/OneDrive/UCSF/JBD_Lab_Rotation/synet_data/targets/JBD30"
    #query = "/Users/matt/OneDrive/UCSF/JBD_Lab_Rotation/synet_data/queries/recognition_seed.fas"
    #blast_results = "/Users/matt/OneDrive/UCSF/JBD_Lab_Rotation/synet_data/out.xml"
    fasta_file = "/Users/matt/OneDrive/UCSF/JBD_Lab_Rotation/synet_data/fna_files/JBD30.fna"
    gff_files = glob.glob("/Users/matt/OneDrive/UCSF/JBD_Lab_Rotation/synet_data/gff_files/*.gff") #this returns all gff files in gff_file/ directory
    recognition_array = ["GeneID:14515639", "GeneID:14515640"]#ACRF1, ACA The idea is that this will be a fasta file with headers with accessions
    gff_db = "/Users/matt/OneDrive/UCSF/JBD_Lab_Rotation/synet_data/gff_db/out_db"

    ################################ -- BLAST -- ##################################
    #blastp_cline = NcbiblastpCommandline(query=query, db=blast_db, evalue=0.1,outfmt=5, out=blast_results)
    #blastp_cline() 
    ###############################################################################
    
    #print(recognitionBuild(fasta_file))

    for accession in recognition_array:
        target_gene = accession
        window = 800
        for gff_file in gff_files:
            db_name = os.path.join("/Users/matt/OneDrive/UCSF/JBD_Lab_Rotation/synet_data/gff_db/", os.path.basename(gff_file[:-4]+".SQL"))
            
            try:

                try:
                    #Step1 extract locus sequences in native genome from the probe ID (dbxref:###)
                    #return a fasta file of sequences found
                    db_xref_list = gff_db_query(db_name, accession, window)
                    print("locus neighbors of {} in {}".format(accession, os.path.basename(gff_file)))
                    print(recognitionBuild(fasta_file, db_xref_list))
                
                    #Step2 blast the fasta file from step 1 and retreive the dbxref from the significant results

                except gffutils.exceptions.FeatureNotFoundError: 
                    print("{} not found in {}".format(accession, os.path.basename(gff_file)))

            except ValueError:
                print("{} db not found, creating now...".format(os.path.basename(gff_file)))
                create_gff_db(gff_file, db_name)
