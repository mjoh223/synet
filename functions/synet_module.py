import gffutils
import re
from pyfaidx import Fasta
import itertools
import csv
import glob
from pathlib2 import Path
from Bio.Blast.Applications import NcbitblastxCommandline
from Bio.Blast import NCBIXML
import os
import numpy as np
from Bio import pairwise2

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

def faidx(fasta_file):
    #Purpose: keep a list of genes of interest and a pointer to what's been used in the BLAST step
    #Input: Genes to be added to recognition list and pointer start
    #Return: Genes to be BlASTed and the pointer end
    genes = Fasta(fasta_file, read_long_names=True, key_function = lambda x: str(re.findall("db_xref=(\w*:\w*)", x)[0])) #clean headers to only be the db_dxref geneid
    #for key in db_xref_list:
    #    print(genes[key].name)
    #    print(genes[key])
    return genes

if __name__ == "__main__":
    gff_files = glob.glob("/Users/matt/OneDrive/UCSF/JBD_Lab_Rotation/synet_data/gff_files/*.gff") #this returns all gff files in gff_file/ directory
    fna_files = glob.glob("/Users/matt/OneDrive/UCSF/JBD_Lab_Rotation/synet_data/fna_files/*.fna")
    recognition_seed = ["GeneID:14515639", "GeneID:14515640"]#ACRF1, ACA The idea is that this will be a fasta file with headers with accessions
    gff_db = "/Users/matt/OneDrive/UCSF/JBD_Lab_Rotation/synet_data/gff_db/out_db"
    synet_db = "/Users/matt/OneDrive/UCSF/JBD_Lab_Rotation/synet_data/targets/tmp" 
    blast_out = "/Users/matt/OneDrive/UCSF/JBD_Lab_Rotation/synet_data/blast_xml/out.xml"
    ################################ -- BLAST -- ##################################
    #blastp_cline = NcbiblastpCommandline(query=query, db=blast_db, evalue=0.1,outfmt=5, out=blast_results)
    #blastp_cline() 
    ###############################################################################
    
    #print(recognitionBuild(fasta_file))
    window = 3000
    f= open("/Users/matt/OneDrive/UCSF/JBD_Lab_Rotation/synet_data/queries/recognitionlist.fa","w+")
    recognition_seed_rawseq = []
    for acr in recognition_seed:
        #get sequence for acr -- DONE --
        for fna in fna_files:
            genes = faidx(fna)
            try:
                #write to fasta file
                f.write(">{}\n".format(genes[acr].name))
                f.write("{}\n".format(str(genes[acr])))
                recognition_seed_rawseq.append([acr, str(genes[acr]) ])
            except KeyError:
                continue
    f.close()
    #blast sequence -- DONE --

    blastx_cline = NcbitblastxCommandline(query="/Users/matt/OneDrive/UCSF/JBD_Lab_Rotation/synet_data/queries/recognitionlist.fa", db=synet_db, evalue=0.00001, outfmt=5, out=blast_out)
    blastx_cline()
    #get many loci -- DONE --
    result_handle = open(blast_out)
    blast_records = NCBIXML.parse(result_handle)
    EVT = 0.00001
    locus_anchor = []
    for blast_record in blast_records:
        for alignment in blast_record.alignments:
            anchor = re.findall("db_xref=(\w*:\w*)", alignment.title)[0]
            locus_anchor.append(anchor)

    loci_list = []
    network_list = []
    for gff_file in gff_files:
        db_name = os.path.join("/Users/matt/OneDrive/UCSF/JBD_Lab_Rotation/synet_data/gff_db/", os.path.basename(gff_file[:-4]+".SQL"))
        for anchor in locus_anchor:
            try:
                many_loci = gff_db_query(db_name, anchor, window)
                row = list(itertools.combinations(many_loci, 2))
                network_list.append(many_loci)
                loci_list.append(many_loci)
            except gffutils.exceptions.FeatureNotFoundError:
                continue
    ##working on this to create a network data structure
    mylist = [list(itertools.combinations(loci, 2)) for loci in loci_list]
    network_array = np.asarray(mylist)
    #print(network_array.reshape(2,-2))
    #np.savetxt("foo.csv", np.asarray(network_list).flatten(), delimiter="\t")
    #align all members in loci to recognitionlist (init with seed members) -- DONE --
    recognitionlist = []
    recognitionlist.append(recognition_seed)
    recognitionlist = sum(recognitionlist, [])
    for fna in fna_files:
        genes = faidx(fna)
        for loci in loci_list:
            for member in loci:
                try:
                    query = str(genes[member])
                    for ref in recognition_seed_rawseq:
                        aln_score = pairwise2.align.localxx(query, ref[1], score_only=True)/max(len(ref[1]),len(query))
                        print(aln_score)
                        #if a significant alignment, assign the same ID name -- working on it --
                        #if aln_score > .8:
                           # print([member, ref[0], aln_score])
                except:
                    continue

    for loci in loci_list:
        for l in loci:
            if l not in recognitionlist:
                recognitionlist.append(l)
    print(recognitionlist)









        #for gff_file in gff_files:

         #   try:
          #      try:
          #          #locus = gff_db_query(db_name, acr, window)
          #          print("locus neighbors of {} found in {}:".format(acr, os.path.basename(gff_file)))
          #          for fna_file in fna_files:
          #              seqs = recognitionBuild(fna_file, locus)
          #              print(seqs)
                    #take the seeds and align them to the native locus members
                    #reference = seqs[acr]
                    #for xref in db_xref_list:
                    #    query = seqs[xref]
                    #    alignments = pairwise2.align.globalxx(str(query), str(reference), score_only=True)
                    #    print( [xref, acr, alignments] ) #want to eventually use protein sequences
                    #Step2 blast the fasta file from step 1 and retreive the dbxref from the significant results
         #       except gffutils.exceptions.FeatureNotFoundError: 
         #           print("{} not found in {}".format(acr, os.path.basename(gff_file)))
         #   except ValueError:
         #       print("{} db not found, creating now...".format(os.path.basename(gff_file)))
        #      create_gff_db(gff_file, db_name)
