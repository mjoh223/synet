from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
from functions.synet_module import create_gff_db
import glob
import os
#purpose: blast sequences in recognition array and retrieve the locus
#input: list of members in the recognition array to be used as query
#output: fasta headers (>name) for the targets that met threshold (which will go to the next step for gff db lookup)

query = "/Users/matt/OneDrive/UCSF/JBD_Lab_Rotation/synet_data/queries/recognition_seed.fas"
target= "/Users/matt/OneDrive/UCSF/JBD_Lab_Rotation/synet_data/targets/synet_blastdb"
out_xml = "/Users/matt/OneDrive/UCSF/JBD_Lab_Rotation/synet_data/blast_xml/out.xml"

blastn_cline = NcbiblastnCommandline(query=query, db=target, evalue=0.01,outfmt=5, out="out.xml")
blastn_cline()

result_handle = open("out.xml")
blast_record = NCBIXML.read(result_handle)

E_VALUE_THRESH = 0.01

for alignment in blast_record.alignments:
	for hsp in alignment.hsps:
		if hsp.expect < E_VALUE_THRESH:
                       print("****Alignment****")
			print("sequence:", alignment.title)
			print("length:", alignment.length)
			print("e value:", hsp.expect)
			print(hsp.query[0:75] + "...")
			print(hsp.match[0:75] + "...")
			print(hsp.sbjct[0:75] + "...")
