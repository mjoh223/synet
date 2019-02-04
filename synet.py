from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
from functions.synet_module import gff_parse
#target = "/netapp/home/mjohnson/targets/mouse"
#query = "/netapp/home/mjohnson/queries/sample.fas"

#blastn_cline = NcbiblastnCommandline(query=query, db=target, evalue=0.01,outfmt=5, out="out.xml")

#blastn_cline()

#result_handle = open("out.xml")
#blast_record = NCBIXML.read(result_handle)

#E_VALUE_THRESH = 0.01

#for alignment in blast_record.alignments:
#	for hsp in alignment.hsps:
#		if hsp.expect < E_VALUE_THRESH:
#                       print("****Alignment****")
#			print("sequence:", alignment.title)
#			print("length:", alignment.length)
#			print("e value:", hsp.expect)
#			print(hsp.query[0:75] + "...")
#			print(hsp.match[0:75] + "...")
#			print(hsp.sbjct[0:75] + "...")
gff_filename = "/Users/matt/OneDrive/UCSF/JBD_Lab_Rotation/metadata/GCF_000006765.1_ASM676v1_genomic.gff"

print(gff_parse(gff_filename))
