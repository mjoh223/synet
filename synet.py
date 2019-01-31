from Bio.Blast.Applications import NcbiblastpCommandline
target = "/netapp/home/mjohnson/synet/targets/mouse"
query = "/netapp/home/mjohnson/synet/queries/4TNV.fas"

blastp_cline = NcbiblastpCommandline(query=query, db=target, evalue=0.001,outfmt=5, out="out.xml")

blastp_cline()
