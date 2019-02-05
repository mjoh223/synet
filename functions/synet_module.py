#from BCBio import GFF
#from BCBio.GFF import GFFExaminer
#from Bio import SeqIO
import gffutils
from pyfaidx import Fasta
import glob
from pathlib2 import Path
from Bio.Blast.Applications import NcbiblastpCommandline
import os

def create_gff_db(gff_file, db_name):
    db = gffutils.create_db(gff_file, db_name, id_spec=["Name"], merge_strategy="create_unique")

def gff_db(db_name, target_gene, window, fasta_file):
    db = gffutils.FeatureDB(db_name)
    gene = db[target_gene]#NP_248700.1
    start = gene.start
    stop = gene.end
    seqid = gene[0] #NC_002516.2
    #probably do a check here to see if CDS is in GFF as a featuretype
    neighborFeatures = list(db.region(region=(seqid, start-window, stop+window), featuretype="CDS", completely_within=False))
    return neighborFeatures #[example.attributes["Name"], example.start, example.end]

def extract_fasta(locus_ids, fasta_file):
    fasta = Fasta(Fasta_file)
    for ids in locus_ids:
        seq = fasta[ids]
    return seq


if __name__ == "__main__":
    blast_db = "/Users/matt/OneDrive/UCSF/JBD_Lab_Rotation/targets/pae_blastdb"
    query = "/Users/matt/OneDrive/UCSF/JBD_Lab_Rotation/queries/aca1.fas"
    blastp_cline = NcbiblastpCommandline(query=query, db=blast_db, evalue=0.1,outfmt=5, out="out.xml")
    blastp_cline() 
    
    fasta_file = "/Users/matt/OneDrive/UCSF/JBD_Lab_Rotation/metadata/GCF_000006765.1_ASM676v1_protein.faa"

    gff_files = glob.glob("/Users/matt/OneDrive/UCSF/JBD_Lab_Rotation/metadata/*.gff")
    target_gene = "YP_007392343.1"#"NP_248723.1" add something here to not care about the .1
    window = 800
    for gff_file in gff_files:
        db_file = Path(os.path.basename(gff_file))
        db_name = os.path.basename(gff_file)
        try:
            my_abs_path = db_file.resolve(strict=True)
        except:# FileNotFoundError:
            print("Creating Database")
            create_gff_db(gff_file=gff_file, db_name=db_name)
            print("\tCreated db: {}".format(db_name))
        else:
            print("Database Found")
            try:
                features = gff_db(db_name=db_name, target_gene=target_gene, window=window, fasta_file=fasta_file)
                headers = [f["Name"] for f in features]
                flat_headers = [x for sub_list in headers for x in sub_list]
                format_list = [len(features), target_gene, window]
                print("Found {} genes around probe {} using a window of {}bp.".format(*format_list))
                print([f[6] for f in features])
                print([f["Parent"] for f in features])
            except gffutils.exceptions.FeatureNotFoundError:
                format_list = [target_gene, db_name]
                print("{} not found in {}".format(*format_list))





    
