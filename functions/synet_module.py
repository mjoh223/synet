#from BCBio import GFF
#from BCBio.GFF import GFFExaminer
#from Bio import SeqIO
import gffutils
from pyfaidx import Fasta
import os.path
from pathlib2 import Path

def create_gff_db(gff_file, db_name):
    db = gffutils.create_db(gff_file, db_name, id_spec=["Name"], merge_strategy=create_unique)

def gff_db(db_name, target_gene, window, fasta_file):
    db = gffutils.FeatureDB(db_name)
    gene = db[target_gene]#NP_248700.1
    start = gene.start
    stop = gene.end
    tmp = list(db.region(region=('NC_002516.2', start-window, stop+window), featuretype="CDS", completely_within=False))
    return tmp #[example.attributes["Name"], example.start, example.end]

def extract_fasta(locus_ids, fasta_file):
    fasta = Fasta(Fasta_file)
    for ids in locus_ids:
        seq = fasta[ids]
    return seq


if __name__ == "__main__":
    fasta_file = "/Users/matt/OneDrive/UCSF/JBD_Lab_Rotation/metadata/GCF_000006765.1_ASM676v1_protein.faa"
    gff_file = "/Users/matt/OneDrive/UCSF/JBD_Lab_Rotation/metadata/GCF_000006765.1_ASM676v1_genomic.gff"
    my_file = Path("pae_db") 
    try:
        my_abs_path = my_file.resolve(strict=True)
    except FileNotFoundError:
        print("Creating Database\n")
        create_gff_db(gff_file=gff_file, db_name="pae_db")
    else:
        print("Database Found\n")
        features = gff_db(db_name="pae_db", target_gene="NP_248723.1", window=10000, fasta_file=fasta_file)
    headers = [f["Name"] for f in features]
    flat_headers = [x for sub_list in headers for x in sub_list]
    print(features[1])








    
