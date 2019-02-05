#from BCBio import GFF
#from BCBio.GFF import GFFExaminer
#from Bio import SeqIO
import gffutils
from pyfaidx import Fasta

def gff_parse(gff_file, seq_dict):
    in_handle = open(gff_file)
    rec_list = []
    for rec in GFF.parse(in_handle, base_dict = seq_dict):
        rec_list.append(rec)
    in_handle.close()
    return rec_list

def fasta_parse(fasta_file):
    in_seq_handle = open(fasta_file)
    seq_dict = SeqIO.to_dict(SeqIO.parse(in_seq_handle, "fasta"))
    in_seq_handle.close()
    return seq_dict

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
    #seq_dict = fasta_parse(fasta_file)
    #print(gff_parse(gff_file, seq_dict)[30])
    features = gff_db(db_name="pae_db", target_gene="NP_248723.1", window=10000, fasta_file=fasta_file)
    headers = [f["Name"] for f in features]
    flat_headers = [x for sub_list in headers for x in sub_list]
    print(features[1])








    
