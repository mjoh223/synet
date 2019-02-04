from BCBio import GFF
from BCBio.GFF import GFFExaminer

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

if __name__ == "__main__":
    fasta_file = ""
    gff_file = ""
    seq_dict = fasta_parse(fasta_file)
    gff_parse(gff_file, seq_dict)
