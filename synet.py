from Bio.Blast import NCBIWWW, NCBIXML
result_handle = NCBIWWW.qblast("blastn", "nt", "8332116")
blast_records = NCBIXML.parse(result_handle)
for blast_record in blast_records:
    print(list(blast_record))
