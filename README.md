# synet, written by Matthew Johnson 2019
anti-CRISPR "guilt by association" network

Phase 1:
for every target genome, retrive its genomic sequence (fasta) and corresponding gff3 file which provides the metadata.

Pae reference genome data and metadata: https://www.ncbi.nlm.nih.gov/genome/?term=Pseudomonas%20aeruginosa%5BOrganism%5D&cmd=DetailsSearc

Definitions:
1) r is the members of the recognition array
2) D is the refseq Database (D_gff are the gff files and D_fna are the sequences that correspond to each other)

Steps:
1) BLAST all r members against D_fna and return the sequences (S) of the hits above threshold e-value 
2) lookup S in all D_gff and return locus members and add to r array
3) repeat

Scripts:
1) create_BLASTDB.py
2) create_database.py
3) BLAST_database.py

################ -- psuedo code -- ################
#attempt 1
for sequence in [recognition array]:
	[putative sequences] = blast(sequence)
	for seq in [putative sequences]:
		[locus members] = gffLookup(seq)
		[recognition array].append([locus members])

#attempt 2
def blastGff(recognition_array, pointer):
	pointer = len(recognition_array) #2
	for seq in recognition_array:
		putative_seqs_list = blast(seq) #18 in list
		for seq in putative_seqs_list: 
			locus_list = gffLookup(seq) #for each putative seq look in the corresponding gff file to extract locus
			recognition_array.append(locus_list) #18 x 5
		pointer = len(recognition_array) #90+2
		return recognition_array #array of now size 92
			
		blastGFF(recognition_array, pointer) 

#attempt 3
for seq in seed:
	For rec in recognition array:
		locus_members_list = gffLookup(rec) #A, B, C -> A1, A2, B1, B2, C1, C2, C3
		putative_new_members = blast(locus_members_list)
		check to see of new members are % identity to the recognition array

#attempt 4
for acr in recognitionseed
	get sequence of acr
		blast sequence
		get many loci
		align all members in loci to recognitionlist (init with seed members in it)
			if a significant alignment, assign the same ID name
		save to networklist: AxB, AxC, BxC (named by IDs)
		append loci recognitionlist

	#get locus members of acr

