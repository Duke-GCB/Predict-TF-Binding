from TF_enrichment import apply_model_to_seqs

kmers = [3]
seqlist = ['GCCCGCAGAGCGGAAGGCGGGATGGCTGGGGGCGGG',
           'GGGCTCAGCGCCGACTGCGCGCCTCTGCCCGCGAAA',
           'CGCCATAGCGACGGCGCCGCAATTTAGGAGCGTGCT',
           'TCCCGCCTTCCGCTCTGCGGGCGGCAGCCGGGCTGG']
modelfile = 'ELK1_100nM_Bound_filtered_normalized_GGAA_1a2a3mer_format.model'
results = apply_model_to_seqs(seqlist, modelfile)
print results

