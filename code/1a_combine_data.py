# This document merges protein & rna data from the Marguerat et al paper
# Both files contain a unique system_id per record, 
# where the id of the protein matches the id of the rna (1 on 1)
# Author: Markus de Ruijter

import csv
import urllib
import urllib2
import requests

path_data = '../data/mrnas_p/'
file_protein = 'Marguerat_protein_proliferation.csv'
file_rna = 'Marguerat_RNA_proliferation.csv'
file_protein_seqs = 'pep.fa'
file_rna_cds_seqs = 'Schizosaccharomyces_pombe.ASM294v2.29.cds.all.fa'
file_rna_5utr_seqs = '5UTR.fa'
file_rna_3utr_seqs = '3UTR.fa'
file_out = 'marguerat_rna_vs_prot_proliferation.csv'

data_p = []
data_r = []

def fasta_to_dict(file_name, id_seperator):
	f = open(file_name).readlines()
	parsing = False

	d = []
	key = ''
	seq = ''
	for line in f:
		if line.startswith(">") and parsing:
			d.append({'key':key, 'seq':seq})
			key = ''
			seq = ''
			parsing = False

		if line.startswith(">") and not parsing:
			key = ((line.split(id_seperator))[0])[1:] # Take first part of the id and remove '>'
			parsing = True
			continue
		if parsing:
			seq = seq+line.strip()

	return d


def export_list_of_dict(listOdicts):
	keys = listOdicts[0].keys()
	with open(path_data+file_out,'wb') as out_file:
		dict_writer = csv.DictWriter(out_file, keys)
		dict_writer.writeheader()
		dict_writer.writerows(listOdicts)


def main():
	print 'running combine_data.py'

	lines_p_file = open(path_data+file_protein)
	lines_p = csv.reader(lines_p_file, dialect='excel-tab')

	lines_r_file = open(path_data+file_rna)
	lines_r = csv.reader(lines_r_file, dialect='excel-tab')
	lines_r.next()	# Skip header

	protein_seqs = fasta_to_dict(path_data+file_protein_seqs, '|')
	rna_cds_seqs = fasta_to_dict(path_data+file_rna_cds_seqs, ' ')
	rna_5utr_seqs = fasta_to_dict(path_data+file_rna_5utr_seqs, '|')
	rna_3utr_seqs = fasta_to_dict(path_data+file_rna_3utr_seqs, '|')

	for line in lines_r:
		lines_p_file.seek(0) # reset iterator for the protein file, otherwise we can't search it again
		sys_id = line[0]
		if not sys_id.startswith('SPAC1556.06'):	# Only the RNA file has a transcript number for each transcript
			sys_id_rna = sys_id+'.1'

		avg_count_rna = line[4]
		avg_count_pro = next((item for item in lines_p if item[0] == sys_id), {4:-1})[4]
		seq_rna_cds = next((item for item in rna_cds_seqs if item['key'] == sys_id_rna), {'seq':''})['seq']
		seq_rna_5utr = next((item for item in rna_5utr_seqs if item['key'] == sys_id), {'seq':''})['seq']
		seq_rna_3utr = next((item for item in rna_3utr_seqs if item['key'] == sys_id), {'seq':''})['seq']
		seq_prot = next((item for item in protein_seqs if item['key'] == sys_id), {'seq':''})['seq']
		seq_prot_length = len(seq_prot)
		seq_5utr_length = len(seq_rna_5utr)
		seq_3utr_length = len(seq_rna_3utr)
		translated = 1 if (avg_count_pro > 0) else 0

		# We're only interested in RNA's with a coding region
		if seq_rna_cds == '':
			continue

		data_r.append({'sys_id':line[0]
			,'avg_count_rna':avg_count_rna
			,'avg_count_pro':avg_count_pro
			,'seq_rna_cds':seq_rna_cds
			,'seq_rna_5utr':seq_rna_5utr
			,'seq_rna_3utr':seq_rna_3utr
			,'seq_prot':seq_prot
			,'seq_prot_length':seq_prot_length
			,'seq_5utr_length':seq_5utr_length
			,'seq_3utr_length':seq_3utr_length
			,'translated':translated})
		print sys_id

	export_list_of_dict(data_r)



if __name__ == '__main__':
	main()