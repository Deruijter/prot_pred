# 3UTR.fa
# 5UTR.fa
# chromosome1.exon.coords
# chromosome2.exon.coords
# chromosome3.exon.coords
# Marguerat_protein_proliferation.csv
# Marguerat_RNA_proliferation.csv
# pep.fa
# 	pombe_protein_halflife.csv	https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4526151/bin/NIHMS640572-supplement-3.xlsx (TO CSV)
# 	pombe_rna_halflife.csv	https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4770384/bin/MSB-12-857-s003.xlsx (TO CSV)
# Schizosaccharomyces_pombe.ASM294v2.29.cds.all.fa
# 	s_pombe_annot.ASM294v2.29.modified.bed
# 	s_pombe_annot.ASM294v2.29.modified.gff3



import urllib
import gzip
import xlrd
import csv
import sys
import os

def extract(compressed_file_name, uncompressed_file_name):
	with gzip.open(compressed_file_name, 'rb') as compressed_f:	
		with open(uncompressed_file_name, 'w') as out_file:
			out_file.write( compressed_f.read() )

def xlsx_to_csv(xlsx_file_name, csv_file_name):
    xlsx_file = xlrd.open_workbook(xlsx_file_name)
    sheet_names = xlsx_file.sheet_names()
    xlsx_sheet = xlsx_file.sheet_by_name(sheet_names[0])
    csv_file = open(csv_file_name, 'w')
    wr = csv.writer(csv_file, quoting=csv.QUOTE_ALL)

    for rownum in range(xlsx_sheet.nrows):
        wr.writerow(xlsx_sheet.row_values(rownum))
    csv_file.close()

number_of_files = 21

print('Retrieving and extracting or converting {} data files'.format(number_of_files))
print '0/{}'.format(number_of_files),
sys.stdout.flush()


urllib.urlretrieve("ftp://ftp.pombase.org/pombe/genome_sequence_and_features/feature_sequences/UTR/OLD/20160512/3UTR.fa.gz", "/tmp/3UTR.fa.gz")
extract('/tmp/3UTR.fa.gz', '../data/mrnas_p/3UTR.fa')
print '\r1/{}'.format(number_of_files),
sys.stdout.flush()
urllib.urlretrieve("ftp://ftp.pombase.org/pombe/genome_sequence_and_features/feature_sequences/UTR/OLD/20160512/5UTR.fa.gz", "/tmp/5UTR.fa.gz")
extract('/tmp/5UTR.fa.gz', '../data/mrnas_p/5UTR.fa')
print '\r2/{}'.format(number_of_files),
sys.stdout.flush()


urllib.urlretrieve("ftp://ftp.pombase.org/pombe/genome_sequence_and_features/Exon_Coordinates/OLD/20160512/chromosome1.exon.coords", "../data/mrnas_p/chromosome1.exon.coords")
print '\r3/{}'.format(number_of_files),
sys.stdout.flush()
urllib.urlretrieve("ftp://ftp.pombase.org/pombe/genome_sequence_and_features/Exon_Coordinates/OLD/20160512/chromosome2.exon.coords", "../data/mrnas_p/chromosome2.exon.coords")
print '\r4/{}'.format(number_of_files),
sys.stdout.flush()
urllib.urlretrieve("ftp://ftp.pombase.org/pombe/genome_sequence_and_features/Exon_Coordinates/OLD/20160512/chromosome3.exon.coords", "../data/mrnas_p/chromosome3.exon.coords")
print '\r5/{}'.format(number_of_files),
sys.stdout.flush()


urllib.urlretrieve('https://curation.pombase.org/dumps/releases/pombase-chado-v60-2016-05-12/pombe-embl/external_data/Quantitative_gene_expression_data/Marguerat_RNA_proliferation.txt','../data/mrnas_p/Marguerat_RNA_proliferation.csv')
print '\r6/{}'.format(number_of_files),
sys.stdout.flush()
urllib.urlretrieve('https://curation.pombase.org/dumps/releases/pombase-chado-v60-2016-05-12/pombe-embl/external_data/Quantitative_gene_expression_data/Marguerat_protein_proliferation.txt','../data/mrnas_p/Marguerat_protein_proliferation.csv')
print '\r7/{}'.format(number_of_files),
sys.stdout.flush()
urllib.urlretrieve('https://curation.pombase.org/dumps/releases/pombase-chado-v60-2016-05-12/pombe-embl/external_data/Quantitative_gene_expression_data/Marguerat_RNA_quiescence.txt','../data/mrnas_p/Marguerat_RNA_quiescence.csv')
print '\r8/{}'.format(number_of_files),
sys.stdout.flush()
urllib.urlretrieve('https://curation.pombase.org/dumps/releases/pombase-chado-v60-2016-05-12/pombe-embl/external_data/Quantitative_gene_expression_data/Marguerat_protein_quiescence.txt','../data/mrnas_p/Marguerat_protein_quiescence.csv')
print '\r9/{}'.format(number_of_files),
sys.stdout.flush()


urllib.urlretrieve('ftp://ftp.pombase.org/pombe/genome_sequence_and_features/feature_sequences/OLD/20160512/pep.fa.gz','/tmp/pep.fa.gz')
extract('/tmp/pep.fa.gz', '../data/mrnas_p/pep.fa')
print '\r10/{}'.format(number_of_files),
sys.stdout.flush()


urllib.urlretrieve('ftp://ftp.pombase.org/pombe/genome_sequence_and_features/OLD/20151112/Schizosaccharomyces_pombe.ASM294v2.29.cds.fa.gz','/tmp/Schizosaccharomyces_pombe.ASM294v2.29.cds.fa.gz')
extract('/tmp/Schizosaccharomyces_pombe.ASM294v2.29.cds.fa.gz', '../data/mrnas_p/Schizosaccharomyces_pombe.ASM294v2.29.cds.all.fa')
print '\r11/{}'.format(number_of_files),
sys.stdout.flush()


extract('../data/htseq_counts_sub_rna_bt2.txt.tar.gz', '../data/htseq_counts_sub_rna_bt2.txt')
print '\r12/{}'.format(number_of_files),
sys.stdout.flush()
extract('../data/htseq_counts_marguerat_p.txt.tar.gz', '../data/htseq_counts_marguerat_p.txt')
print '\r13/{}'.format(number_of_files),
sys.stdout.flush()
extract('../data/htseq_counts_eser.txt.tar.gz', '../data/htseq_counts_eser.txt')
print '\r14/{}'.format(number_of_files),
sys.stdout.flush()


urllib.urlretrieve('http://www.cell.com/cms/attachment/2078098308/2070960585/mmc3.xlsx','/tmp/pombe_protein_halflife.xlsx')
xlsx_to_csv('/tmp/pombe_protein_halflife.xlsx', '../data/mrnas_p/pombe_protein_halflife.csv')
print '\r15/{}'.format(number_of_files),
sys.stdout.flush()
urllib.urlretrieve('https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4770384/bin/MSB-12-857-s003.xlsx','/tmp/pombe_rna_halflife.xlsx')
xlsx_to_csv('/tmp/pombe_rna_halflife.xlsx', '../data/mrnas_p/pombe_rna_halflife.csv')
print '\r16/{}'.format(number_of_files),
sys.stdout.flush()


extract('../data/mrnas_p/s_pombe_annot.ASM294v2.29.modified.bed.tar.gz', '../data/mrnas_p/s_pombe_annot.ASM294v2.29.modified.bed')
print '\r17/{}'.format(number_of_files),
sys.stdout.flush()


extract('../data/mrnas_p/sub_rna_bt2_minus.bedgraph.tar.gz', '../data/mrnas_p/sub_rna_bt2_minus.bedgraph')
print '\r18/{}'.format(number_of_files),
sys.stdout.flush()
extract('../data/mrnas_p/sub_rna_bt2_plus.bedgraph.tar.gz', '../data/mrnas_p/sub_rna_bt2_plus.bedgraph')
print '\r19/{}'.format(number_of_files),
sys.stdout.flush()
extract('../data/mrnas_p/sub_rpf_bt2_minus.bedgraph.tar.gz', '../data/mrnas_p/sub_rpf_bt2_minus.bedgraph')
print '\r20/{}'.format(number_of_files),
sys.stdout.flush()
extract('../data/mrnas_p/sub_rpf_bt2_plus.bedgraph.tar.gz', '../data/mrnas_p/sub_rpf_bt2_plus.bedgraph')
print '\r21/{}'.format(number_of_files)
sys.stdout.flush()
# NOTE: Modified version of this file is already in the project folder:
# urllib.urlretrieve('ftp://ftp.pombase.org/pombe/genome_sequence_and_features/OLD/20160512/schizosaccharomyces_pombe.chr.gff3','/tmp/s_pombe_annot.ASM294v2.29.gff3')

# Modified using
# with open('/tmp/schizosaccharomyces_pombe.chr.gff3', 'r') as file :
#   filedata = file.read()
# filedata = filedata.replace('Name=SPAC1556.06.2:exon:2;Parent=transcript:SPAC1556.06.1;constitutive=1;ensembl_end_phase=0;ensembl_phase=0;exon_id=SPAC1556.06.2:exon:2;rank=2;version=1', 'Name=SPAC1556.06.1:exon:2;Parent=transcript:SPAC1556.06.1;constitutive=1;ensembl_end_phase=0;ensembl_phase=0;exon_id=SPAC1556.06.1:exon:2;rank=2;version=1')
# with open('../data/mrnas_p/s_pombe_annot.ASM294v2.29.modified.gff3', 'w') as file:
#   file.write(filedata)
# os.system('sox input.wav -b 24 output.aiff rate -v -L -b 90 48k')

# BED FILE MADE USING:
# grep 'ID=transcript:' ../data/mrnas_p/s_pombe_annot.ASM294v2.29.modified.gff3 > /tmp/s_pombe_annot.ASM294v2.29.modified.gff3
# sed -i 's/III/chr3/g' /tmp/s_pombe_annot.ASM294v2.29.modified.gff3
# sed -i 's/II/chr2/g' /tmp/s_pombe_annot.ASM294v2.29.modified.gff3
# sed -i 's/I/chr1/g' /tmp/s_pombe_annot.ASM294v2.29.modified.gff3
# (using bedops:)
# convert2bed --input=gff --output=bed < /tmp/s_pombe_annot.ASM294v2.29.modified.gff3 > ../data/mrnas_p/s_pombe_annot.ASM294v2.29.modified.bed

print('Done.')