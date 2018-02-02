require(rstudioapi)
setwd(paste(dirname(getActiveDocumentContext()$path), '/..', sep=''))

require(stringi)
require(dplyr)

mrnas_p = readRDS('data/mrnas_p.RObject')
permuted_pro_per_rna = readRDS('./data/permuted_pro_per_rna.RObject')
# REMOVE LINE BELOW IF YOU WANT TO CREATE ALL PERMUTATION DATA SETS
permuted_pro_per_rna = data.frame('X0'=permuted_pro_per_rna$X0) # X0 is the non-permuted (original) data



##############################################################################################
##### GENERAL FUNCTIONS ######################################################################

# GENERATE DATA FRAME WITH CODON | AA | PROPERTY | UNIQUE SYMBOL
# (N:nonpolar, P:polar, B:basic, A:acidic, X:stopcodon)
CalcCodonsExtended = function(){
	bases = c('T','C','A','G')
	codon_info = as.data.frame(matrix(ncol=4,nrow=0))
	colnames(codon_info) = c('codon','aa','prop','uniq')
	aas = 'FFLLSSSSYYXXCCXWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
	props = 'NNNNPPPPPPXXPPXNNNNNNNNNBBPPBBBBNNNNPPPPPPBBPPBBNNNNNNNNAAAANNNN'# N:nonpolar, P:polar, B:basic, A:acidic, X:stopcodon
	uniq_symbol = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789abcdefghijklmnopqrstuvwxyz.-'
	i = 1
	for(b1 in bases){
		for(b2 in bases){
			for(b3 in bases){
				codon = sprintf('%s%s%s',b1,b2,b3)
				codon_info = rbind(codon_info, data.frame(codon=codon, aa=substr(aas,i,i), prop=substr(props,i,i), uniq=substr(uniq_symbol,i,i)))
				i = i+1
			}
		}
	}
	
	codon_info[,'codon_rna'] = gsub('T','U',codon_info[,'codon'])
	codon_info[,'anticodon'] = stri_reverse(chartr('ATGC','TACG',codon_info[,'codon']))
	
	return(codon_info)
}
codon_info = CalcCodonsExtended()

# From http://stats.stackexchange.com/a/164830/114397
# Modified it a bit though
FindPeaks <- function (x, m = 1){
	shape <- diff(sign(diff(x, na.pad = FALSE)))
	pks <- sapply(which(shape < 0), FUN = function(i){
		i = i+1
		z <- i - m
		z <- ifelse(z > 0, z, 1)
		w <- i + m
		w <- ifelse(w < length(x), w, length(x))
		if(all(x[c(z : (i-1))] < x[i]) & all(x[c((i + 1) : w)] <= x[i])) return(i) else return(numeric(0))
	})
	pks <- unlist(pks)
	pks
}

# Calculate how much a variable is "per codon" (e.g. translation rate per codon). 
# Uses codon RATIO in sequence, not the absolute amount.
CalcVarPerCodonRatio_sub = function(rna_codons, var){
	t = matrix(nrow=64, ncol=1)
	t[,1] = 0
	rownames(t) = codon_info[,'codon']
	
	t2 = table(rna_codons)
	t2 = t2 / sum(t2)
	t2 = t2 * var
	t = as.matrix(c(t2, t[!rownames(t)%in%rownames(t2),]))
	
	t = t[order(rownames(t))]
	return(t)
}
CalcVarPerCodonRatio = function(d, var_name=NULL, codon_subset=NULL){
	if(is.null(var_name)){
		# Used for counting/summing the codon ratios
		var_name = 'ratio_sum'
		temp = apply(d, 1, function(x) {
			codons = unlist(x$codons)
			if(!is.null(codon_subset)){
				codons = codons[codon_subset[1]:codon_subset[2]]
			}
			CalcVarPerCodonRatio_sub(codons, 1)
		})
	} else {
		temp = apply(d, 1, function(x) {
			codons = unlist(x$codons)
			if(!is.null(codon_subset)){
				codons = codons[codon_subset[1]:codon_subset[2]]
			}
			CalcVarPerCodonRatio_sub(codons, as.numeric(x[var_name]))
		})
	}
	
	rownames(temp) = codon_info[order(codon_info['codon']), 'codon']
	temp = rowSums(temp)
	
	temp = as.data.frame(temp)
	colnames(temp) = c(var_name)
	temp$codon = rownames(temp)
	return(temp)
}

GetStructureData = function(d_struct, w=14){
  # Get delta G info using RNALL (still need to find out whether to use RNALL or RNALfold but lets start with RNALL for now)
  f = '3rd_party/Rnall2.0/temp.fa'
  
  for(row_num in 1:nrow(mrnas_p)){
    #row_num = 5
    #w=14
    #seq_subset='init'
    #sys_id = mrna$sys_id
    print(row_num)
    
    mrna = mrnas_p[row_num,]
    
    write(sprintf('>%s\n%s%s%s', 'temp', mrna$seq_rna_5utr, mrna$seq_rna_cds, mrna$seq_rna_3utr), file=f)
    system(sprintf('3rd_party/Rnall2.0/Rnall2.0.exe -i 3rd_party/Rnall2.0/temp.fa -o 3rd_party/Rnall2.0/temp.out -g -w %s -l %s', w, w), wait=T)
    dg = c()
    
    # TRY AGAIN IF THE FILE WASN'T WRITTEN PROPERLY
    i = 0;
    while(i<10){
      dg = tryCatch(read.table('3rd_party/Rnall2.0/temp.out', header=F, sep='\t', skip=1), error=function(e) {return(FALSE)})
      if(dg==F){
        print('oh noes, file isnt ready!')
        i = i+1
        Sys.sleep(0.3)
        system(sprintf('3rd_party/Rnall2.0/Rnall2.0.exe -i 3rd_party/Rnall2.0/temp.fa -o 3rd_party/Rnall2.0/temp.out -g -w %s -l %s', w, w))
        if(i == 10){
          stop('Error reading or writing file')
        }
      } else {
        i = 10
      }
    }
    
    d_struct = rbind(d_struct, data.frame('sys_id'=mrna$sys_id, 'struct'=dg, 'w'=w))
    
  }
  return(d_struct)
}



##############################################################################################
##### DATA SET COMPONENTS ####################################################################

# B. NUCLEOTIDE DATA
AppendNucleotideData = function(d, seq_subset='full'){
	nucleotides = c('A','T','G','C')
	
	for(row_num in 1:nrow(mrnas_p)){
		mrna = mrnas_p[row_num,]
		
		seq = switch(seq_subset,
					 full={sprintf('%s%s%s', mrna$seq_rna_5utr, mrna$seq_rna_cds, mrna$seq_rna_3utr)},
					 utr5={mrna$seq_rna_5utr},
					 cds ={mrna$seq_rna_cds},
					 utr3={mrna$seq_rna_5utr},
					 cds_s={substring(mrna$seq_rna_cds, 4, 30)}) # first 9 codons after start codon
		
		count_total = nchar(sprintf('%s', seq))
		for(n in nucleotides){
			count_n = length(gregexpr(n, sprintf('%s', seq), fixed=T)[[1]]) 
			d[row_num, sprintf('nt_total_ratio_%s_%s', seq_subset, tolower(n))] = count_n / count_total
		}
	}
	
	# some 5' 3' UTRs are 0 length
	d[d[,sprintf('nt_total_ratio_%s_a', seq_subset)]==Inf,sprintf('nt_total_ratio_%s_a', seq_subset)] = median(d[d[,sprintf('nt_total_ratio_%s_a', seq_subset)]!=Inf,sprintf('nt_total_ratio_%s_a', seq_subset)])
	d[d[,sprintf('nt_total_ratio_%s_t', seq_subset)]==Inf,sprintf('nt_total_ratio_%s_t', seq_subset)] = median(d[d[,sprintf('nt_total_ratio_%s_t', seq_subset)]!=Inf,sprintf('nt_total_ratio_%s_t', seq_subset)])
	d[d[,sprintf('nt_total_ratio_%s_g', seq_subset)]==Inf,sprintf('nt_total_ratio_%s_g', seq_subset)] = median(d[d[,sprintf('nt_total_ratio_%s_g', seq_subset)]!=Inf,sprintf('nt_total_ratio_%s_g', seq_subset)])
	d[d[,sprintf('nt_total_ratio_%s_c', seq_subset)]==Inf,sprintf('nt_total_ratio_%s_c', seq_subset)] = median(d[d[,sprintf('nt_total_ratio_%s_c', seq_subset)]!=Inf,sprintf('nt_total_ratio_%s_c', seq_subset)])
	
	return(d)
}

# C. CODON DATA
AppendCodonDataNoStop = function(d, seq_subset='cds'){
	codon_info_no_stop = codon_info[codon_info$prop!='X',]
	for(row_num in 1:nrow(mrnas_p)){
		mrna = mrnas_p[row_num,]
		
		codon_subset = switch(seq_subset,
							  cds={NULL},
							  cds_s={c(2, 10)}, # first 9 codons after start codon
							  cds_s2={c(2, 51)}) # 51 codons as recommended by tuller 2011
		
		codons = CalcVarPerCodonRatio(mrna, NULL, codon_subset)
		for(ci in codon_info_no_stop$codon){
			d[row_num, sprintf('codon_total_ratio_%s_%s', seq_subset, tolower(ci))] = codons[codons$codon==ci,'ratio_sum']
		}
	}
	return(d)
}

# C. CODON DATA
AppendCodonDataNoStopNoStart = function(d, seq_subset='cds'){
	codon_info_no_stop = codon_info[codon_info$prop!='X',]
	for(row_num in 1:nrow(mrnas_p)){
		mrna = mrnas_p[row_num,]
		
		codon_subset = switch(seq_subset,
							  cds={c(2, length(mrna$codons[[1]]))},
							  cds_s={c(2, 10)}, # first 9 codons after start codon
							  cds_s2={c(2, 51)}) # 51 codons as recommended by tuller 2011
		
		codons = CalcVarPerCodonRatio(mrna, NULL, codon_subset)
		for(ci in codon_info_no_stop$codon){
			d[row_num, sprintf('codon_total_ratio_%s_%s', seq_subset, tolower(ci))] = codons[codons$codon==ci,'ratio_sum']
		}
	}
	return(d)
}

# D. STRUCTURE DATA
AppendStructureData = function(d, seq_subset='full', w=14){
  
  filename = sprintf('data/d_struct_split%s.RObject',w)
  d_struct_split = data.frame('sys_id'=character(),'struct'=character(),'w'=integer())
  if(file.exists(filename)){
    d_struct_split = readRDS(filename)
  } else {
    d_struct_split = GetStructureData(d_struct_split, 14)
    d_struct_split = split(d_struct_split, d_struct_split$sys_id, drop=T)
    saveRDS(d_struct_split, file=filename)
  }
  
	# Get delta G info using RNALL (still need to find out whether to use RNALL or RNALfold but lets start with RNALL for now)
	f = '3rd_party/Rnall2.0/temp.fa'
	
	temp = paste(rep(c('G','C'), each=w), collapse='')
	write(sprintf('>%s\n%s', 'test', temp), file=f)
	system(sprintf('3rd_party/Rnall2.0/Rnall2.0.exe -i 3rd_party/Rnall2.0/temp.fa -o 3rd_party/Rnall2.0/temp.out -g -w %s -l %s', w, w))
	dg = read.csv('3rd_party/Rnall2.0/temp.out', header=F, sep='\t', skip=1)
	colnames(dg) = c('window_size','center_position','delta_g')
	
	dg_window_max = max(dg$delta_g) 
	dg_window_min = min(dg$delta_g) # minimum free energy for this window size
	
	for(row_num in 1:nrow(mrnas_p)){
		mrna = mrnas_p[row_num,]
		
		dg = d_struct_split[as.character(mrna$sys_id)][[1]]
		
		colnames(dg) = c('sys_id','window_size','center_position','delta_g')
		dg = switch(seq_subset,
					full={dg},
					utr5={ dg[findInterval(dg$center_position, c(1, mrna$seq_5utr_length+1))==1,] },
					cds ={ dg[findInterval(dg$center_position, c(mrna$seq_5utr_length+1, mrna$seq_5utr_length+(mrna$seq_prot_length*3+3+1)))==1,]},
					utr3={ dg[findInterval(dg$center_position, c(mrna$seq_5utr_length+(mrna$seq_prot_length*3+3+1)))==1,]},
					cds_s ={ dg[findInterval(dg$center_position, c(mrna$seq_5utr_length+1, mrna$seq_5utr_length+1+30))==1,]}, # Include the start codon here because structure spans multiple bases
					init = {dg[findInterval(dg$center_position, c(mrna$seq_5utr_length+1-5, mrna$seq_5utr_length+1))==1,] }, # last 5 positions of 5' utr
					cds_s2={ dg[findInterval(dg$center_position, c(mrna$seq_5utr_length+1, mrna$seq_5utr_length+1+102))==1,]} # As mentioned by tuller 2011
		)
		
		dg_upper_75 = dg[dg$delta_g>dg_window_max*0.75,]
		dg_upper_50 = dg[dg$delta_g>dg_window_max*0.50,]
		dg_upper_25 = dg[dg$delta_g>dg_window_max*0.25,]
		dg_upper_10 = dg[dg$delta_g>dg_window_max*0.10,]
		dg_lower_75 = dg[dg$delta_g<dg_window_min*0.75,]
		dg_lower_50 = dg[dg$delta_g<dg_window_min*0.50,]
		dg_lower_25 = dg[dg$delta_g<dg_window_min*0.25,]
		dg_lower_10 = dg[dg$delta_g<dg_window_min*0.10,]
		
		dg_peaks = dg[FindPeaks(dg$delta_g,m=round(w/4,0)),]
		dg_valleys = dg[FindPeaks(-dg$delta_g,m=round(w/4,0)),]
		max_peaks = round(nchar(paste(c(mrna$seq_rna_5utr, mrna$seq_rna_cds, mrna$seq_rna_3utr), collapse='')) / (w/4), 0)
		dg_upper_75_peaks = nrow(dg_peaks[dg_peaks$delta_g > dg_window_max*0.75,]) / max_peaks
		dg_upper_50_peaks = nrow(dg_peaks[dg_peaks$delta_g > dg_window_max*0.50,]) / max_peaks
		dg_upper_25_peaks = nrow(dg_peaks[dg_peaks$delta_g > dg_window_max*0.25,]) / max_peaks
		dg_upper_10_peaks = nrow(dg_peaks[dg_peaks$delta_g > dg_window_max*0.10,]) / max_peaks
		dg_lower_75_peaks = nrow(dg_valleys[dg_valleys$delta_g < dg_window_min*0.75,]) / max_peaks
		dg_lower_50_peaks = nrow(dg_valleys[dg_valleys$delta_g < dg_window_min*0.50,]) / max_peaks
		dg_lower_25_peaks = nrow(dg_valleys[dg_valleys$delta_g < dg_window_min*0.25,]) / max_peaks
		dg_lower_10_peaks = nrow(dg_valleys[dg_valleys$delta_g < dg_window_min*0.10,]) / max_peaks
		
		dg_sum = summary(dg$delta_g)
		total_rows = nrow(dg)
		d[row_num, 'lowest'] = dg_sum['Min.']
		d[row_num, 'highest'] = dg_sum['Max.']
		d[row_num, 'mean'] = dg_sum['Mean']
		d[row_num, 'median'] = dg_sum['Median']
		d[row_num, 'q1'] = dg_sum['1st Qu.']
		d[row_num, 'q2'] = dg_sum['3rd Qu.']
		d[row_num, 'upper_10'] = nrow(dg_upper_10) / total_rows
		d[row_num, 'upper_25'] = nrow(dg_upper_25) / total_rows
		d[row_num, 'upper_50'] = nrow(dg_upper_50) / total_rows
		d[row_num, 'upper_75'] = nrow(dg_upper_75) / total_rows
		d[row_num, 'lower_10'] = nrow(dg_lower_10) / total_rows
		d[row_num, 'lower_25'] = nrow(dg_lower_25) / total_rows
		d[row_num, 'lower_50'] = nrow(dg_lower_50) / total_rows
		d[row_num, 'lower_75'] = nrow(dg_lower_75) / total_rows
		d[row_num, 'upper_75_peaks'] = dg_upper_75_peaks
		d[row_num, 'upper_50_peaks'] = dg_upper_50_peaks
		d[row_num, 'upper_25_peaks'] = dg_upper_25_peaks
		d[row_num, 'upper_10_peaks'] = dg_upper_10_peaks
		d[row_num, 'lower_75_peaks'] = dg_lower_75_peaks
		d[row_num, 'lower_50_peaks'] = dg_lower_50_peaks
		d[row_num, 'lower_25_peaks'] = dg_lower_25_peaks
		d[row_num, 'lower_10_peaks'] = dg_lower_10_peaks
	}
	
	d[!is.finite(d$highest), 'highest'] = median(d[is.finite(d$highest), 'highest'])
	d[!is.finite(d$mean), 'mean'] = median(d[is.finite(d$mean), 'mean'])
	d[!is.finite(d$median), 'median'] = median(d[is.finite(d$median), 'median'])
	d[!is.finite(d$q1), 'q1'] = median(d[is.finite(d$q1), 'q1'])
	d[!is.finite(d$q2), 'q2'] = median(d[is.finite(d$q2), 'q2'])
	d[!is.finite(d$upper_10), 'upper_10'] = median(d[is.finite(d$upper_10), 'upper_10'])
	d[!is.finite(d$upper_25), 'upper_25'] = median(d[is.finite(d$upper_25), 'upper_25'])
	d[!is.finite(d$upper_50), 'upper_50'] = median(d[is.finite(d$upper_50), 'upper_50'])
	d[!is.finite(d$upper_75), 'upper_75'] = median(d[is.finite(d$upper_75), 'upper_75'])
	d[!is.finite(d$lower_10), 'lower_10'] = median(d[is.finite(d$lower_10), 'lower_10'])
	d[!is.finite(d$lower_25), 'lower_25'] = median(d[is.finite(d$lower_25), 'lower_25'])
	d[!is.finite(d$lower_50), 'lower_50'] = median(d[is.finite(d$lower_50), 'lower_50'])
	d[!is.finite(d$lower_75), 'lower_75'] = median(d[is.finite(d$lower_75), 'lower_75'])
	# Note "lowest" is used as reference column to see whether there is any data at all
	d[!is.finite(d$lowest), 'upper_75_peaks'] = median(d[is.finite(d$lowest), 'upper_75_peaks'])
	d[!is.finite(d$lowest), 'upper_50_peaks'] = median(d[is.finite(d$lowest), 'upper_50_peaks'])
	d[!is.finite(d$lowest), 'upper_25_peaks'] = median(d[is.finite(d$lowest), 'upper_25_peaks'])
	d[!is.finite(d$lowest), 'upper_10_peaks'] = median(d[is.finite(d$lowest), 'upper_10_peaks'])
	d[!is.finite(d$lowest), 'lower_75_peaks'] = median(d[is.finite(d$lowest), 'lower_75_peaks'])
	d[!is.finite(d$lowest), 'lower_50_peaks'] = median(d[is.finite(d$lowest), 'lower_50_peaks'])
	d[!is.finite(d$lowest), 'lower_25_peaks'] = median(d[is.finite(d$lowest), 'lower_25_peaks'])
	d[!is.finite(d$lowest), 'lower_10_peaks'] = median(d[is.finite(d$lowest), 'lower_10_peaks'])
	d[!is.finite(d$lowest), 'lowest'] = median(d[is.finite(d$lowest), 'lowest'])
	
	row_prefix = sprintf('dg_%s_%s', seq_subset, w)
	colnames(d) = c(colnames(d)[1:(length(colnames(d))-22)], paste(row_prefix, c(colnames(d)[(length(colnames(d))-21):length(colnames(d))]), sep='_') )
	
	return(d)
}

# E. 5'/3' ORF DATA
AppendUorfData = function(d){
	for(row_num in 1:nrow(mrnas_p)){
		mrna = mrnas_p[row_num,]
		
		d[row_num, 'orfs_5utr_fr0_atgs'] = mrna$orfs_5utr_fr0_atgs
		d[row_num, 'orfs_5utr_fr1_atgs'] = mrna$orfs_5utr_fr1_atgs
		d[row_num, 'orfs_5utr_fr2_atgs'] = mrna$orfs_5utr_fr2_atgs
		d[row_num, 'orfs_3utr_fr0_atgs'] = mrna$orfs_3utr_fr0_atgs
		d[row_num, 'orfs_3utr_fr1_atgs'] = mrna$orfs_3utr_fr1_atgs
		d[row_num, 'orfs_3utr_fr2_atgs'] = mrna$orfs_3utr_fr2_atgs
		
	}
	
	return(d)
}
AppendUorfData2 = function(d){
	for(row_num in 1:nrow(mrnas_p)){
		mrna = mrnas_p[row_num,]
		
		d[row_num, 'orfs_5utr_atgs'] = mrna$orfs_5utr_fr0_atgs + mrna$orfs_5utr_fr1_atgs + mrna$orfs_5utr_fr2_atgs
		d[row_num, 'orfs_3utr_atgs'] = mrna$orfs_3utr_fr0_atgs + mrna$orfs_3utr_fr1_atgs + mrna$orfs_3utr_fr2_atgs
		
	}
	
	return(d)
}

# F. HALF LIFE DATA
AppendHalflifeData = function(d, rsr=T, rhl=T, phl=T){
	for(row_num in 1:nrow(mrnas_p)){
		mrna = mrnas_p[row_num,]
		
		if(rsr){ d[row_num, 'synthesistime'] = mrna$synthesistime }
		if(rhl){ d[row_num, 'halflife'] = mrna$halflife }
		if(phl){ d[row_num, 'p_halflife'] = mrna$p_halflife}
	}
	return(d)
}

# G. AMINO ACID DATA
AppendAminoAcidDataNoStop = function(d, seq_subset='cds'){
	codon_info_no_stop = codon_info[codon_info$prop!='X',]
	for(row_num in 1:nrow(mrnas_p)){
		mrna = mrnas_p[row_num,]
		
		codon_subset = switch(seq_subset,
							  cds={NULL},
							  cds_s={c(2, 10)}, # first 9 codons after start codon
							  cds_s2={c(2,39)}) # First 39 codons as mentioned by Tuller 2011
		
		codons = CalcVarPerCodonRatio(mrna, NULL, codon_subset)
		
		for(amino_acid in unique(codon_info_no_stop$aa)){
			d[row_num, sprintf('amino_acid_ratio_%s_%s', seq_subset, tolower(amino_acid))] = sum(codons[codons$codon%in%codon_info_no_stop[codon_info_no_stop$aa==amino_acid,'codon'],'ratio_sum'])
		}
	}
	return(d)
}
# G. AMINO ACID DATA
AppendAminoAcidDataNoStopNoStart = function(d, seq_subset='cds'){
	codon_info_no_stop = codon_info[codon_info$prop!='X',]
	for(row_num in 1:nrow(mrnas_p)){
		mrna = mrnas_p[row_num,]
		
		codon_subset = switch(seq_subset,
							  cds={c(2, length(mrna$codons[[1]]))},
							  cds_s={c(2, 10)}, # first 9 codons after start codon
							  cds_s2={c(2,39)}) # First 39 codons as mentioned by Tuller 2011
		
		codons = CalcVarPerCodonRatio(mrna, NULL, codon_subset)
		
		for(amino_acid in unique(codon_info_no_stop$aa)){
			d[row_num, sprintf('amino_acid_ratio_%s_%s', seq_subset, tolower(amino_acid))] = sum(codons[codons$codon%in%codon_info_no_stop[codon_info_no_stop$aa==amino_acid,'codon'],'ratio_sum'])
		}
	}
	return(d)
}
# G. AMINO ACID PROPERTY DATA
AppendAminoAcidPropDataNoStop = function(d, seq_subset='cds'){
	codon_info_no_stop = codon_info[codon_info$prop!='X',]
	for(row_num in 1:nrow(mrnas_p)){
		mrna = mrnas_p[row_num,]
		
		codon_subset = switch(seq_subset,
							  cds={NULL},
							  cds_s={c(2, 10)}, # first 9 codons after start codon
							  cds_s2={c(2,39)}) # First 39 codons as mentioned by Tuller 2011
		
		codons = CalcVarPerCodonRatio(mrna, NULL, codon_subset)
		
		for(amino_acid_prop in unique(codon_info_no_stop$prop)){
			d[row_num, sprintf('amino_acid_prop_ratio_%s_%s', seq_subset, tolower(amino_acid_prop))] = sum(codons[codons$codon%in%codon_info_no_stop[codon_info_no_stop$prop==amino_acid_prop,'codon'],'ratio_sum'])
		}
	}
	return(d)
}
# G. AMINO ACID PROPERTY DATA
AppendAminoAcidPropDataNoStopNoStart = function(d, seq_subset='cds'){
	codon_info_no_stop = codon_info[codon_info$prop!='X',]
	for(row_num in 1:nrow(mrnas_p)){
		mrna = mrnas_p[row_num,]
		
		codon_subset = switch(seq_subset,
							  cds={c(2, length(mrna$codons[[1]]))},
							  cds_s={c(2, 10)}, # first 9 codons after start codon
							  cds_s2={c(2,39)}) # First 39 codons as mentioned by Tuller 2011
		
		codons = CalcVarPerCodonRatio(mrna, NULL, codon_subset)
		
		for(amino_acid_prop in unique(codon_info_no_stop$prop)){
			d[row_num, sprintf('amino_acid_prop_ratio_%s_%s', seq_subset, tolower(amino_acid_prop))] = sum(codons[codons$codon%in%codon_info_no_stop[codon_info_no_stop$prop==amino_acid_prop,'codon'],'ratio_sum'])
		}
	}
	return(d)
}
# G. AMINO ACID DATA POLY PROLINE
AppendAminoAcidPprData = function(d){
	for(row_num in 1:nrow(mrnas_p)){
		mrna = mrnas_p[row_num,]
		
		# Get [G|D|A|P] PP [P|W|N|D|G] from AA sequence
		matches = gregexpr('(?=([G|D|A|P]PP|PP[P|W|N|D|G]))',mrna$seq_prot, perl=T)[[1]]
		pps = length(matches[matches>-1 & matches<34]) # Only get those at the start of the cds
		
		d[row_num, 'pp_count'] = pps
	}
	return(d)
}

# I. RPF DATA
AppendRpfData = function(d, rpf_per_rna=T, rpf_regression=T) {
  for(row_num in 1:nrow(mrnas_p)){
    mrna = mrnas_p[row_num,]
    
    if(rpf_per_rna) { d[row_num, 'rpf_per_rna'] = mrna$rpf_per_rna }
    if(rpf_regression) { d[row_num, 'rpf_regression'] = mrna$rpf_regression }
  }
  return(d)
}



##############################################################################################
##### INDIVIDUAL DATA SETS ###################################################################

# BASE DATA SET
data_set_empty = mrnas_p[,c('pro_per_rna','sys_id')]

# A1 - seq lengths
CreateA1 = function(){
  data_set_a1 = mrnas_p[,c('pro_per_rna','sys_id','seq_5utr_length','seq_prot_length','seq_3utr_length')]
  
  for(i in 1:length(permuted_pro_per_rna)){
    postfix = ''
    if(i > 1) {
      postfix = sprintf('_p%s',(i-1))
      data_set_a1$pro_per_rna = permuted_pro_per_rna[,i]
    }
    
    write.table(data_set_a1, sprintf('./feature_sets/data_set_a1%s.csv',postfix), sep='\t', row.names=F)
  }
}

# B - Nucleotides
CreateBSets = function(){
  data_set_b1 = AppendNucleotideData(data_set_empty)
  data_set_b2 = AppendNucleotideData(data_set_empty, 'utr5')
  data_set_b3 = AppendNucleotideData(data_set_empty, 'cds')
  data_set_b4 = AppendNucleotideData(data_set_empty, 'utr3')
  data_set_b5 = AppendNucleotideData(data_set_empty, 'cds_s')
  
  for(i in 1:length(permuted_pro_per_rna)){
    postfix = ''
    if(i > 1) {
      postfix = sprintf('_p%s',(i-1))
      data_set_b1$pro_per_rna = permuted_pro_per_rna[,i]
      data_set_b2$pro_per_rna = permuted_pro_per_rna[,i]
      data_set_b3$pro_per_rna = permuted_pro_per_rna[,i]
      data_set_b4$pro_per_rna = permuted_pro_per_rna[,i]
      data_set_b5$pro_per_rna = permuted_pro_per_rna[,i]
    }
    
    write.table(data_set_b1, sprintf('./feature_sets/data_set_b1%s.csv',postfix), sep='\t', row.names=F)
    write.table(data_set_b2, sprintf('./feature_sets/data_set_b2%s.csv',postfix), sep='\t', row.names=F)
    write.table(data_set_b3, sprintf('./feature_sets/data_set_b3%s.csv',postfix), sep='\t', row.names=F)
    write.table(data_set_b4, sprintf('./feature_sets/data_set_b4%s.csv',postfix), sep='\t', row.names=F)
    write.table(data_set_b5, sprintf('./feature_sets/data_set_b5%s.csv',postfix), sep='\t', row.names=F)
  }
}

# C - Codons
CreateCSets = function(){
  data_set_c5 = AppendCodonDataNoStop(data_set_empty, 'cds_s')
  data_set_c6 = AppendCodonDataNoStop(data_set_empty, 'cds_s2')
  data_set_c7 = AppendCodonDataNoStopNoStart(data_set_empty)
  
  for(i in 1:length(permuted_pro_per_rna)){
    postfix = ''
    if(i > 1) {
      postfix = sprintf('_p%s',(i-1))
      data_set_c5$pro_per_rna = permuted_pro_per_rna[,i]
      data_set_c6$pro_per_rna = permuted_pro_per_rna[,i]
      data_set_c7$pro_per_rna = permuted_pro_per_rna[,i]
    }
    
    write.table(data_set_c5, sprintf('./feature_sets/data_set_c5%s.csv',postfix), sep='\t', row.names=F)
    write.table(data_set_c6, sprintf('./feature_sets/data_set_c6%s.csv',postfix), sep='\t', row.names=F)
    write.table(data_set_c7, sprintf('./feature_sets/data_set_c7%s.csv',postfix), sep='\t', row.names=F)
  }
}

# D14 - Struct window size 14
CreateD14Sets = function(){
  
  data_set_db1 = AppendStructureData(data_set_empty, 'full', 14)
  data_set_db2 = AppendStructureData(data_set_empty, 'utr5', 14)
  data_set_db3 = AppendStructureData(data_set_empty, 'cds', 14)
  data_set_db4 = AppendStructureData(data_set_empty, 'utr3', 14)
  data_set_db5 = AppendStructureData(data_set_empty, 'cds_s', 14)
  data_set_db6 = AppendStructureData(data_set_empty, 'init', 14)
  data_set_db7 = AppendStructureData(data_set_empty, 'cds_s2', 14)
  
  for(i in 1:length(permuted_pro_per_rna)){
    postfix = ''
    if(i > 1) {
      postfix = sprintf('_p%s',(i-1))
      data_set_db1$pro_per_rna = permuted_pro_per_rna[,i]
      data_set_db2$pro_per_rna = permuted_pro_per_rna[,i]
      data_set_db3$pro_per_rna = permuted_pro_per_rna[,i]
      data_set_db4$pro_per_rna = permuted_pro_per_rna[,i]
      data_set_db5$pro_per_rna = permuted_pro_per_rna[,i]
      data_set_db6$pro_per_rna = permuted_pro_per_rna[,i]
      data_set_db7$pro_per_rna = permuted_pro_per_rna[,i]
    }
    
    write.table(data_set_db1, sprintf('./feature_sets/data_set_db1%s.csv',postfix), sep='\t', row.names=F)
    write.table(data_set_db2, sprintf('./feature_sets/data_set_db2%s.csv',postfix), sep='\t', row.names=F)
    write.table(data_set_db3, sprintf('./feature_sets/data_set_db3%s.csv',postfix), sep='\t', row.names=F)
    write.table(data_set_db4, sprintf('./feature_sets/data_set_db4%s.csv',postfix), sep='\t', row.names=F)
    write.table(data_set_db5, sprintf('./feature_sets/data_set_db5%s.csv',postfix), sep='\t', row.names=F)
    write.table(data_set_db6, sprintf('./feature_sets/data_set_db6%s.csv',postfix), sep='\t', row.names=F)
    write.table(data_set_db7, sprintf('./feature_sets/data_set_db7%s.csv',postfix), sep='\t', row.names=F)
  }
}

# E - SORFs
CreateESets = function(){
  data_set_e1 = AppendUorfData(data_set_empty)
  data_set_e2 = AppendUorfData2(data_set_empty)
  
  for(i in 1:length(permuted_pro_per_rna)){
    postfix = ''
    if(i > 1) {
      postfix = sprintf('_p%s',(i-1))
      data_set_e1$pro_per_rna = permuted_pro_per_rna[,i]
      data_set_e2$pro_per_rna = permuted_pro_per_rna[,i]
    }
    
    write.table(data_set_e1, sprintf('./feature_sets/data_set_e1%s.csv',postfix), sep='\t', row.names=F)
    write.table(data_set_e2, sprintf('./feature_sets/data_set_e2%s.csv',postfix), sep='\t', row.names=F)
  }
}

# F - Half life
CreateFSets = function(){
  data_set_f1 = AppendHalflifeData(data_set_empty)
  data_set_f2 = AppendHalflifeData(data_set_empty, T, F, F) # rna synth
  data_set_f3 = AppendHalflifeData(data_set_empty, F, T, F) # rna hl
  data_set_f4 = AppendHalflifeData(data_set_empty, F, F, T) # pro hl
  
  for(i in 1:length(permuted_pro_per_rna)){
    postfix = ''
    if(i > 1) {
      postfix = sprintf('_p%s',(i-1))
      data_set_f1$pro_per_rna = permuted_pro_per_rna[,i]
      data_set_f2$pro_per_rna = permuted_pro_per_rna[,i]
      data_set_f3$pro_per_rna = permuted_pro_per_rna[,i]
      data_set_f4$pro_per_rna = permuted_pro_per_rna[,i]
    }
    
    write.table(data_set_f1, sprintf('./feature_sets/data_set_f1%s.csv',postfix), sep='\t', row.names=F)
    write.table(data_set_f2, sprintf('./feature_sets/data_set_f2%s.csv',postfix), sep='\t', row.names=F)
    write.table(data_set_f3, sprintf('./feature_sets/data_set_f3%s.csv',postfix), sep='\t', row.names=F)
    write.table(data_set_f4, sprintf('./feature_sets/data_set_f4%s.csv',postfix), sep='\t', row.names=F)
  }
}

# G - Amino Acids
CreateGSets = function(){
  data_set_g5 = AppendAminoAcidPprData(data_set_empty)
  data_set_g9 = AppendAminoAcidDataNoStop(data_set_empty, 'cds_s')
  data_set_gA = AppendAminoAcidDataNoStop(data_set_empty, 'cds_s2')
  data_set_gC = AppendAminoAcidPropDataNoStop(data_set_empty, 'cds_s')
  data_set_gD = AppendAminoAcidPropDataNoStop(data_set_empty, 'cds_s2')
  data_set_gE = AppendAminoAcidDataNoStopNoStart(data_set_empty)
  data_set_gF = AppendAminoAcidPropDataNoStopNoStart(data_set_empty)
  
  for(i in 1:length(permuted_pro_per_rna)){
    postfix = ''
    if(i > 1) {
      postfix = sprintf('_p%s',(i-1))
      data_set_g5$pro_per_rna = permuted_pro_per_rna[,i]
      data_set_g9$pro_per_rna = permuted_pro_per_rna[,i]
      data_set_gA$pro_per_rna = permuted_pro_per_rna[,i]
      data_set_gC$pro_per_rna = permuted_pro_per_rna[,i]
      data_set_gD$pro_per_rna = permuted_pro_per_rna[,i]
      data_set_gE$pro_per_rna = permuted_pro_per_rna[,i]
      data_set_gF$pro_per_rna = permuted_pro_per_rna[,i]
    }
    
    write.table(data_set_g5, sprintf('./feature_sets/data_set_g5%s.csv',postfix), sep='\t', row.names=F)
    write.table(data_set_g9, sprintf('./feature_sets/data_set_g9%s.csv',postfix), sep='\t', row.names=F)
    write.table(data_set_gA, sprintf('./feature_sets/data_set_gA%s.csv',postfix), sep='\t', row.names=F)
    write.table(data_set_gC, sprintf('./feature_sets/data_set_gC%s.csv',postfix), sep='\t', row.names=F)
    write.table(data_set_gD, sprintf('./feature_sets/data_set_gD%s.csv',postfix), sep='\t', row.names=F)
    write.table(data_set_gE, sprintf('./feature_sets/data_set_gE%s.csv',postfix), sep='\t', row.names=F)
    write.table(data_set_gF, sprintf('./feature_sets/data_set_gF%s.csv',postfix), sep='\t', row.names=F)
  }
}

# H1 - RNA abundance
CreateH1 = function(){
  data_set_h1 = mrnas_p[,c('pro_per_rna','sys_id', 'avg_count_rna')]
  
  for(i in 1:length(permuted_pro_per_rna)){
    postfix = ''
    if(i > 1) {
      postfix = sprintf('_p%s',(i-1))
      data_set_h1$pro_per_rna = permuted_pro_per_rna[,i]
    }
    
    write.table(data_set_h1, sprintf('./feature_sets/data_set_h1%s.csv',postfix), sep='\t', row.names=F)
  }
}

# I - RPF
CreateISets = function(){
  data_set_i1 = AppendRpfData(data_set_empty)
  data_set_i2 = AppendRpfData(data_set_empty, T, F) # rpf per rna
  data_set_i3 = AppendRpfData(data_set_empty, F, T) # rpf regression
  
  for(i in 1:length(permuted_pro_per_rna)){
    postfix = ''
    if(i > 1) {
      postfix = sprintf('_p%s',(i-1))
      data_set_i1$pro_per_rna = permuted_pro_per_rna[,i]
      data_set_i2$pro_per_rna = permuted_pro_per_rna[,i]
      data_set_i3$pro_per_rna = permuted_pro_per_rna[,i]
    }
    
    write.table(data_set_i1, sprintf('./feature_sets/data_set_i1%s.csv',postfix), sep='\t', row.names=F)
    write.table(data_set_i2, sprintf('./feature_sets/data_set_i2%s.csv',postfix), sep='\t', row.names=F)
    write.table(data_set_i3, sprintf('./feature_sets/data_set_i3%s.csv',postfix), sep='\t', row.names=F)
  }
}

# FULL
CreateFullSet = function(){
  # This data set is used for the automatic feature selection
  data_set_a1 = mrnas_p[,c('pro_per_rna','sys_id','seq_5utr_length','seq_prot_length','seq_3utr_length')]
  data_set_ultra_mega = AppendNucleotideData(data_set_a1, 'full')
  data_set_ultra_mega = AppendNucleotideData(data_set_ultra_mega, 'utr5')
  data_set_ultra_mega = AppendNucleotideData(data_set_ultra_mega, 'cds')
  data_set_ultra_mega = AppendNucleotideData(data_set_ultra_mega, 'utr3')
  data_set_ultra_mega = AppendNucleotideData(data_set_ultra_mega, 'cds_s')
  
  data_set_ultra_mega = AppendCodonDataNoStopNoStart(data_set_ultra_mega, 'cds')
  data_set_ultra_mega = AppendCodonDataNoStop(data_set_ultra_mega, 'cds_s')
  data_set_ultra_mega = AppendCodonDataNoStop(data_set_ultra_mega, 'cds_s2')
  
  data_set_ultra_mega = AppendStructureData(data_set_ultra_mega, 'full', 14)
  data_set_ultra_mega = AppendStructureData(data_set_ultra_mega, 'utr5', 14)
  data_set_ultra_mega = AppendStructureData(data_set_ultra_mega, 'cds', 14)
  data_set_ultra_mega = AppendStructureData(data_set_ultra_mega, 'utr3', 14)
  data_set_ultra_mega = AppendStructureData(data_set_ultra_mega, 'cds_s', 14)
  data_set_ultra_mega = AppendStructureData(data_set_ultra_mega, 'init', 14)
  data_set_ultra_mega = AppendStructureData(data_set_ultra_mega, 'cds_s2', 14)
  
  data_set_ultra_mega = AppendStructureData(data_set_ultra_mega, 'full', 28)
  data_set_ultra_mega = AppendStructureData(data_set_ultra_mega, 'utr5', 28)
  data_set_ultra_mega = AppendStructureData(data_set_ultra_mega, 'cds', 28)
  data_set_ultra_mega = AppendStructureData(data_set_ultra_mega, 'utr3', 28)
  data_set_ultra_mega = AppendStructureData(data_set_ultra_mega, 'cds_s', 28)
  data_set_ultra_mega = AppendStructureData(data_set_ultra_mega, 'init', 28)
  data_set_ultra_mega = AppendStructureData(data_set_ultra_mega, 'cds_s2', 28)
  
  data_set_ultra_mega = AppendStructureData(data_set_ultra_mega, 'full', 50)
  data_set_ultra_mega = AppendStructureData(data_set_ultra_mega, 'utr5', 50)
  data_set_ultra_mega = AppendStructureData(data_set_ultra_mega, 'cds', 50)
  data_set_ultra_mega = AppendStructureData(data_set_ultra_mega, 'utr3', 50)
  data_set_ultra_mega = AppendStructureData(data_set_ultra_mega, 'cds_s', 50)
  data_set_ultra_mega = AppendStructureData(data_set_ultra_mega, 'init', 50)
  data_set_ultra_mega = AppendStructureData(data_set_ultra_mega, 'cds_s2', 50)
  
  data_set_ultra_mega = AppendUorfData(data_set_ultra_mega)
  data_set_ultra_mega = AppendUorfData2(data_set_ultra_mega)
  data_set_ultra_mega = AppendHalflifeData(data_set_ultra_mega)
  
  data_set_ultra_mega = AppendAminoAcidDataNoStopNoStart(data_set_ultra_mega)
  data_set_ultra_mega = AppendAminoAcidDataNoStop(data_set_ultra_mega, 'cds_s')
  data_set_ultra_mega = AppendAminoAcidPropDataNoStopNoStart(data_set_ultra_mega)
  data_set_ultra_mega = AppendAminoAcidPropDataNoStop(data_set_ultra_mega, 'cds_s')
  data_set_ultra_mega = AppendAminoAcidPprData(data_set_ultra_mega)  #PolyProline
  data_set_ultra_mega = AppendAminoAcidDataNoStop(data_set_ultra_mega, 'cds_s2')
  data_set_ultra_mega = AppendAminoAcidPropDataNoStop(data_set_ultra_mega, 'cds_s2')
  
  data_set_ultra_mega[,'avg_count_rna'] = mrnas_p[,'avg_count_rna'] # H1
  data_set_ultra_mega = AppendRpfData(data_set_ultra_mega)
  
  write.table(data_set_ultra_mega, './feature_sets/data_set_ultra_mega2.csv', sep='\t', row.names=F)
  data_set_ultra_mega = read.table('./feature_sets/data_set_ultra_mega2.csv', sep='\t', header = T)
}

CreateA1()
CreateBSets()
CreateCSets()
CreateD14Sets()
CreateESets()
CreateFSets()
CreateGSets()
CreateH1()
CreateISets()
CreateFullSet()



##############################################################################################
##### (MANUAL) COMBINED DATA SETS ############################################################


# seq, nuc, cod, aac, aap, ppr, str, orf, abn
data_set_x0 = mrnas_p[,c('pro_per_rna','sys_id','seq_5utr_length','seq_prot_length','seq_3utr_length')]
data_set_x0 = AppendNucleotideData(data_set_x0, 'cds')
data_set_x0 = AppendCodonDataNoStopNoStart(data_set_x0)
data_set_x0 = AppendAminoAcidDataNoStopNoStart(data_set_x0)
data_set_x0 = AppendAminoAcidPropDataNoStopNoStart(data_set_x0)
data_set_x0 = AppendAminoAcidPprData(data_set_x0)
data_set_x0 = AppendStructureData(data_set_x0, 'cds', 14)
data_set_x0 = AppendUorfData2(data_set_x0)
data_set_x0$avg_count_rna = mrnas_p$avg_count_rna
write.table(data_set_x0, './feature_sets/data_set_x0_111111111.csv', sep='\t', row.names=F)

# seq, nuc, cod, aac, aap, ---, str, ---, ---
data_set_x0 = mrnas_p[,c('pro_per_rna','sys_id','seq_5utr_length','seq_prot_length','seq_3utr_length')]
data_set_x0 = AppendNucleotideData(data_set_x0, 'cds')
data_set_x0 = AppendCodonDataNoStopNoStart(data_set_x0)
data_set_x0 = AppendAminoAcidDataNoStopNoStart(data_set_x0)
data_set_x0 = AppendAminoAcidPropDataNoStopNoStart(data_set_x0)
#data_set_x0 = AppendAminoAcidPprData(data_set_x0)
data_set_x0 = AppendStructureData(data_set_x0, 'cds', 14)
#data_set_x0 = AppendUorfData2(data_set_x0)
#data_set_x0$avg_count_rna = mrnas_p$avg_count_rna
write.table(data_set_x0, './feature_sets/data_set_x0_111110100.csv', sep='\t', row.names=F)

# seq, nuc, cod, ---, ---, ---, str, ---, ---
data_set_x0 = mrnas_p[,c('pro_per_rna','sys_id','seq_5utr_length','seq_prot_length','seq_3utr_length')]
data_set_x0 = AppendNucleotideData(data_set_x0, 'cds')
data_set_x0 = AppendCodonDataNoStopNoStart(data_set_x0)
#data_set_x0 = AppendAminoAcidDataNoStopNoStart(data_set_x0)
#data_set_x0 = AppendAminoAcidPropDataNoStopNoStart(data_set_x0)
#data_set_x0 = AppendAminoAcidPprData(data_set_x0)
data_set_x0 = AppendStructureData(data_set_x0, 'cds', 14)
#data_set_x0 = AppendUorfData2(data_set_x0)
#data_set_x0$avg_count_rna = mrnas_p$avg_count_rna
write.table(data_set_x0, './feature_sets/data_set_x0_111000100.csv', sep='\t', row.names=F)

# seq, nuc, ---, aac, ---, ---, str, ---, ---
data_set_x0 = mrnas_p[,c('pro_per_rna','sys_id','seq_5utr_length','seq_prot_length','seq_3utr_length')]
data_set_x0 = AppendNucleotideData(data_set_x0, 'cds')
#data_set_x0 = AppendCodonDataNoStopNoStart(data_set_x0)
data_set_x0 = AppendAminoAcidDataNoStopNoStart(data_set_x0)
#data_set_x0 = AppendAminoAcidPropDataNoStopNoStart(data_set_x0)
#data_set_x0 = AppendAminoAcidPprData(data_set_x0)
data_set_x0 = AppendStructureData(data_set_x0, 'cds', 14)
#data_set_x0 = AppendUorfData2(data_set_x0)
#data_set_x0$avg_count_rna = mrnas_p$avg_count_rna
write.table(data_set_x0, './feature_sets/data_set_x0_110100100.csv', sep='\t', row.names=F)

# seq, nuc, cod, ---, aap, ---, str, ---, ---
data_set_x0 = mrnas_p[,c('pro_per_rna','sys_id','seq_5utr_length','seq_prot_length','seq_3utr_length')]
data_set_x0 = AppendNucleotideData(data_set_x0, 'cds')
data_set_x0 = AppendCodonDataNoStopNoStart(data_set_x0)
data_set_x0 = AppendAminoAcidDataNoStopNoStart(data_set_x0)
#data_set_x0 = AppendAminoAcidPropDataNoStopNoStart(data_set_x0)
#data_set_x0 = AppendAminoAcidPprData(data_set_x0)
data_set_x0 = AppendStructureData(data_set_x0, 'cds', 14)
#data_set_x0 = AppendUorfData2(data_set_x0)
#data_set_x0$avg_count_rna = mrnas_p$avg_count_rna
write.table(data_set_x0, './feature_sets/data_set_x0_111010100.csv', sep='\t', row.names=F)

# seq, nuc, ---, aac, aap, ---, str, ---, ---
data_set_x0 = mrnas_p[,c('pro_per_rna','sys_id','seq_5utr_length','seq_prot_length','seq_3utr_length')]
data_set_x0 = AppendNucleotideData(data_set_x0, 'cds')
#data_set_x0 = AppendCodonDataNoStopNoStart(data_set_x0)
data_set_x0 = AppendAminoAcidDataNoStopNoStart(data_set_x0)
data_set_x0 = AppendAminoAcidPropDataNoStopNoStart(data_set_x0)
#data_set_x0 = AppendAminoAcidPprData(data_set_x0)
data_set_x0 = AppendStructureData(data_set_x0, 'cds', 14)
#data_set_x0 = AppendUorfData2(data_set_x0)
#data_set_x0$avg_count_rna = mrnas_p$avg_count_rna
write.table(data_set_x0, './feature_sets/data_set_x0_110110100.csv', sep='\t', row.names=F)

# seq, nuc, cod, aac, ---, ---, str, ---, ---
data_set_x0 = mrnas_p[,c('pro_per_rna','sys_id','seq_5utr_length','seq_prot_length','seq_3utr_length')]
data_set_x0 = AppendNucleotideData(data_set_x0, 'cds')
data_set_x0 = AppendCodonDataNoStopNoStart(data_set_x0)
data_set_x0 = AppendAminoAcidDataNoStopNoStart(data_set_x0)
#data_set_x0 = AppendAminoAcidPropDataNoStopNoStart(data_set_x0)
#data_set_x0 = AppendAminoAcidPprData(data_set_x0)
data_set_x0 = AppendStructureData(data_set_x0, 'cds', 14)
#data_set_x0 = AppendUorfData2(data_set_x0)
#data_set_x0$avg_count_rna = mrnas_p$avg_count_rna
write.table(data_set_x0, './feature_sets/data_set_x0_111100100.csv', sep='\t', row.names=F)




#### DATA SETS WITH EXPERIMENTAL DATA (BASED ON BEST NO-EXPERIMENTAL-DATA DATA SET)####

best_no_exp_data_data_set = mrnas_p[,c('pro_per_rna','sys_id','seq_5utr_length','seq_prot_length','seq_3utr_length')]
best_no_exp_data_data_set = AppendNucleotideData(best_no_exp_data_data_set, 'cds')
best_no_exp_data_data_set = AppendAminoAcidDataNoStopNoStart(best_no_exp_data_data_set)
best_no_exp_data_data_set = AppendAminoAcidPropDataNoStopNoStart(best_no_exp_data_data_set)
best_no_exp_data_data_set = AppendStructureData(best_no_exp_data_data_set, 'cds', 14)

# RPF, RSR, RHL, PHL
data_set_x1 = AppendRpfData(best_no_exp_data_data_set)
data_set_x1 = AppendHalflifeData(data_set_x1)
write.table(data_set_x1, './feature_sets/data_set_x0_110110100-1111.csv', sep='\t', row.names=F)

# RPF, ---, ---, ---
data_set_x1 = AppendRpfData(best_no_exp_data_data_set)
write.table(data_set_x1, './feature_sets/data_set_x0_110110100-1000.csv', sep='\t', row.names=F)

# ---, RSR, ---, ---
data_set_x1 = AppendHalflifeData(best_no_exp_data_data_set, T, F, F)
write.table(data_set_x1, './feature_sets/data_set_x0_110110100-0100.csv', sep='\t', row.names=F)

# ---, ---, RHL, ---
data_set_x1 = AppendHalflifeData(best_no_exp_data_data_set, F, T, F)
write.table(data_set_x1, './feature_sets/data_set_x0_110110100-0010.csv', sep='\t', row.names=F)

# ---, ---, ---, PHL
data_set_x1 = AppendHalflifeData(best_no_exp_data_data_set, F, F, T)
write.table(data_set_x1, './feature_sets/data_set_x0_110110100-0001.csv', sep='\t', row.names=F)



##############################################################################################
##### (AUTOMATIC) COMBINED DATA SETS #########################################################

# See the machine learning python script






