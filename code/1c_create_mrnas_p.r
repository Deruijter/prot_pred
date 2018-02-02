require(rstudioapi)
setwd(paste(dirname(getActiveDocumentContext()$path), '/..', sep=''))

require(dplyr)
require(stringi)

GetMrnaLength = function(mrna){
  return(mrna$seq_5utr_length+mrna$seq_3utr_length+mrna$seq_prot_length*3+3)
}

FilterMrnasFromMargueratData = function(d){
  d[,'seq_rna_cds'] = as.character(d[,'seq_rna_cds'])
  d[,'avg_count_rna'] = as.double(d[,'avg_count_rna'])
  d = d[d$avg_count_rna > 0,]  # Remove 0 and lower RNA values (they screw up the plots, plus I don't think we can asume to find protein if there is no RNA found for that gene)
  d = d[d$avg_count_pro > 0,]  # Remove 0 and lower proteins (there is a high likelihood that mRNAs WERE translated but MS didn't find them)
  d = d[nchar(d$seq_rna_cds) %% 3 == 0,] # Remove RNA's that don't have a sequence length that's a multiple of 3 (i.e. have a "broken" codon)
  # We're not going to use the mitochondrial mRNA or tRNA
  d = d[grep('^SPMIT',d$sys_id,invert=T),]
  d = d[d$avg_count_rna >= 1,]
  d$pro_per_rna = d$avg_count_pro/d$avg_count_rna
  d$codons = strsplit(d$seq_rna_cds, '(?<=...)', perl=TRUE)
  return(d)
}

AppendSorfData = function(mrnas_p){
  primes = c('5','3')
  for(p in primes){
    sorfs = read.csv(sprintf('data/mrnas_p/%sUTR_orfs.out',p), sep='\t', header=T)
    sorfs$sys_id = lapply(as.character(sorfs$GeneName), function(x){return(  substring(strsplit(x,'|',fixed=T)[[1]][1],2)  )})
    sorfs$sys_id = as.character(sorfs$sys_id)
    sorfs$GeneName = NULL
    #write.table(sorfs, file='data_processed/5UTR_orfs_cleaned.out', sep='\t',row.names=F)
    sorfs.atg_distinct = sorfs %>%
      group_by(sys_id,Frame,ATG.pos) %>%
      dplyr::summarise(avg_uORF_len = mean(uORF_len)
                       , stops = n())
    sorfs.sys_id_distinct = sorfs.atg_distinct %>%
      group_by(sys_id,Frame) %>%
      dplyr::summarise(avg_uORF_len = mean(avg_uORF_len)
                       , avg_stops_count = mean(stops)
                       , atgs = n())
    sorfs.sys_id_distinct = as.data.frame(sorfs.sys_id_distinct)
    
    mrnas_p[,sprintf('orfs_%sutr_fr0_atgs',p)] = as.numeric(apply(mrnas_p, 1, function(x) {  return(  sorfs.sys_id_distinct[sorfs.sys_id_distinct$Frame==0 & sorfs.sys_id_distinct$sys_id==x$sys_id,'atgs'])  }))
    mrnas_p[,sprintf('orfs_%sutr_fr1_atgs',p)] = as.numeric(apply(mrnas_p, 1, function(x) {  return(  sorfs.sys_id_distinct[sorfs.sys_id_distinct$Frame==-1 & sorfs.sys_id_distinct$sys_id==x$sys_id,'atgs'])  }))
    mrnas_p[,sprintf('orfs_%sutr_fr2_atgs',p)] = as.numeric(apply(mrnas_p, 1, function(x) {  return(  sorfs.sys_id_distinct[sorfs.sys_id_distinct$Frame==-2 & sorfs.sys_id_distinct$sys_id==x$sys_id,'atgs'])  }))
    mrnas_p[is.na(mrnas_p[,sprintf('orfs_%sutr_fr0_atgs',p)]),sprintf('orfs_%sutr_fr0_atgs',p)] = 0
    mrnas_p[is.na(mrnas_p[,sprintf('orfs_%sutr_fr1_atgs',p)]),sprintf('orfs_%sutr_fr1_atgs',p)] = 0
    mrnas_p[is.na(mrnas_p[,sprintf('orfs_%sutr_fr2_atgs',p)]),sprintf('orfs_%sutr_fr2_atgs',p)] = 0
  }
  return(mrnas_p)
}

AppendHalfLifeData = function(mrnas_p){
  data_hl = read.csv('data/mrnas_p/pombe_rna_halflife.csv')
  data_hl$sys_id = substring(data_hl$map.ids, 6)
  data_hl2 = data_hl[nchar(data_hl$sys_id)>0 & !grepl('gene',data_hl$sys_id),]
  data_hl2 = data_hl2[data_hl2$TU.class!='multicistronic',]
  mrnas_p$halflife = as.numeric(apply(mrnas_p, 1, function(x) {return(data_hl2[data_hl2$sys_id == x['sys_id'], 'half.life'])}))
  mrnas_p$synthesistime = as.numeric(apply(mrnas_p, 1, function(x) {return(data_hl2[data_hl2$sys_id == x['sys_id'], 'synthesis.time'])}))
  mrnas_p[is.na(mrnas_p$halflife),'halflife'] = median(mrnas_p[!is.na(mrnas_p$halflife),'halflife'])
  mrnas_p[is.na(mrnas_p$synthesistime),'synthesistime'] = median(mrnas_p[!is.na(mrnas_p$synthesistime),'synthesistime'])
  
  data_hlp = read.table('data/mrnas_p/pombe_protein_halflife.csv', sep=',', header=T, quote ="\"")
  mrnas_p$p_halflife = as.character(apply(mrnas_p, 1, function(x) {return( as.character(data_hlp[data_hlp$ENSG==x['sys_id'], 't1.2..min.'] ))}))
  tmp = mrnas_p[!mrnas_p$p_halflife %in% c('character(0)','n.d.'),'p_halflife']
  max(as.numeric(mrnas_p[!mrnas_p$p_halflife %in% c('character(0)','n.d.'),'p_halflife']))
  mrnas_p[mrnas_p$p_halflife=='character(0)','p_halflife'] = median(as.numeric(mrnas_p[!mrnas_p$p_halflife %in% c('character(0)','n.d.'),'p_halflife']))
  mrnas_p[mrnas_p$p_halflife=='n.d.','p_halflife'] = max(as.numeric(mrnas_p[!mrnas_p$p_halflife %in% c('character(0)','n.d.'),'p_halflife']))
  mrnas_p$p_halflife = as.numeric(mrnas_p$p_halflife)
  return(mrnas_p)
}

#### RPF DATA ####
Rpf._GetCoverage = function(){
  # This will take a while since the files are quite big
  d_plus = read.csv('data/mrnas_p/sub_rpf_bt2_plus.bedgraph', header=F, sep='\t', stringsAsFactors=F)
  d_minus = read.csv('data/mrnas_p/sub_rpf_bt2_minus.bedgraph', header=F, sep='\t', stringsAsFactors=F)
  d_plus_rna = read.csv('data/mrnas_p/sub_rna_bt2_plus.bedgraph', header=F, sep='\t', stringsAsFactors=F)
  d_minus_rna = read.csv('data/mrnas_p/sub_rna_bt2_minus.bedgraph', header=F, sep='\t', stringsAsFactors=F)
  
  colnames(d_plus) = colnames(d_minus) = colnames(d_plus_rna) = colnames(d_minus_rna) = c('chromosome','pos','count')

  # Split per chromosome so we can search faster later on
  Rpf.d_plus <<- Rpf.d_minus <<- Rpf.d_plus_rna <<- Rpf.d_minus_rna <<- list()
  for(c in unique(d_plus_rna$chromosome)){
    Rpf.d_plus[[c]] <<- d_plus[d_plus$count > 0 & d_plus$chromosome == c,]
    Rpf.d_minus[[c]] <<- d_minus[d_minus$count > 0 & d_minus$chromosome == c,]
    Rpf.d_plus_rna[[c]] <<- d_plus_rna[d_plus_rna$count > 0 & d_plus_rna$chromosome == c,]
    Rpf.d_minus_rna[[c]] <<- d_minus_rna[d_minus_rna$count > 0 & d_minus_rna$chromosome == c,]
  }
}

Rpf._GetExonCoords = function(mrnas_p){
  exon_cds_coords = read.table('data/mrnas_p/chromosome1.exon.coords', sep='\t',header=F,stringsAsFactors=F)
  colnames(exon_cds_coords) = c('sys_id','start','end','strand')
  exon_cds_coords$chromosome = 'I'
  temp_exon_cds_coords = read.table('data/mrnas_p/chromosome2.exon.coords', sep='\t',header=F,stringsAsFactors=F)
  colnames(temp_exon_cds_coords) = c('sys_id','start','end','strand')
  temp_exon_cds_coords$chromosome = 'II'
  exon_cds_coords = rbind(exon_cds_coords, temp_exon_cds_coords)
  temp_exon_cds_coords = read.table('data/mrnas_p/chromosome3.exon.coords', sep='\t',header=F,stringsAsFactors=F)
  colnames(temp_exon_cds_coords) = c('sys_id','start','end','strand')
  temp_exon_cds_coords$chromosome = 'III'
  exon_cds_coords = rbind(exon_cds_coords, temp_exon_cds_coords)
  
  mrnas_p$cds_exon_start_stops = apply(mrnas_p, 1, function(x){
    temp = exon_cds_coords[exon_cds_coords$sys_id == x['sys_id'],]
    return(sort(c(temp$start, temp$end)))
  })
  mrnas_p$chromosome = apply(mrnas_p, 1, function(x){
    temp = exon_cds_coords[exon_cds_coords$sys_id == x['sys_id'],]
    return(temp$chromosome[1])
  })
  mrnas_p$strand = apply(mrnas_p, 1, function(x){
    temp = exon_cds_coords[exon_cds_coords$sys_id == x['sys_id'],]
    return(ifelse(temp$strand[1]>0,'+','-'))
  })
  
  # add 1 to all exon stop values
  for(sys_id in mrnas_p[,'sys_id']){
    start_stops = unlist(mrnas_p[mrnas_p$sys_id==sys_id,'cds_exon_start_stops'])
    if(length(start_stops) > 0) {
      start_stops[seq(2,length(start_stops),2)] = start_stops[seq(2,length(start_stops),2)] + 1
      mrnas_p$cds_exon_start_stops[mrnas_p$sys_id==sys_id][1] = list(start_stops)
    }
  }
  return(mrnas_p)
}

Rpf._CalcRibosomeCov = function(mrna){
  start_ends = unlist(mrna$cds_exon_start_stops)
  
  # IN CASE WE HAVE MULTIPLE EXONS (I ASSUME THESE VALUES ARE THE SAME FOR EACH EXON OF A TRANSCRIPT)
  t_id = mrna$sys_id[1]
  t_strand = mrna$strand[1]
  t_chrom = mrna$chromosome[1]
  
  positions = seq(start_ends[1],start_ends[length(start_ends)],1)
  positions = positions[findInterval(positions, start_ends) %% 2 == 1]
  
  cov = data.frame(pos=positions)
  
  if(t_strand == '+'){
    d_plus = Rpf.d_plus[[t_chrom]]
    d_plus = d_plus[d_plus$pos >= start_ends[1] 
                                & d_plus$pos < start_ends[length(start_ends)],]
    temp = merge(x = cov, y = d_plus[d_plus$chromosome==t_chrom,c('pos','count')], by = "pos", all.x = TRUE)
    temp[is.na(temp)] = 0
    return(temp$count)
  } else {
    d_minus = Rpf.d_minus[[t_chrom]]
    d_minus = d_minus[d_minus$pos >= start_ends[1] 
                                & d_minus$pos < start_ends[length(start_ends)],]
    temp = merge(x = cov, y = d_minus[d_minus$chromosome==t_chrom,c('pos','count')], by = "pos", all.x = TRUE)
    temp[is.na(temp)] = 0
    return(rev(temp$count))
  }
}

Rpf._CalcRnaCov = function(mrna){
  start_ends = unlist(mrna$cds_exon_start_stops)
  
  # IN CASE WE HAVE MULTIPLE EXONS (I ASSUME THESE VALUES ARE THE SAME FOR EACH EXON OF A TRANSCRIPT)
  t_id = mrna$sys_id[1]
  t_strand = mrna$strand[1]
  t_chrom = mrna$chromosome[1]
  
  positions = seq(start_ends[1],start_ends[length(start_ends)],1)
  positions = positions[findInterval(positions, start_ends) %% 2 == 1]

  cov = data.frame(pos=positions)
  
  if(t_strand == '+'){
    d_plus_rna = Rpf.d_plus_rna[[t_chrom]]
    d_plus_rna = d_plus_rna[d_plus_rna$pos >= start_ends[1] 
                                & d_plus_rna$pos < start_ends[length(start_ends)],]
    cov = merge(x = cov, y = d_plus_rna[d_plus_rna$chromosome==t_chrom,c('pos','count')], by = "pos", all.x = TRUE)
    cov[is.na(cov)] = 0
    return(cov$count)
  } else {
    d_minus_rna = Rpf.d_minus_rna[[t_chrom]]
    d_minus_rna = d_minus_rna[d_minus_rna$pos >= start_ends[1] 
                                & d_minus_rna$pos < start_ends[length(start_ends)],]
    temp = merge(x = cov, y = d_minus_rna[d_minus_rna$chromosome==t_chrom,c('pos','count')], by = "pos", all.x = TRUE)
    temp[is.na(temp)] = 0
    return(rev(temp$count))
  }
}

AppendRpfData = function(mrnas_p){
  # The reason for using the coverage is that we can also compute the regression of ribosome density over the transcript (and other things in the future)
  Rpf._GetCoverage()
  mrnas_p = Rpf._GetExonCoords(mrnas_p)
  
  i=1;
  for(sys_id in mrnas_p$sys_id){
    print(sprintf('%s - %s', i, sys_id))
    mrna = mrnas_p[mrnas_p$sys_id==sys_id,]
    
    if(mrna$strand %in% c('-','+')){
      rpf_cov = Rpf._CalcRibosomeCov(mrna)
      rna_cov = Rpf._CalcRnaCov(mrna)
      cds = (mrna$seq_5utr_length+1):(GetMrnaLength(mrna)-mrna$seq_3utr_length)
      rpf_cov = rpf_cov[cds] # CDS only
      rna_cov = rna_cov[cds] # CDS only
      rna_cov[is.na(rna_cov)] = 0 # for some reason the lengths don't always match 100% so the last 2 or 3 positions can be NA
      sum_rna_cov = sum(rna_cov)
      if(sum_rna_cov == 0){
        sum_rna_cov = 32 # read length, i.e. coverage is 1 (prevents inf values later on)
      }
      m1 = lm(rpf_cov~seq(1, length(rpf_cov))) # Don't sort! We need the order as appearing on the transcript
      mrnas_p[mrnas_p$sys_id==sys_id, 'rpf_regression'] = m1$coefficients[2]
      mrnas_p[mrnas_p$sys_id==sys_id, 'rpf_per_rna'] = sum(rpf_cov)/sum_rna_cov
    }
    i = i+1
  }
  mrnas_p[is.na(mrnas_p$rpf_per_rna),'rpf_per_rna'] = median(mrnas_p[!is.na(mrnas_p$rpf_per_rna),'rpf_per_rna'])
  mrnas_p[is.na(mrnas_p$rpf_regression),'rpf_regression'] = median(mrnas_p[!is.na(mrnas_p$rpf_regression),'rpf_regression'])
  return(mrnas_p)
}

s = 'proliferation' # In case we want to do the quiescent sample as well some time
mrnas_p = read.csv(sprintf('data/mrnas_p/marguerat_rna_vs_prot_%s.csv',s))  # 5143 lines / 11 vars
mrnas_p = FilterMrnasFromMargueratData(mrnas_p)
mrnas_p = AppendSorfData(mrnas_p)
mrnas_p = AppendHalfLifeData(mrnas_p)
mrnas_p = AppendRpfData(mrnas_p) # Will take a few minutes

saveRDS(mrnas_p, file='data/mrnas_p.RObject')




