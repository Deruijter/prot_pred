require(rstudioapi)
setwd(paste(dirname(getActiveDocumentContext()$path), '/..', sep=''))

require(plyr)
require(dplyr)
library(stringi)
require(plotly) # for heatmap
require(webshot) # for exporting
library(gridExtra)
library(reshape)
library(ggplot2)
library(extrafont)
library(gdata)
library(viridis)

colfuncBPC <- colorRampPalette(c("#999999", '#115599', "black"))
#font_import(pattern="DejaVuSansMono.ttf", prompt=F) # Works native on Ubuntu but not supported on windows/IOS
font_import(pattern="LUCON.TTF", prompt=F) # Lucida Console font, might have to change this in Windows/IOS (font needs to be installed for ubuntu)
loadfonts()

#### GENERAL FUNCTIONS ####
# WIKIPEDIA GENERAL R2 FORMULA
CalcR2 = function(true, pred){ 
  RSS = sum((true-pred)^2) 
  ESS = sum(((pred-mean(true))^2))
  TSS = sum((true-mean(true))^2)
  return(1 - RSS / TSS  )
}
PrettifyFeatureNames = function(fn){
  fn = gsub('amino_acid_ratio_',replacement = 'AA ', fn)
  fn = gsub('amino_acid_prop_ratio_',replacement = 'AAP ', fn)
  fn = gsub('dg_full_14_', replacement='S14 full_', fn)
  fn = gsub('dg_utr5_14_', replacement='S14 utr5_', fn)
  fn = gsub('dg_utr3_14_', replacement='S14 utr3_', fn)
  fn = gsub('dg_cds_s2_14_', replacement='S14 cds_s2_', fn)
  fn = gsub('dg_cds_s_14_', replacement='S14 cds_s_', fn)
  fn = gsub('dg_cds_14_', replacement='S14 cds_', fn)
  fn = gsub('dg_init_14_', replacement='S14 init_', fn)
  fn = gsub('dg_full_28_', replacement='S28 full_', fn)
  fn = gsub('dg_utr5_28_', replacement='S28 utr5_', fn)
  fn = gsub('dg_utr3_28_', replacement='S28 utr3_', fn)
  fn = gsub('dg_cds_s2_28_', replacement='S28 cds_s2_', fn)
  fn = gsub('dg_cds_s_28_', replacement='S28 cds_s_', fn)
  fn = gsub('dg_cds_28_', replacement='S28 cds_', fn)
  fn = gsub('dg_init_28_', replacement='S28 init_', fn)
  fn = gsub('dg_full_50_', replacement='S50 full_', fn)
  fn = gsub('dg_utr5_50_', replacement='S50 utr5_', fn)
  fn = gsub('dg_utr3_50_', replacement='S50 utr3_', fn)
  fn = gsub('dg_cds_s2_50_', replacement='S50 cds_s2_', fn)
  fn = gsub('dg_cds_s_50_', replacement='S50 cds_s_', fn)
  fn = gsub('dg_cds_50_', replacement='S50 cds_', fn)
  fn = gsub('dg_init_50_', replacement='S50 init_', fn)
  
  fn = gsub('codon_total_ratio_',replacement = 'COD ', fn)
  fn = gsub('nt_total_ratio_',replacement = 'NUC ', fn)
  fn = gsub('nt_total_ratio_',replacement = 'NUC ', fn)
  fn = gsub('75$', replacement='75%', fn)
  fn = gsub('50$', replacement='50%', fn)
  fn = gsub('25$', replacement='25%', fn)
  fn = gsub('10$', replacement='10%', fn)
  fn = gsub('75_peaks$', replacement='75% peaks', fn)
  fn = gsub('50_peaks$', replacement='50% peaks', fn)
  fn = gsub('25_peaks$', replacement='25% peaks', fn)
  fn = gsub('10_peaks$', replacement='10% peaks', fn)
  
  fn = gsub('50_', replacement='w=50 ', fn)
  fn = gsub('28_', replacement='w=28 ', fn)
  fn = gsub('14_', replacement='w=14 ', fn)
  
  fn = gsub('q2', replacement='q3', fn)
  
  fn = gsub('seq_',replacement = 'SEQ ', fn)
  fn = gsub('cds_s2_',replacement = 'cds3 ', fn)
  fn = gsub('cds_s_',replacement = 'cds2 ', fn)
  fn = gsub('cds_',replacement = 'cds1 ', fn)
  fn = gsub('init_',replacement = 'init ', fn)
  fn = gsub('full_',replacement = 'full ', fn)
  fn = gsub('utr5_',replacement = '5utr ', fn)
  fn = gsub('utr3_',replacement = '3utr ', fn)
  
  fn = gsub('orfs_5utr_fr0_atgs',replacement = 'UO1 5utr fr0 atgs', fn)
  fn = gsub('orfs_5utr_fr1_atgs',replacement = 'UO1 5utr fr1 atgs', fn)
  fn = gsub('orfs_5utr_fr2_atgs',replacement = 'UO1 5utr fr2 atgs', fn)
  fn = gsub('orfs_3utr_fr0_atgs',replacement = 'UO1 3utr fr0 atgs', fn)
  fn = gsub('orfs_3utr_fr1_atgs',replacement = 'UO1 3utr fr1 atgs', fn)
  fn = gsub('orfs_3utr_fr2_atgs',replacement = 'UO1 3utr fr2 atgs', fn)
  fn = gsub('orfs_5utr_atgs',replacement = 'UO2 5utr atgs', fn)
  fn = gsub('orfs_3utr_atgs',replacement = 'UO2 3utr atgs', fn)
  
  fn = gsub('p_halflife',replacement = 'HLF protein halflife', fn)
  fn = gsub('^halflife$',replacement = 'HLF rna halflife', fn)
  fn = gsub('synthesistime',replacement = 'HLF rna synthesistime', fn)
  
  fn = gsub('avg_count_rna',replacement = 'ABN avg rna count', fn)
  
  fn = gsub('rpf_per_rna',replacement = 'RPF rpf/rna', fn)
  fn = gsub('rpf_regression',replacement = 'RPF rpf regression', fn)
  
  fn = gsub('pp_count',replacement = 'PPR poly proline regions', fn)
  
  fn = gsub('_', replacement=' ', fn)
  return(fn)
}
# For feature importance plot
PrettifyFeatureNamesFi = function(fn){
  fn[grep(pattern='amino_acid_ratio_.*', fn, perl=T)] = toupper(fn[grep(pattern='amino_acid_ratio_.*', fn, perl=T)])
  fn = gsub('amino_acid_ratio_',replacement = 'Amino acid ', fn, ignore.case=T)
  fn = gsub('amino_acid_prop_ratio_',replacement = 'Amino acid prop. ', fn)
  fn = gsub('dg_cds_14_', replacement='dG cds_', fn)
  fn = gsub('nt_total_ratio_',replacement = 'Nucleotide ', fn)
  
  fn = gsub('75$', replacement='75%', fn)
  fn = gsub('50$', replacement='50%', fn)
  fn = gsub('25$', replacement='25%', fn)
  fn = gsub('10$', replacement='10%', fn)
  fn = gsub('75_peaks$', replacement='75% peaks', fn)
  fn = gsub('50_peaks$', replacement='50% peaks', fn)
  fn = gsub('25_peaks$', replacement='25% peaks', fn)
  fn = gsub('10_peaks$', replacement='10% peaks', fn)
  
  fn = gsub('q2', replacement='q3', fn)
  fn = gsub('cds_',replacement = '', fn, ignore.case=T)
  
  fn = gsub('p_halflife',replacement = 'Protein half life', fn)
  fn = gsub('^halflife$',replacement = 'RNA half life', fn)
  fn = gsub('synthesistime',replacement = 'RNA synthesis time', fn)
  
  fn = gsub('avg_count_rna',replacement = 'RNA abundance', fn)
  
  fn = gsub('rpf_per_rna',replacement = 'RPF/RNA', fn)
  fn = gsub('rpf_regression',replacement = 'RPF regression', fn)
  
  fn = gsub('_', replacement=' ', fn)
  return(fn)
}

PrettifyFeatureSetNames = function(s){
  old = c('a1',
          'b1','b2','b3','b4','b5',
          'c5','c6','c7',
          'db1','db2','db3','db4','db5','db6','db7',
          'e1', 'e2',
          'f1','f2','f3','f4',
          'g5','g9','gA','gC','gD','gE','gF',
          'h1',
          'i1','i2','i3',
          'linear_model', 
          'x0_110110100', 'x0_110110100-0100','x0_110110100-0010','x0_110110100-0001','x0_110110100-1000','x0_110110100-1111')
  new = c('Sequence  -  ', 
          'Nucleotide full', 'Nucleotide 5utr', 'Nucleotide cds1', 'Nucleotide 3utr', 'Nucleotide cds2',
          'Codon cds2', 'Codon cds3', 'Codon cds1',
          'Struct W14 full', 'Struct W14 5utr', 'Struct W14 cds1', 'Struct W14 3utr', 'Struct W14 cds2', 'Struct W14 init', 'Struct W14 cds3',
          'sORF1  -  ', 'sORF2  -  ',
          'Halflife full', 'Halflife r.sr', 'Halflife r.hl', 'Halflife p.hl',
          'Proline R.  -  ', 
          'Amino Acid cds2','Amino Acid cds3',
          'A.A. prop. cds2','A.A. prop. cds3','Amino Acid cds1','A.A. prop. cds1',
          'RNA abund.  -  ',
          'Ribo-seq full', 'Ribo-seq /rna', 'Ribo-seq regr',
          'Linear     ',
          'COMP       ', 
          'COMP rna.sr',
          'COMP rna.hl',
          'COMP pro.hl',
          'COMP rpf   ',
          'COMP exp   '
  )
  
  return(new[match(s, old)])
}
# BAR PLOT COMPARISSON BETWEEN DATA SETS
BarPlotComparison = function(data_frames, var_values, var_labels, xlab, ylab, data_frame_names, title='', order_by_label=F, ylim=c(NA,NA), show_mean=F, ttl=T, xlb=T, ylb=T, lgnd=T, strlen=NA){#, colrs=plot_colors){
  data_frames = rev(data_frames)
  data_frame_names = rev(data_frame_names)
  x = unlist(data_frames[[1]][order(data_frames[[1]][,var_labels]),var_labels])
  to_plot = data.frame(x=x)
  for(i in 1:length(data_frames)){
    to_plot[,data_frame_names[[i]]] = data_frames[[i]][order(data_frames[[i]][,var_labels]),var_values]
  }
  
  order_by = 'value'
  if(order_by_label){ order_by = 'variable' }
  
  melted<-melt(to_plot, id="x")
  if(!is.na(strlen)){
    melted$x = stri_pad(melted$x, strlen, ' ', side='left')
  }
  melted$x = factor(melted$x, levels=(melted[melted$variable==tail(data_frame_names,1),])[order(melted[melted$variable==tail(data_frame_names, 1),order_by]),'x'])
  p = ggplot(melted, aes(x=x,y=value,fill=variable)) + 
    geom_bar(stat="identity",position = "dodge", width=0.8)+
    labs(title = title, x=xlab, y=ylab) +
    theme(axis.text.y = element_text(size = rel(1.5),family='Lucida Console')
          , axis.text.x = element_text(size = rel(1.5))
          , axis.title.x = element_text(size = rel(2))
          , axis.title.y = element_text(size = rel(2))
          , plot.title = element_text(size = rel(2))
          , legend.position = "top"
          , legend.title = element_blank()
          , legend.text = element_text(size = rel(1.2))
          , legend.key.size = unit(0.1, "in")
          , panel.background = element_rect(fill = 'white')
          , panel.grid.minor = element_line(color='lightgray')
          , panel.grid.major = element_line(color='lightgray')
    ) + 
    scale_y_continuous(limits = ylim) +
    coord_flip() +
    scale_fill_manual(values = colfuncBPC(length(data_frames))) +
    guides(fill = guide_legend(reverse = TRUE))
  if(show_mean){
    p = p + geom_hline(aes(yintercept=mean( (data_frames[[1]])[,var_values] ) ), color='orange', size=1)
  }
  if(!ttl){
    p = p + theme(plot.title=element_blank())
  }
  if(!xlb){
    p = p + theme(axis.title.x=element_blank())
  }
  if(!ylb){
    p = p + theme(axis.title.y=element_blank())
  }
  if(!lgnd){
    p = p + theme(legend.position="none")
  }
  return(p)
}
split_to_matrix = function(x,split){
  temp = strsplit(as.character(x),split)
  length_col = length(temp[[1]])
  return(matrix(unlist(temp), ncol=length_col, byrow=T))
}
calc_rpkm = function(d){
  return(d[,'reads'] / ((d[,'total_exon_length']/1000) * (sum(d$reads)/1000000)))
}


# #### ML SCORES FILE ####
# # Read ml scores file and make the output ready for analysis in R
# ml_scores = read.table('ml_output/ml_scores.csv',header=F,sep=',', stringsAsFactors=F)
# colnames(ml_scores) = c('learner','data_set','testset','training_score','test_score','fi')
# ml_scores$base_learner = apply(ml_scores, 1, function(x) return(substring(x['learner'],1,3)))
# ml_scores$learner_params = apply(ml_scores, 1, function(x) return(substring(x['learner'],5)))
# 
# # Extract only best training results
# ml_scores_filtered = as.data.frame(ml_scores %>%
#                      group_by(base_learner, data_set, testset) %>%
#                      summarize(max_train = max(training_score)))
# 
# # removed feature sets:
# removed_feature_sets = c('a1_p','c1','c2','c3','c4',
#                          'da1','da2','da3','da4','da5','da6','da7',
#                          'dc1','dc2','dc3','dc4','dc5','dc6','dc7',
#                          'dd1','dd2','dd3','dd4','dd5','dd6','dd7',
#                          'g1','g2','g3','g4','g6','g7','g8','gB',
#                          'knn','nnr','rfr','rir', 'rir_p',
#                          'rfr1','rfr2','rfr3','rir1','rir2','rir3',
#                          'ultra_mega','ultra_mega_rfe_et','ultra_mega_rfe_ls','ultra_mega_rfe_rr',
#                          'ultra_mega_sfm_en0.25_0.05','ultra_mega_sfm_en0.70_0.05','ultra_mega_sfm_en0.95_0.05',
#                          'ultra_mega_sfm_et','ultra_mega_sfm_ls','ultra_mega_var_vt0.05',
#                          'x0_111111111', 'x0_111110100','x0_111100100','x0_111010100','x0_111000100','x0_110100100')
# 
# ml_scores_filtered = ml_scores_filtered[!ml_scores_filtered$data_set %in% removed_feature_sets,]
# 
# for(i in 1:nrow(ml_scores_filtered)){
#   print(i)
#   ml_scores_filtered[i, 'test_score'] = mean(ml_scores[ml_scores$base_learner==ml_scores_filtered[i,'base_learner'] & ml_scores$data_set==ml_scores_filtered[i,'data_set'] & ml_scores$testset==ml_scores_filtered[i,'testset'] & ml_scores$training_score==ml_scores_filtered[i,'max_train'] & ml_scores$test_score!=-1, 'test_score'] )
#   ml_scores_filtered[i, 'fi'] = ml_scores[ml_scores$base_learner==ml_scores_filtered[i,'base_learner'] & ml_scores$data_set==ml_scores_filtered[i,'data_set'] & ml_scores$testset==ml_scores_filtered[i,'testset'] & ml_scores$training_score==ml_scores_filtered[i,'max_train'] & ml_scores$test_score!=-1, 'fi'][[1]]
# }
# ml_scores_filtered$fi = c(apply(ml_scores_filtered, 1, function(x) return(  as.numeric(unlist(stri_extract_all_regex(x['fi'], "-?[0-9]\\.[0-9]*(e-[0-9]*)?")))  ))) # Make feature importance readible by R
# 
# 
# 
# 
# 
# 
# # Read ml scores file and make the output ready for analysis in R
# ml_scores = read.table('ml_output/ml_scores.csv',header=F,sep=',', stringsAsFactors=F)
# colnames(ml_scores) = c('learner','data_set','testset','training_score','test_score','fi')
# ml_scores$base_learner = apply(ml_scores, 1, function(x) return(substring(x['learner'],1,3)))
# ml_scores$learner_params = apply(ml_scores, 1, function(x) return(substring(x['learner'],5)))
# 
# # Extract only best training results
# ml_scores_filtered2 = as.data.frame(ml_scores %>%
#                                      group_by(base_learner, data_set, testset) %>%
#                                      summarize(max_train = max(training_score)))
# 
# # removed feature sets:
# removed_feature_sets = c('a1_p','c1','c2','c3','c4',
#                          'da1','da2','da3','da4','da5','da6','da7',
#                          'dc1','dc2','dc3','dc4','dc5','dc6','dc7',
#                          'dd1','dd2','dd3','dd4','dd5','dd6','dd7',
#                          'g1','g2','g3','g4','g6','g7','g8','gB',
#                          'knn','nnr','rfr','rir', 'rir_p',
#                          'rfr1','rfr2','rfr3','rir1','rir2','rir3',
#                          'ultra_mega','ultra_mega_rfe_et','ultra_mega_rfe_ls','ultra_mega_rfe_rr',
#                          'ultra_mega_sfm_en0.25_0.05','ultra_mega_sfm_en0.70_0.05','ultra_mega_sfm_en0.95_0.05',
#                          'ultra_mega_sfm_et','ultra_mega_sfm_ls','ultra_mega_var_vt0.05',
#                          'x0_111111111', 'x0_111110100','x0_111100100','x0_111010100','x0_111000100','x0_110100100',
#                          'best1','best2','best3','rfr_b1','rfr_b2','rfr_b3','rir_b1','rir_b2','rir_b3',
#                          'ultra_mega2_rfe_et','ultra_mega2_rfe_ls','ultra_mega2_rfe_rr','ultra_mega2_sfm_en0.25_0.05',
#                          'ultra_mega2_sfm_en0.70_0.05','ultra_mega2_sfm_en0.95_0.05','ultra_mega2_sfm_et',
#                          'ultra_mega2_sfm_ls','ultra_mega2_var_vt0.05', 'ultra_mega2')
# 
# ml_scores_filtered2 = ml_scores_filtered2[!ml_scores_filtered2$data_set %in% removed_feature_sets,]
# 
# for(i in 1:nrow(ml_scores_filtered2)){
#   ml_scores_filtered2[i, 'test_score'] = mean(ml_scores[ml_scores$base_learner==ml_scores_filtered2[i,'base_learner'] & ml_scores$data_set==ml_scores_filtered2[i,'data_set'] & ml_scores$testset==ml_scores_filtered2[i,'testset'] & ml_scores$training_score==ml_scores_filtered2[i,'max_train'] & ml_scores$test_score!=-1, 'test_score'] )
#   ml_scores_filtered2[i, 'fi'] = ml_scores[ml_scores$base_learner==ml_scores_filtered2[i,'base_learner'] & ml_scores$data_set==ml_scores_filtered2[i,'data_set'] & ml_scores$testset==ml_scores_filtered2[i,'testset'] & ml_scores$training_score==ml_scores_filtered2[i,'max_train'] & ml_scores$test_score!=-1, 'fi'][[1]]
# }
# ml_scores_filtered2$fi = c(apply(ml_scores_filtered2, 1, function(x) return(  as.numeric(unlist(stri_extract_all_regex(x['fi'], "-?[0-9]\\.[0-9]*(e-[0-9]*)?")))  ))) # Make feature importance readible by R
# 
# ml_scrs = rbind(ml_scores_filtered2, ml_scores_filtered)
# ml_scores_filtered = ml_scrs
# saveRDS(ml_scores_filtered, file = 'data/ml_scores_filtered2.RObject')
ml_scores_filtered = readRDS(file = 'data/ml_scores_filtered2.RObject')



# Fig XXX and Fig XXX
CreateHeatmap = function(){
  scr = readRDS('ml_scores_filtered2.RObject')
  scr = scr[!grepl('\\_p[0-9]{1,2}', scr$data_set),] # Remove all results from the permutation test
  scr[scr$test_score < 0, 'test_score'] = 0 # Set the cut off for the plot at 0 (so as good as, or worse than random)
  
  learner_data_set_means = as.data.frame(scr %>%
                                           group_by(base_learner, data_set) %>%
                                           summarize(mean_score = mean(test_score)
                                                     , stdev = sd(test_score)
                                                     , count = n()))
  learner_data_set_means$data_set = PrettifyFeatureSetNames(learner_data_set_means$data_set)
  mat = cbind(
    learner_data_set_means[learner_data_set_means$base_learner=='knn','mean_score'],
    learner_data_set_means[learner_data_set_means$base_learner=='rir','mean_score'],
    learner_data_set_means[learner_data_set_means$base_learner=='svm','mean_score'], 
    learner_data_set_means[learner_data_set_means$base_learner=='nnr','mean_score'],
    learner_data_set_means[learner_data_set_means$base_learner=='rfr','mean_score'])
  
  # Reorder so the heatmap is a bit more organized
  order = c(1, #seq
            2,4,6,3,5, # nuc
            9,8,7, #cod
            28,25,24, # aa
            29,27,26, # aap
            23, # poly proline
            31,32,33, # rpf
            10,12,16,14,11,13,15, # struct
            19,22,21,20, # halflife
            30, # rna abdnc
            17,18, # sorf
            34,35,36,37,38,39 # comp
  )
  mat = mat[rev(order),]
  
  
  data_set_names = rev(as.list(unique(learner_data_set_means$data_set))[order])
  
  # Make the data set names a bit nicer
  # data_set_names = list('Sequence  -  ',
  #                       'Nucleotide full', 'Nucleotide 5utr', 'Nucleotide cds1', 'Nucleotide 3utr', 'Nucleotide cds2',
  #                       'Manual1 codon', 'Manual1 amino', 'Manual1 c & a',
  #                       'Codon cds2', 'Codon cds3', 'Codon cds1',
  #                       'Struct W14 full', 'Struct W14 5utr', 'Struct W14 cds1', 'Struct W14 3utr', 'Struct W14 cds2', 'Struct W14 init', 'Struct W14 cds3',
  #                       'sORF1  -  ', 'sORF2  -  ',
  #                       'Halflife  -  ',
  #                       'Proline R.  -  ',
  #                       'Amino Acid cds2','Amino Acid cds3',
  #                       'A.A. prop. cds2','A.A. prop. cds3','Amino Acid cds1','A.A. prop. cds1',
  #                       'RNA abund.  -  ','Ribo-seq cds1',
  #                       'Manual2 codon','Manual2 amino','Manual2 c & a',
  #                       'Manual3 codon','Manual3 amino','Manual3 c & a',
  #                       'ALL DATA',
  #                       'RFE ExtraTrees', 'RFE LinearSVR ', 'RFE RidgeRegr.',
  #                       'SFM EN L1=0.25','SFM EN L1=0.70','SFM EN L1=0.95','SFM ExtraTrees','SFM LinearSVR ',
  #                       'Var. thresh. 5%'
  # )
  
  # STYLE
  f <- list(
    family = "Lucida Console",
    size = 12
  )
  x = list(
    title = "Machine Learners",
    titlefont = list(size=16, color='black')
  )
  y <- list(
    title = "Feature Sets",
    titlefont = list(size=16, color='black')
  )
  colfunc <- colorRampPalette(c("#FFFFFF", "#000000"))
  
  
  
  p_all = plot_ly(x=list('KNN','RIR','SVM','NNR','RFR'), z = mat, type = "heatmap", 
                  y=data_set_names,
                  colors=colfunc(100)[round(100 * min(mat) / max(mat)):round((max(unlist(mat)) / max(unlist(mat))) * 100)],
                  width=400, height=(20*nrow(mat))+20, size=0.3) %>%
    layout(margin=list(b=20, l=140, r=0, t=0, pad=0, autoexpand=T), font=f)
  print(p_all)
  export(p_all, file="figures/heatmap_full.pdf", vwidth=400, vheight=(20*nrow(mat))+20, delay=0)
  
  
  
  comp_only = seq(1,6)
  m = mat[comp_only,]
  d = data_set_names[comp_only]
  p_all = plot_ly(x=list('KNN','RIR','SVM','NNR','RFR'), z = m, type = "heatmap", 
                  y=d,
                  colors=colfunc(100)[round(100 * min(m) / max(m)):round((max(unlist(m)) / max(unlist(mat))) * 100)],
                  width=400, height=(20*nrow(m))+20, size=0.3) %>%
    layout(margin=list(b=20, l=140, r=0, t=0, pad=0, autoexpand=T), font=f)
  print(p_all)
  export(p_all, file="figures/heatmap_comp.pdf", vwidth=400, vheight=(20*nrow(m))+20, delay=0)
  
  
  
  # COMPLETE HEAT MAP
  individual_only = c(7,8,9,
                      10,11,12,13,
                      19,20,
                      21,22,23,
                      24,
                      27,
                      30,
                      33,
                      37,38,
                      39
                      )
  #d = rev(list('COMP         ','COMP + rna.st','COMP + rna.hl', 'COMP + pro.hl', 'COMP + rpf   ', 'COMP + exp   '))
  d = data_set_names[individual_only]
  m = mat[individual_only,]
  p_all = plot_ly(x=list('KNN','RIR','SVM','NNR','RFR'), z = m, type = "heatmap", 
                  y=d,
                  colors=colfunc(100)[round(100 * min(m) / max(m)):round((max(unlist(m)) / max(unlist(mat))) * 100)],
                  width=400, height=(20*nrow(m))+20, size=0.3) %>%
    layout(margin=list(b=20, l=140, r=0, t=0, pad=0, autoexpand=T), font=f)
  print(p_all)
  export(p_all, file="figures/heatmap_individual_filtered.pdf", vwidth=400, vheight=(20*nrow(m))+20, delay=0)
}

# Fig XXX 
GetRpRatioPredictions = function() {
  learner= 'rfr'
  
  predictions = data.frame()
  training_ids = read.table('./data/training_ids.csv', sep='\t', header=T, stringsAsFactors=F)
  mrnas_p = readRDS('./data/mrnas_p.RObject')
  feature_sets = c('linear_model', 'x0_110110100', 'x0_110110100-0100','x0_110110100-0010','x0_110110100-0001','x0_110110100-1000','x0_110110100-1111')
  for(fs in feature_sets) {
    pred_fs = data.frame()
    #fs='linear_model'
    for(i in 1:10){
      train_ids = training_ids[training_ids[,sprintf('set_%s',i)]==T,'sys_id']
      test_ids = training_ids[training_ids[,sprintf('set_%s',i)]==F,'sys_id']
      
      if(fs=='linear_model'){
        mrnas_train = mrnas_p[mrnas_p$sys_id %in% train_ids,]
        mrnas_test = mrnas_p[mrnas_p$sys_id %in% test_ids,]
        
        m2 = lm(mrnas_train$avg_count_pro ~ mrnas_train$avg_count_rna)
        
        real_values = mrnas_test$pro_per_rna
        pred_values = rep(m2$coefficients[2], length(real_values))
        pred_temp = data.frame(mrnas_test$sys_id, real_values, pred_values)
        
        pred_temp[,4] = apply(pred_temp, 1, function(x) {  return(mrnas_p[mrnas_p$sys_id==x[1], 'avg_count_pro'])   }  )
        pred_temp[,5] = apply(pred_temp, 1, function(x) {  return(mrnas_p[mrnas_p$sys_id==x[1], 'avg_count_rna'])   }  )
        pred_temp[,6] = pred_temp[,5] * pred_temp[,3]
        
        pred_fs = rbind(pred_fs, pred_temp)
        
      } else {
        pred_temp = read.table(sprintf('./ml_output/predictions/ml_predictions_%s_%s_%s.csv', learner,fs,i),header=F,sep=',', stringsAsFactors=F)
        pred_temp = cbind(as.data.frame(test_ids[1:951]), pred_temp, stringsAsFactors=F) # ORDER STAYED THE SAME
        
        # DE-SCALE THE DATA
        descaled_real = (pred_temp[,2] * sd(mrnas_p[mrnas_p$sys_id %in% train_ids, 'pro_per_rna'])) + mean(mrnas_p[mrnas_p$sys_id %in% train_ids, 'pro_per_rna'])
        m = lm(mrnas_p[mrnas_p$sys_id %in% test_ids[1:951], 'pro_per_rna'] ~ descaled_real)
        descaled_real = descaled_real * m$coefficients[2] + m$coefficients[1] # should match the real rp_ratio in the mrnas_p data set
        pred_temp[,2] = descaled_real
        
        descaled_pred = (pred_temp[,3] * sd(mrnas_p[mrnas_p$sys_id %in% train_ids, 'pro_per_rna'])) + mean(mrnas_p[mrnas_p$sys_id %in% train_ids, 'pro_per_rna'])
        descaled_pred = descaled_pred * m$coefficients[2] + m$coefficients[1] # De-scale using only the training data!
        pred_temp[,3] = as.numeric(as.character(descaled_pred))
        
        pred_temp[,4] = apply(pred_temp, 1, function(x) {  return(mrnas_p[mrnas_p$sys_id==x[1], 'avg_count_pro'])   }  )
        pred_temp[,5] = apply(pred_temp, 1, function(x) {  return(mrnas_p[mrnas_p$sys_id==x[1], 'avg_count_rna'])   }  )
        pred_temp[,6] = pred_temp[,5] * pred_temp[,3]
        
        pred_fs = rbind(pred_fs, pred_temp)
      }
    }
    
    colnames(pred_fs) = c('sys_id', 'real','pred', 'avg_count_pro', 'avg_count_rna', 'pred_pro')
    
    pred_fs = as.data.frame(pred_fs %>%
                              group_by(sys_id) %>%
                              summarize(real = mean(real)
                                        , pred = mean(pred)
                                        , avg_count_pro = mean(avg_count_pro)
                                        , avg_count_rna = mean(avg_count_rna)
                                        , pred_pro = mean(pred_pro)
                                        ))
    
    pred_fs[,'log_real'] = log(pred_fs$real)
    pred_fs[,'log_pred'] = log(pred_fs$pred)
    pred_fs[,'regression'] = predict(lm(pred_fs$log_real ~ pred_fs$log_pred))
    pred_fs[,'d'] = densCols(pred_fs$log_real, pred_fs$log_pred, colramp = colorRampPalette(viridis(10)))
    
    pred_fs[,'col'] = sprintf('%s | r² %.2f',PrettifyFeatureSetNames(fs), round(CalcR2(pred_fs[,2], pred_fs[,3]),2))
    predictions = rbind(predictions, pred_fs)
  }
  
  predictions$col = factor(predictions$col, levels=sort(unique(predictions$col))[c(7,3,5,4,2,6,1)])
  predictions = predictions[order(predictions$regression),]
  return(predictions)
}
GetProteinPredictions = function() {
  learner= 'rfr'
  
  predictions = data.frame()
  training_ids = read.table('./data/training_ids.csv', sep='\t', header=T, stringsAsFactors=F)
  mrnas_p = readRDS('./data/mrnas_p.RObject')
  feature_sets = c('linear_model', 'x0_110110100', 'x0_110110100-0100','x0_110110100-0010','x0_110110100-0001','x0_110110100-1000','x0_110110100-1111')
  for(fs in feature_sets) {
    pred_fs = data.frame()
    #fs='linear_model'
    for(i in 1:10){
      train_ids = training_ids[training_ids[,sprintf('set_%s',i)]==T,'sys_id']
      test_ids = training_ids[training_ids[,sprintf('set_%s',i)]==F,'sys_id']
      
      if(fs=='linear_model'){
        mrnas_train = mrnas_p[mrnas_p$sys_id %in% train_ids,]
        mrnas_test = mrnas_p[mrnas_p$sys_id %in% test_ids,]
        
        m2 = lm(mrnas_train$avg_count_pro ~ mrnas_train$avg_count_rna)
        
        real_values = mrnas_test$pro_per_rna
        pred_values = rep(m2$coefficients[2], length(real_values))
        pred_temp = data.frame(mrnas_test$sys_id, real_values, pred_values)
        
        pred_temp[,4] = apply(pred_temp, 1, function(x) {  return(mrnas_p[mrnas_p$sys_id==x[1], 'avg_count_pro'])   }  )
        pred_temp[,5] = apply(pred_temp, 1, function(x) {  return(mrnas_p[mrnas_p$sys_id==x[1], 'avg_count_rna'])   }  )
        pred_temp[,6] = pred_temp[,5] * pred_temp[,3]
        
        pred_temp[,7] = '#FF5577'
        pred_temp[,8] = '#FF5577'
        
        pred_fs = rbind(pred_fs, pred_temp)
        
      } else {
        pred_temp = read.table(sprintf('./ml_output/ml_predictions_%s_%s_%s.csv', learner,fs,i),header=F,sep=',', stringsAsFactors=F)
        pred_temp = cbind(as.data.frame(test_ids[1:951]), pred_temp, stringsAsFactors=F) # ORDER STAYED THE SAME
        
        # DE-SCALE THE DATA
        descaled_real = (pred_temp[,2] * sd(mrnas_p[mrnas_p$sys_id %in% train_ids, 'pro_per_rna'])) + mean(mrnas_p[mrnas_p$sys_id %in% train_ids, 'pro_per_rna'])
        m = lm(mrnas_p[mrnas_p$sys_id %in% test_ids[1:951], 'pro_per_rna'] ~ descaled_real)
        descaled_real = descaled_real * m$coefficients[2] + m$coefficients[1] # should match the real rp_ratio in the mrnas_p data set
        pred_temp[,2] = descaled_real
        
        mrnas_test = mrnas_p[mrnas_p$sys_id %in% test_ids,]
        pred_temp[,2] = mrnas_test$pro_per_rna[1:951]
        
        descaled_pred = (pred_temp[,3] * sd(mrnas_p[mrnas_p$sys_id %in% train_ids, 'pro_per_rna'])) + mean(mrnas_p[mrnas_p$sys_id %in% train_ids, 'pro_per_rna'])
        descaled_pred = descaled_pred * m$coefficients[2] + m$coefficients[1] # De-scale using only the training data!
        pred_temp[,3] = as.numeric(as.character(descaled_pred))
        
        pred_temp[,4] = apply(pred_temp, 1, function(x) {  return(mrnas_p[mrnas_p$sys_id==x[1], 'avg_count_pro'])   }  )
        pred_temp[,5] = apply(pred_temp, 1, function(x) {  return(mrnas_p[mrnas_p$sys_id==x[1], 'avg_count_rna'])   }  )
        pred_temp[,6] = pred_temp[,5] * pred_temp[,3]
        
        pred_temp[,7] = '#006666'
        pred_temp[,8] = '#006666'
        
        pred_fs = rbind(pred_fs, pred_temp)
      }
    }
    
    colnames(pred_fs) = c('sys_id', 'real','pred', 'avg_count_pro', 'avg_count_rna', 'pred_pro', 'clr', 'clr2')
    
    pred_fs = as.data.frame(pred_fs %>%
                              group_by(sys_id, real, avg_count_pro, avg_count_rna, clr, clr2) %>%
                              summarize(#real = mean(real)
                                pred = mean(pred)
                                #, avg_count_pro = mean(avg_count_pro)
                                #, avg_count_rna = mean(avg_count_rna)
                                , pred_pro = mean(pred_pro)
                              ))
    
    pred_fs[,'log_real'] = log(pred_fs$real)
    pred_fs[,'log_pred'] = log(pred_fs$pred)
    pred_fs[,'log_avg_count_pro'] = log(pred_fs$avg_count_pro)
    pred_fs[,'log_pred_pro'] = log(pred_fs$pred_pro)
    pred_fs[,'regression'] = predict(lm(pred_fs$log_real ~ pred_fs$log_pred))
    pred_fs[,'d'] = densCols(pred_fs$log_real, pred_fs$log_pred, colramp = colorRampPalette(viridis(10)))
    
    pred_fs[,'dens'] = densCols(pred_fs$log_avg_count_pro, pred_fs$log_pred_pro, colramp = colorRampPalette(viridis(10)))
    
    pred_fs[,'col'] = sprintf('%s | r² %.2f',PrettifyFeatureSetNames(fs), round(CalcR2(pred_fs$log_avg_count_pro, pred_fs$log_pred_pro),2))
    predictions = rbind(predictions, pred_fs)
  }
  
  predictions$col = factor(predictions$col, levels=sort(unique(predictions$col))[c(7,3,5,4,2,6,1)])
  predictions = predictions[order(predictions$regression),]
  return(predictions)
}
CreateRpRatioPrediction = function(){
  # Predictions for the best learner / feature set combination
  predictions = GetRpRatioPredictions()
  
  p = ggplot(predictions) + # export size 6 x 4 inch
    geom_point(aes(log_pred, log_real, color=d), size = 2, show.legend = T) +
    scale_color_identity() +
    geom_abline(aes(intercept=0, slope=1, color='#000000AA'), size=1.5, show.legend = T) + 
    facet_wrap(~ col, ncol=7) + 
    theme_bw() + 
    theme(legend.position="top"
          , title = element_text(size = rel(1.5))
          , legend.text = element_text(size = rel(1.4))
          , axis.text = element_text(size=rel(1.4))
          , strip.text = element_text(size=rel(1.4))
    ) + 
    labs(title='RP ratio predictions (log scale)', x='Predicted RP ratio', y='Real RP ratio')
  
  pdf(file = "figures/rpratio_predictions.pdf", 17, 4)
  print(p)
  dev.off()
  print(p)
}

# Fig XXX
CreateProteinAbundancePredictions = function(){
  # training_ids = read.table('./data/training_ids2.csv', sep='\t', header=T, stringsAsFactors=F)
  # mrnas_p = readRDS('./data/mrnas_p.RObject')
  # 
  # # Get predictions made for the best learner/feature set combination
  # learner= 'rfr'
  # ds='x0_111111111'
  # ds='x0_110110100'
  # #ds='rir_b2'
  # predictions_detailed = mrnas_p[,c('sys_id','avg_count_rna','avg_count_pro','pro_per_rna')]
  # predictions_detailed$pro_per_rna_pred = list(c())
  # predictions_detailed$avg_count_pro_pred = list(c())
  # for(i in c(1:10)){
  #   test_ids = training_ids[training_ids[,sprintf('set_%s',i)]==F,'sys_id']
  #   train_ids = training_ids[training_ids[,sprintf('set_%s',i)]==T,'sys_id']
  #   predictions = read.table(sprintf('./ml_output/ml_predictions_%s_%s_%s.csv', learner,ds,i),header=F,sep=',', stringsAsFactors=F)
  #   predictions = cbind(as.data.frame(test_ids[1:951]), predictions, stringsAsFactors=F) # ORDER STAYED THE SAME
  #   
  #   # DE-SCALE THE DATA
  #   descaled_real = (predictions[,2] * sd(mrnas_p[mrnas_p$sys_id %in% train_ids, 'pro_per_rna'])) + mean(mrnas_p[mrnas_p$sys_id %in% train_ids, 'pro_per_rna'])
  #   m = lm(mrnas_p[mrnas_p$sys_id %in% test_ids[1:951], 'pro_per_rna'] ~ descaled_real)
  #   descaled_real = descaled_real * m$coefficients[2] + m$coefficients[1] # should match the real rp_ratio in the mrnas_p data set
  #   predictions[,2] = descaled_real
  #   
  #   descaled_pred = (predictions[,3] * sd(mrnas_p[mrnas_p$sys_id %in% train_ids, 'pro_per_rna'])) + mean(mrnas_p[mrnas_p$sys_id %in% train_ids, 'pro_per_rna'])
  #   descaled_pred = descaled_pred * m$coefficients[2] + m$coefficients[1] # De-scale using only the training data!
  #   predictions[,3] = as.numeric(as.character(descaled_pred))
  #   
  #   apply(predictions, 1, function(x) {
  #     # For future reference: http://stackoverflow.com/questions/42907524/overwrite-append-to-multiple-value-column-in-r-data-frame
  #     predictions_detailed$pro_per_rna_pred[predictions_detailed$sys_id==x[1]] <<- list(c(unlist(predictions_detailed[predictions_detailed$sys_id==x[1],]$pro_per_rna_pred), as.numeric(x[3])))
  #     predictions_detailed$avg_count_pro_pred[predictions_detailed$sys_id==x[1]] <<- list(c(unlist(predictions_detailed[predictions_detailed$sys_id==x[1],]$avg_count_pro_pred), (predictions_detailed[predictions_detailed$sys_id==x[1],]$avg_count_rna * as.numeric(x[3]))))
  #   })
  # }
  # 
  # predictions_detailed$avg_count_pro_pred_avg = apply(predictions_detailed, 1, function(x){  
  #   mean(as.numeric(unlist(x['avg_count_pro_pred'],',')))
  # })
  # predictions_detailed = predictions_detailed[!is.na(predictions_detailed$avg_count_pro_pred_avg),] # Some RNA's did not have predictions
  
  predictions = GetProteinPredictions()
  
  contour = predictions
  contour$plottype = 'Real conc.'
  density = predictions
  density$plottype = 'Density'
  contour$plottype = factor(contour$plottype, levels=c('Real conc.','Density'))
  density$plottype = factor(density$plottype, levels=c('Real conc.','Density'))
  
  p <- ggplot(data=contour) + 
    facet_grid(plottype ~ col, scales='free', as.table=T) +
    geom_point(data=contour, aes(log_pred_pro, log_avg_count_pro), colour='#DDDDDD', size = 2, show.legend = T) +
    scale_color_identity() +
    scale_fill_identity() +
    scale_alpha_identity() +
    geom_point(data=contour, aes(log_avg_count_pro, log_avg_count_pro), show.legend=F) + 
    stat_density2d(data=contour, aes(log_pred_pro, log_avg_count_pro, colour=clr), bins=10, size=0.7, show.legend=F) +
    geom_density(data=density, aes(log_avg_count_pro), fill='gray', colour='black', size=0.7) +
    geom_density(data=density, aes(log_pred_pro, colour=clr, fill=clr, alpha=0.3), size=0.7, show.legend=F) +
    labs(x='Protein abundance (log scale)', y='',title='Protein abundance predictions') +
    theme_bw() + 
    theme(legend.position="top"
          , title = element_text(size = rel(1.5))
          , legend.text = element_text(size = rel(1.4))
          , axis.text = element_text(size=rel(1.4))
          , strip.text = element_text(size=rel(1.4))
    )
  
  pdf(file = "figures/protein_abund_predictions.pdf", 17, 6)
  print(p)
  dev.off()
  print(p)
  
  return()
  
  
  #### OLD !!
  
  #linear_model = predict(lm(predictions_detailed$avg_count_pro ~ predictions_detailed$avg_count_rna))
  
  # DENSITY PLOT DATA
  counts_real = data.frame(concentration=log(predictions_detailed$avg_count_pro))
  counts_real$type = '3. perfect pred.'
  counts_pred_lin = data.frame(concentration=log(linear_model))
  counts_pred_lin$type = '2. linear model '
  counts_pred_adv = data.frame(concentration=log(predictions_detailed$avg_count_pro_pred_avg))
  counts_pred_adv$type = '1. RFR model '
  counts = rbind(counts_real, counts_pred_lin, counts_pred_adv)
  counts$plottype='density'
  
  # CONTOUR PLOT DATA
  df1 = data.frame(x=log(predictions_detailed$avg_count_pro_pred_avg), y = log(predictions_detailed$avg_count_pro))
  df2 = data.frame(x=log(linear_model), y = log(predictions_detailed$avg_count_pro))
  df3 = data.frame(x=log(predictions_detailed$avg_count_pro), y = log(predictions_detailed$avg_count_pro))
  df1$type='1. RFR model '
  df2$type='2. linear model '
  df3$type='3. perfect pred.'
  dfa = rbind(df1,df2)
  dfb = df3
  dfa$plottype = dfb$plottype = 'real conc.'
    
  r2_linear = CalcR2(log(predictions_detailed$avg_count_pro), log(linear_model))
  r2_rfr = CalcR2(log(predictions_detailed$avg_count_pro), log(predictions_detailed$avg_count_pro_pred_avg))
  txt = data.frame(x=4, y= 13.5, type='1. RFR model ', text=sprintf('R² %.2f',round(r2_rfr,2)), plottype='real conc.')
  txt = rbind(txt, data.frame(x=6, y= 13.5, type='2. linear model ', text=sprintf('R² %.2f',round(r2_linear,2)), plottype='real conc.'))
  p <- ggplot(data=dfa) + # export size: 6 x 6 inch
    facet_grid(plottype~., scales='free', as.table=F) +
    geom_point(data=dfb, aes(x,y,color=type), show.legend=F) + 
    stat_density2d(data=dfa, aes(x,y,color=type), size=0.7, show.legend=F) + 
    geom_text(data=txt, aes(x,y,color=type, label=text), size=5, show.legend=F) + 
    geom_density(data=counts, aes(concentration, fill = type, colour=type), size=0.7) +
    labs(x='Protein concentration (log scale)', y='',title='Protein concentrations predictions') +
    scale_colour_manual("", values = c('#006666','#FF5577','#000000')) +
    scale_fill_manual("", values = c('#00BBBB','#FF557777','#00000055')) +#00BBBB
    theme(legend.position="top"
          , title = element_text(size = rel(1.5))
          , legend.text = element_text(size = rel(1.4))
          , axis.text = element_text(size=rel(1.4))
          , strip.text = element_text(size=rel(1.4))
    )
  print(p)
}

# Fig XXX and Fig XXX
CreateFeatureVsRPRatio = function(){
  
  data_set_ultra_mega = read.table('./feature_sets/data_set_x0_110110100-1111.csv', sep='\t', header = T)
  features = c('seq_prot_length','rpf_per_rna','halflife','p_halflife',
               'amino_acid_ratio_cds_r', 'amino_acid_ratio_cds_g','amino_acid_ratio_cds_s','amino_acid_ratio_cds_c',
               'codon_total_ratio_cds_agg', 'codon_total_ratio_cds_ggt','codon_total_ratio_cds_tca','codon_total_ratio_cds_tgt')
  titles = c('Protein length','RPF/RNA','RNA half-life','Prot. half-life',
             'Arginine','Glycine','Serine','Cysteine',
             'AGG (R)', 'GGT (G)', 'TCA (S)', 'TGT (C)')
  
  features = c('rpf_per_rna','rpf_regression','amino_acid_ratio_cds_r',
               'halflife', 'amino_acid_prop_ratio_cds_n','dg_cds_14_q1')
  titles = c('RPF/RNA','RPF regression','Arginine',
             'RNA half life','Non polar A.A.','Q1 delta G')
  cols = 1:3
  rows = 1:2
  feature_data = data.frame()
  for(row in rows){
    for(col in cols){
      i = length(cols)*(row-1) + col
      f = features[i]
      title = titles[i]
      f1 = log(data_set_ultra_mega$pro_per_rna)
      f2 = data_set_ultra_mega[,f]
      f2 = f2 - min(c(data_set_ultra_mega[data_set_ultra_mega[,f]<0,f],0))
      f2 = f2 + min(f2[f2>0]/2)
      f2 = scale(log(f2))
      f2[!findInterval(f2,c(-10,10))==1] = 0
      c = cor.test(f1, f2, method='spearman')
      title = sprintf('%s\nrho: %s, p: %s', title, round(c$estimate,2), ifelse(c$p.value < 0.01, '<0.01', as.character(round(c$p.value,2))))
      d = densCols(f1, f2, colramp = colorRampPalette(viridis(10))) # colorRampPalette(rev(rainbow(10, end = 4/6))))
      temp = data.frame(x=f1, y=f2, col=title, d = d)
      feature_data = rbind(feature_data, temp)
    }
  }
  feature_data.cor <- ddply(feature_data, .(col), function(val) sprintf("rho=%.2f | p=%s", cor(val$x, val$y, method='spearman'), ifelse(cor.test(val$x, val$y, method='spearman')$p.value<0.01,'<0.01',round(cor.test(val$x, val$y, method='spearman')$p.value,2))))
  
  p <- ggplot(feature_data) + # Export 6x6 inch / 600x600 px
    geom_point(aes(x, y, col = d), size = 1) +
    scale_color_identity() +
    scale_alpha(range = c(0.1, 1)) +
    facet_wrap(~ col, ncol=3) + 
    theme_bw() + 
    theme(legend.position="top"
          , title = element_text(size = rel(1.5))
          , legend.text = element_text(size = rel(1.4))
          , axis.text = element_text(size=rel(1.4))
          , strip.text = element_text(size=rel(0.9))
    ) + 
    labs(x='RP ratio (log scaled)', y='Feature value (centered and scaled)') 
  
  print(p) 
  pdf(file = "figures/feature_vs_rpratio_6.pdf", 6, 4)
  print(p) 
  dev.off()
  
  
  # Another plot for only the most important amino acids (according to the feature importance plot)
  features = c('amino_acid_ratio_cds_r','amino_acid_ratio_cds_h','amino_acid_ratio_cds_q','amino_acid_ratio_cds_m',
               'amino_acid_ratio_cds_l', 'amino_acid_ratio_cds_v','amino_acid_ratio_cds_g','amino_acid_ratio_cds_i')
  titles = c('Arginine','Histidine','Glutamine','Methionine',
             'Leucine','Valine','Glycine','Isoleucine')
  cols = 1:4
  rows = 1:2
  feature_data = data.frame()
  for(row in rows){
    for(col in cols){
      i = length(cols)*(row-1) + col
      f = features[i]
      title = titles[i]
      f1 = log(data_set_ultra_mega$pro_per_rna)
      f2 = data_set_ultra_mega[,f]
      f2 = f2 - min(c(data_set_ultra_mega[data_set_ultra_mega[,f]<0,f],0))
      f2 = f2 + min(f2[f2>0]/2)
      f2 = scale(log(f2))
      f2[!findInterval(f2,c(-10,10))==1] = 0
      c = cor.test(f1, f2, method='spearman')
      title = sprintf('%s\nrho: %s, p: %s', title, round(c$estimate,2), ifelse(c$p.value < 0.01, '<0.01', as.character(round(c$p.value,2))))
      d = densCols(f1, f2, colramp = colorRampPalette(viridis(10))) # colorRampPalette(rev(rainbow(10, end = 4/6))))
      temp = data.frame(x=f1, y=f2, col=title, d = d)
      feature_data = rbind(feature_data, temp)
    }
  }
  feature_data.cor <- ddply(feature_data, .(col), function(val) sprintf("rho=%.2f | p=%s", cor(val$x, val$y, method='spearman'), ifelse(cor.test(val$x, val$y, method='spearman')$p.value<0.01,'<0.01',round(cor.test(val$x, val$y, method='spearman')$p.value,2))))
  
  p <- ggplot(feature_data) + # Export 6x6 inch / 600x600 px
    geom_point(aes(x, y, col = d), size = 1) +
    scale_color_identity() +
    scale_alpha(range = c(0.1, 1)) +
    facet_wrap(~ col, ncol=4) + 
    theme_bw() + 
    theme(legend.position="top"
          , title = element_text(size = rel(1.5))
          , legend.text = element_text(size = rel(1.4))
          , axis.text = element_text(size=rel(1.4))
          , strip.text = element_text(size=rel(0.9))
    ) + 
    labs(x='RP ratio (log scaled)', y='Feature value (centered and scaled)') 
  
  print(p) 
  pdf(file = "figures/feature_vs_rpratio_top_aa.pdf", 6, 4)
  print(p) 
  dev.off()
  
  
  
  # Same but for all features in best 3 (plus unused categories)
  data_set_best3_3 = read.table('./feature_sets//data_set_x0_110110100-1111.csv', sep='\t', header = T)
  temp_adfwe = data.frame(feature = character(), rho = numeric(), pvalue = numeric())
  feature_data = data.frame()
  titles = PrettifyFeatureNames(colnames(data_set_best3_3))
  i = 3
  for(f in colnames(data_set_best3_3)[3:(length(colnames(data_set_best3_3)))]){
    title = titles[i]
    f1 = log(data_set_best3_3$pro_per_rna)
    f2 = data_set_best3_3[,f]
    f2 = f2 - min(c(data_set_best3_3[data_set_best3_3[,f]<0,f],0))
    f2 = f2 + min(f2[f2>0]/2)
    if(min(f2)!=max(f2)){
      f2 = scale(log(f2))
    }
    f2[!findInterval(f2,c(-10,10))==1] = 0
    c = cor.test(f1, f2, method='spearman')
    temp_adfwe = rbind(temp_adfwe, data.frame(feature=title, rho = c$estimate, pvalue=c$p.value))
    title = sprintf('%s\nrho: %s, p: %s', title, round(c$estimate,2), ifelse(c$p.value < 0.01, '<0.01', as.character(round(c$p.value,2))))
    d = densCols(f1, f2, colramp = colorRampPalette(viridis(10)))
    temp = data.frame(x=f1, y=f2, col=title, d=d)
    feature_data = rbind(feature_data, temp)
    i = i + 1
  }
  p <- ggplot(feature_data) + #exported at 1400 x 1800 px
    stat_density2d(aes(x,y), size=0.1, show.legend=F, color='#4444FF') + 
    stat_density2d(geom = "polygon", aes(x, y, alpha=..level..), contour=T, show.legend=F, fill='#5555FF') +
    facet_wrap(~ col, ncol=8) + 
    theme_bw() + 
    theme(axis.title.x = element_text(size = rel(1.5))
          , axis.title.y = element_text(size = rel(1.5))
          , axis.text = element_text(size = rel(1))
          , panel.spacing = unit(0, "lines")) + 
    labs(x='RP ratio (log scaled)', y='Feature value (centered and scaled)')
  
  print(p) 
  pdf(file = "figures/feature_vs_rpratio_allx.pdf", 14, 18)
  print(p) 
  dev.off()
}


GetFeatureImportancePlot = function(ds, xlab=T, ylab=T, lgnd=T, strlen=NA){

  
  scr = readRDS('data/ml_scores_filtered2.RObject')
  scr = scr[!grepl('\\_p[0-9]{1,2}', scr$data_set),] # remove permutation test results
  ds_data1 = data.frame(feature = character(), importance=numeric(), stringsAsFactors=FALSE) 
  ds_data_raw = read.table(sprintf('feature_sets/data_set_%s.csv',ds), header=T,sep='\t',stringsAsFactors=F)
  colnames(ds_data_raw) = PrettifyFeatureNames(colnames(ds_data_raw))
  ds_data1[1:(length(colnames(ds_data_raw))-2),'feature'] = colnames(ds_data_raw)[3:length(colnames(ds_data_raw))]
  ds_data2 = ds_data3 = ds_data1
  ds_data1[1:(length(colnames(ds_data_raw))-2),'importance'] = as.numeric(colSums(matrix(unlist(scr[scr$base_learner=='rfr' & scr$data_set==ds,'fi']), ncol=length(scr[scr$base_learner=='rfr' & scr$data_set==ds,'fi'][[1]]), byrow=T)) / sum(unlist(scr[scr$base_learner=='rfr' & scr$data_set==ds,'fi'])))
  ds_data2[1:(length(colnames(ds_data_raw))-2),'importance'] = abs(as.numeric(colSums(matrix(unlist(scr[scr$base_learner=='svm' & scr$data_set==ds,'fi']), ncol=length(scr[scr$base_learner=='svm' & scr$data_set==ds,'fi'][[1]]), byrow=T)) / sum(unlist(scr[scr$base_learner=='svm' & scr$data_set==ds,'fi'])))) / sum(abs(as.numeric(colSums(matrix(unlist(scr[scr$base_learner=='svm' & scr$data_set==ds,'fi']), ncol=length(scr[scr$base_learner=='svm' & scr$data_set==ds,'fi'][[1]]), byrow=T)) / sum(unlist(scr[scr$base_learner=='svm' & scr$data_set==ds,'fi'])))))
  ds_data3[1:(length(colnames(ds_data_raw))-2),'importance'] = abs(as.numeric(colSums(matrix(unlist(scr[scr$base_learner=='rir' & scr$data_set==ds,'fi']), ncol=length(scr[scr$base_learner=='rir' & scr$data_set==ds,'fi'][[1]]), byrow=T)) / sum(unlist(scr[scr$base_learner=='rir' & scr$data_set==ds,'fi'])))) / sum(abs(as.numeric(colSums(matrix(unlist(scr[scr$base_learner=='rir' & scr$data_set==ds,'fi']), ncol=length(scr[scr$base_learner=='rir' & scr$data_set==ds,'fi'][[1]]), byrow=T)) / sum(unlist(scr[scr$base_learner=='rir' & scr$data_set==ds,'fi'])))))

  return(BarPlotComparison(list(ds_data1, ds_data2, ds_data3)
                           , 'importance'
                           , 'feature'
                           , 'Feature'
                           , 'Importance'
                           , list('RFR','SVM','RIR')
                           , 'title'
                           , F, c(NA,NA), T
                           , ttl=F, xlb=xlab, ylb=ylab, lgnd=lgnd
                           , strlen = strlen))
}
CreateFeatureImportance = function(){
  pdf(file = "figures/feature_importance_full.pdf", 14, 18)
  fi_a1 = GetFeatureImportancePlot('a1', F, F, F, 24)
  fi_b3 = GetFeatureImportancePlot('b3', F, F, T, 11)
  fi_c7 = GetFeatureImportancePlot('c7', T, T, F, NA)
  fi_db7 = GetFeatureImportancePlot('db3', T, F, F, NA)
  fi_f1 = GetFeatureImportancePlot('f1', F, F, F, 24)
  fi_gE = GetFeatureImportancePlot('gE', T, F, F, 11)
  fi_gF = GetFeatureImportancePlot('gF', F, F, F, 11)
  fi_i1 = GetFeatureImportancePlot('i1', F, F, F, 24)
  
  # Fig XXX
  lay <- rbind(c(3,1,2),
               c(3,1,2),
               c(3,5,2),
               c(3,5,7),
               c(3,8,7),
               c(3,8,7),
               c(3,4,6),
               c(3,4,6),
               c(3,4,6),
               c(3,4,6),
               c(3,4,6),
               c(3,4,6),
               c(3,4,6),
               c(3,4,6),
               c(3,4,6),
               c(3,4,6),
               c(3,4,6),
               c(3,4,6))
  print(grid.arrange(fi_a1, fi_b3, fi_c7, fi_db7, fi_f1, fi_gE, fi_gF, fi_i1, layout_matrix=lay)) # export 14x18 inch
  dev.off()
  
  # Fig XXX
  pdf(file = "figures/feature_importance_3.pdf", 9, 4)
  fi_gE = GetFeatureImportancePlot('gE', T, T, F, NA)
  fi_f1 = GetFeatureImportancePlot('f1', F, F, T, NA)
  fi_i1 = GetFeatureImportancePlot('i1', T, F, F, NA)

  fi_gE$data$x = mapvalues(fi_gE$data$x, from=levels(fi_gE$data$x), to=paste('  ', toupper(substr(levels(fi_gE$data$x), 9, 9))))
  fi_f1$data$x = revalue(fi_f1$data$x, c('HLF protein halflife'='Prot. half-life', 'HLF rna halflife'='RNA half-life', 'HLF rna synthesistime'='RNA synth. rate'))
  fi_i1$data$x = revalue(fi_i1$data$x, c('RPF rpf/rna'='RPF / RNA', 'RPF rpf regression'=' RPF Regression'))
  
  lay = rbind(c(1,2),
              c(1,3))
  print(grid.arrange(fi_gE, fi_f1, fi_i1, layout_matrix=lay)) # export 9x4 inch
  dev.off()
  print(grid.arrange(fi_gE, fi_f1, fi_i1, layout_matrix=lay)) # export 9x4 inch
  
  # Fig XXX
  fi_comp_exp = GetFeatureImportancePlot('x0_110110100-1111', T, T, T, NA)

  pdf(file = "figures/feature_importance_comp-exp.pdf", 6, 16)
  print(fi_comp_exp) # export 9x4 inch
  dev.off()
  print(fi_comp_exp) # export 9x4 inch
}

CreateFeatureImportanceAFS = function(){ # For the Automatic Feature Selection methods
  pdf(file = "figures/feature_importance_auto_feat_sel.pdf", 14, 18)
  dfs = read.table('ml_output/fs_features.csv',header=F,sep=',',stringsAsFactors=F)
  data_set_ultra_mega = read.table('./feature_sets/data_set_ultra_mega2.csv', sep='\t', header = T)
  dfs$V3 = substring(dfs$V3, 2, nchar(dfs$V3))
  fn = PrettifyFeatureNames(colnames(data_set_ultra_mega))
  
  # SelectFromModel Elastic Net L1 = 0.95
  sfm_en3 = dfs[dfs$V1 == 'ultra_mega2_sfm_en0.95_0.05',]
  sfm_en3_features = c()
  for(i in 1:nrow(sfm_en3)){
    print(i)
    sfm_en3_features = c(sfm_en3_features, which(strsplit(sfm_en3[i,'V3'], '\\s+', perl=T)[[1]]=='True'))
  }
  temp = table(sfm_en3_features)
  
  asdf = data.frame(feature=character(), amount=numeric(), stringsAsFactors=F)
  temp2 = temp[temp>1]
  for(i in 1:length(temp2)){
    asdf[i,'feature'] = fn[as.numeric(names(temp2)[i])+2]
    asdf[i,'amount'] = temp2[i]
  }
  auto_fi1 = BarPlotComparison(list(asdf)
                               ,'amount'
                               ,'feature'
                               ,'Feature'
                               ,'Iterations'
                               ,list('SFM EN L1:0.95')
                               ,'FI SFM - ElasticNet3'
                               ,F
                               , c(NA,NA)
                               , F
                               , F)
  
  
  # SelectFromModel ExtraTreesRegressor
  sfm_et = dfs[dfs$V1 == 'ultra_mega2_sfm_et',]
  sfm_et_features = c()
  for(i in 1:nrow(sfm_et)){
    print(i)
    sfm_et_features = c(sfm_et_features, which(strsplit(sfm_et[i,'V3'], '\\s+', perl=T)[[1]]=='True'))
  }
  temp = table(sfm_et_features)
  
  asdf = data.frame(feature=character(), amount=numeric(), stringsAsFactors=F)
  temp2 = temp[temp>5]
  for(i in 1:length(temp2)){
    asdf[i,'feature'] = fn[as.numeric(names(temp2)[i])+2]
    asdf[i,'amount'] = temp2[i]
  }
  auto_fi2 = BarPlotComparison(list(asdf),'amount','feature','Feature','Iterations',list('SFM ETR'),'FI SFM ExtraTreesR.',F, c(NA,NA), F, F, T, F)
  
  grid.arrange(auto_fi1, auto_fi2, ncol=2)
  dev.off()
  grid.arrange(auto_fi1, auto_fi2, ncol=2)
}

# Fig XXX
CreatePermutationTest = function(){
  scr = readRDS('ml_scores_filtered2.RObject')
  
  perm_results = data.frame(scores=numeric(), type=character(), set=character(), stringsAsFactors=F)
  # Same order as heatmap
  order = c(1, #seq
            2,4,6,3,5, # nuc
            9,8,7, #cod
            28,25,24, # aa
            29,27,26, # aap
            23, # poly proline
            31,32,33, # rpf
            10,12,16,14,11,13,15, # struct
            19,22,21,20, # halflife
            30, # rna abdnc
            17,18 # sorf
  )
  data_sets = factor(unique(scr[!grepl('\\_p[0-9]{1,2}', scr$data_set),'data_set'])[order], ordered=F)
  
  perm_results=data.frame(scores=numeric(), type=character(), set=factor(levels=levels(data_sets)))
  
  for(ml in unique(scr$base_learner)){
    #ml='rfr'
    for(ds in data_sets){
      #ds='a1'
      if(nrow(scr[scr$data_set==sprintf('%s_p1',ds), ]) > 0){
        scores_original = scr[scr$data_set==ds & scr$base_learner==ml,'test_score']
        scores_permutation = c()
        for(i in 1:10){
          scores_permutation = c(scores_permutation, scr[scr$data_set==sprintf('%s_p%s',ds,i+1) & scr$base_learner==ml, 'test_score'])
        }
        pval = t.test(scores_permutation, scores_original)$p.value
        print(sprintf('%s -- %s ---- %s', ml, ds, pval))
        #hist(c(scores_original, scores_permutation), breaks=40)
        temp = data.frame(scores=scores_permutation, type='perm', set=sprintf('%s %s\np: %s', toupper(ml), PrettifyFeatureSetNames(ds),ifelse(pval < 0.01, '<0.01', as.character(round(pval,2)))))
        temp = rbind(temp, data.frame(scores=scores_original, type='orig', set=sprintf('%s %s\np: %s', toupper(ml), PrettifyFeatureSetNames(ds),ifelse(pval < 0.01, '<0.01', as.character(round(pval,2))))))
        perm_results = rbind(perm_results, temp)
      }
    } 
  }
  perm_results[perm_results$scores < -0.05, 'scores'] = -0.05 # We need some kind of minimum because some values are very large negative
  perm_results$set = gsub(perm_results$set, pattern=' - ', replacement='', fixed=T)
  perm_results = transform(perm_results, set=factor(set,levels=unique(perm_results$set))) # Force same order as we put them in the data set
  
  p = ggplot(perm_results, aes(scores))  + 
    geom_histogram(data=subset(perm_results,type == 'perm'),fill = "black", binwidth=0.01) +
    #geom_histogram(data=subset(perm_results,type == 'orig'),fill = "red", binwidth=0.001) + 
    geom_vline(data=subset(perm_results,type == 'orig'), aes(xintercept=scores), col='red', size=0.3) +
    facet_wrap(~ set, ncol=9) + 
    theme_bw() + 
    lims(x=c(-0.05, 0.21)) + 
    theme(axis.title.x = element_text(size = rel(1.5))
          , axis.title.y = element_text(size = rel(1.5))
          , panel.spacing = unit(0, "lines")
          , axis.text.x = element_text(angle=90))
  
  pdf(file = "figures/permutation_test.pdf", 14, 18)
  print(p)
  dev.off()
  print(p)
}

# Fig XXX
CreateRpkmComparison = function(){
  # GET EXON INFO
  d_annot = read.csv('data/mrnas_p/s_pombe_annot.ASM294v2.29.modified.bed', sep='\t',header=F)
  colnames(d_annot) = c('chromosome', 'start', 'end', 'V4', 'V5', 'strand','V7','region_type','V9','details')
  d_annot = d_annot[d_annot$region_typ=='exon',]
  d_annot = d_annot[d_annot$chromosome %in% c('III','II','I'),]
  details_matrix = split_to_matrix(d_annot$details,';')
  d_annot$transcript_id = substring(details_matrix[,2], 19)
  d_annot$length = d_annot$end - d_annot$start
  detach("package:plyr", unload=TRUE) # plyr and dplyr conflict
  require(dplyr)
  total_exon_lengths_and_counts_per_transcript = as.data.frame(d_annot %>%
                                                                 group_by(transcript_id) %>%
                                                                 summarize(total_exon_length=sum(length),
                                                                           exon_count = n()))
  
  # CREATE BASE DATA SET AND APPEND DATA FROM THE READ COUNTS DATA SETS
  d_ = total_exon_lengths_and_counts_per_transcript
  data_sets = c('marguerat_p','sub_rna_bt2','eser')
  for(ds in data_sets){
    d = read.csv(sprintf('data/htseq_counts_%s.txt',ds), header=F, sep='\t', stringsAsFactors=F)
    d$transcript_id = apply(d, 1, function(x) {return(  strsplit(x['V1'], ':', fixed=T)[[1]][[1]]  )})
    
    # BECAUSE THESE AREN'T IMPORTANT AT THE MOMENT AND SCREW UP OUR PLOTS
    d = d[regexpr('^SPRRNA.*',d$transcript_id)==-1,]
    d = d[regexpr('^SP.TRNA.*',d$transcript_id)==-1,]
    d = d[regexpr('^SPMIT.*',d$transcript_id)==-1,]
    d = d[regexpr('^SPSNORNA.*',d$transcript_id)==-1,]
    d = d[regexpr('^SPSNRNA.*',d$transcript_id)==-1,]
    
    d = d %>%
      group_by(transcript_id) %>%
      summarize(reads = sum(V2))
    d = merge(d,total_exon_lengths_and_counts_per_transcript,by="transcript_id")
    
    d$rpkm = calc_rpkm(d)
    d$rpkm_norm = d$rpkm / sum(d$rpkm)
    
    d_[,sprintf('%s_rpkm',ds)] = apply(d_, 1, function(x) {
      temp = d[d$transcript_id==x['transcript_id'],]
      return(ifelse(nrow(temp)==0,0,as.numeric(temp['rpkm'])))  })
    
    d_[,sprintf('%s_rpkm_norm',ds)] = apply(d_, 1, function(x) {
      temp = d[d$transcript_id==x['transcript_id'],]
      return(ifelse(nrow(temp)==0,0,as.numeric(temp['rpkm_norm'])))  })
  }
  
  c <- cor.test(d_$marguerat_p_rpkm_norm, d_$sub_rna_bt2_rpkm_norm, method='spearman')
  temp1 = d_$marguerat_p_rpkm
  temp1[temp1<=1] = 1
  temp2 = d_$sub_rna_bt2_rpkm
  temp2[temp2<=1] = 1
  m = summary(lm(log(temp1) ~ log(temp2)))
  temp = data.frame(x=d_$marguerat_p_rpkm, y=d_$sub_rna_bt2_rpkm, col=sprintf('Marguerat vs Subtelny\nr^2: %s, p: %s', round(m$r.squared,2), ifelse(m$coefficients[2,4] < 0.01, '< 0.01', as.character(round(m$m$coefficients[2,4],2)))))
  
  c <- cor.test(d_$marguerat_p_rpkm, d_$eser_rpkm, method='spearman')
  temp1 = d_$marguerat_p_rpkm
  temp1[temp1<=1] = 1
  temp2 = d_$eser_rpkm
  temp2[temp2<=1] = 1
  m = summary(lm(log(temp1) ~ log(temp2)))
  temp = rbind(temp, data.frame(x=d_$marguerat_p_rpkm, y=d_$eser_rpkm, col=sprintf('Marguerat vs Eser\nr^2: %s, p: %s', round(m$r.squared,2), ifelse(m$coefficients[2,4] < 0.01, '< 0.01', as.character(round(m$coefficients[2,4],2))))))
  
  c <- cor.test(d_$sub_rna_bt2_rpkm, d_$eser_rpkm, method='spearman')
  temp1 = d_$sub_rna_bt2_rpkm
  temp1[temp1<=0] = 1
  temp2 = d_$eser_rpkm
  temp2[temp2<=0] = 1
  m = summary(lm(log(temp1) ~ log(temp2)))
  temp = rbind(temp, data.frame(x=d_$sub_rna_bt2_rpkm, y=d_$eser_rpkm, col=sprintf('Subtelny vs Eser\nr^2: %s, p: %s', round(m$r.squared,2), ifelse(m$coefficients[2,4] < 0.01, '< 0.01', as.character(round(m$coefficients[2,4],2))))))
  temp$x = log(temp$x)
  temp$y = log(temp$y)
  temp$d = densCols(temp$x, temp$y, colramp = colorRampPalette(viridis(10)))
  temp[is.na(temp$d),'d'] = '#000000'
  
  p <- ggplot(temp) +
    geom_point(aes(x, y, col=d), size = 1) +
    scale_color_identity() +
    facet_wrap(~ col, ncol=8) + 
    theme_bw() + 
    theme(axis.title.x = element_text(size = rel(1.2)), axis.title.y = element_text(size = rel(1.2))
          , axis.text = element_text(size = rel(1.1))
          , strip.text = element_text(size=rel(1.1))
          , panel.spacing = unit(0, "lines")) + 
    labs(x='RPKM (left) (log scaled)', y='RPKM (right) (log scaled)')
  
  
  pdf(file = "figures/rpkm_comparison.pdf", 6, 3)
  print(p)
  dev.off()
  print(p)
  
}

CreateResultsTable = function(){
  
  predictions = GetProteinPredictions()
  
  # [1] COMP exp. | r² 0.62   COMP rpf | r² 0.62    COMP rna.hl | r² 0.58 COMP pro.hl | r² 0.59
  # [5] COMP | r² 0.57        COMP rna.sr | r² 0.57 Linear | r² 0.48     
  # 7 Levels: Linear | r² 0.48 COMP | r² 0.57 COMP rna.sr | r² 0.57 ... COMP exp. | r² 0.62
  
  ### DISTRIBUTIONS COMPARE
  t.test(predictions[predictions$col=='COMP | r² 0.57','log_avg_count_pro'], predictions[predictions$col=='COMP | r² 0.57','log_pred_pro'])$p.value
  t.test(predictions[predictions$col=='COMP | r² 0.57','log_avg_count_pro'], predictions[predictions$col=='COMP rna.sr | r² 0.57','log_pred_pro'])$p.value
  t.test(predictions[predictions$col=='COMP | r² 0.57','log_avg_count_pro'], predictions[predictions$col=='COMP rna.hl | r² 0.58','log_pred_pro'])$p.value
  t.test(predictions[predictions$col=='COMP | r² 0.57','log_avg_count_pro'], predictions[predictions$col=='COMP pro.hl | r² 0.59','log_pred_pro'])$p.value
  t.test(predictions[predictions$col=='COMP | r² 0.57','log_avg_count_pro'], predictions[predictions$col=='COMP rpf | r² 0.62','log_pred_pro'])$p.value
  t.test(predictions[predictions$col=='COMP | r² 0.57','log_avg_count_pro'], predictions[predictions$col=='COMP exp. | r² 0.62','log_pred_pro'])$p.value
  plot(predictions[predictions$col=='COMP | r² 0.57','log_avg_count_pro'], predictions[predictions$col=='COMP exp. | r² 0.62','log_pred_pro'])
  
  
  
  ### RP RATIO COMPARE
  scr = readRDS('ml_scores_filtered.RObject')
  scr = scr[!grepl('\\_p[0-9]{1,2}', scr$data_set),] # Remove all results from the permutation test
  scr = scr[grepl('x0_', scr$data_set),]
  scr = scr[scr$base_learner=='rfr',] 
  
  t.test(scr[scr$data_set=='x0_110110100','test_score'],scr[scr$data_set=='x0_110110100-0100','test_score'])$p.value
  t.test(scr[scr$data_set=='x0_110110100','test_score'],scr[scr$data_set=='x0_110110100-0010','test_score'])$p.value
  t.test(scr[scr$data_set=='x0_110110100','test_score'],scr[scr$data_set=='x0_110110100-0001','test_score'])$p.value
  t.test(scr[scr$data_set=='x0_110110100','test_score'],scr[scr$data_set=='x0_110110100-1000','test_score'])$p.value
  t.test(scr[scr$data_set=='x0_110110100','test_score'],scr[scr$data_set=='x0_110110100-1111','test_score'])$p.value
  for(fs in unique(scr$data_set)){
    p = t.test(scr[scr$data_set=='x0_110110100','test_score'],scr[scr$data_set==fs,'test_score'])$p.value
    fs = PrettifyFeatureSetNames(fs)
    print(paste(fs,"--- PVALUE: ", p))
  } 
  
  
  
  
  #### PROTEIN PREDICTIONS COMPARE (R2)
  learner= 'rfr'
  
  predictions = data.frame()
  training_ids = read.table('./data/training_ids2.csv', sep='\t', header=T, stringsAsFactors=F)
  mrnas_p = readRDS('./data/mrnas_p.RObject')
  feature_sets = c('linear_model', 'x0_110110100', 'x0_110110100-0100','x0_110110100-0010','x0_110110100-0001','x0_110110100-1000','x0_110110100-1111')
  r2s = data.frame()
  
  for(fs in feature_sets) {
    #fs='linear_model'
    for(i in 1:10){
      #i = 1
      #fs='linear_model'
      train_ids = training_ids[training_ids[,sprintf('set_%s',i)]==T,'sys_id']
      test_ids = training_ids[training_ids[,sprintf('set_%s',i)]==F,'sys_id']
      
      if(fs=='linear_model'){
        print('a')
        mrnas_train = mrnas_p[mrnas_p$sys_id %in% train_ids,]
        mrnas_test = mrnas_p[mrnas_p$sys_id %in% test_ids,]
        
        m2 = lm(mrnas_train$avg_count_pro ~ mrnas_train$avg_count_rna)
        
        real_values = mrnas_test$pro_per_rna
        pred_values = rep(m2$coefficients[2], length(real_values))
        pred_temp = data.frame(mrnas_test$sys_id, real_values, pred_values)
        
        colnames(pred_temp) = c('sys_id','real','pred')
        
        pred_temp$avg_count_pro = apply(pred_temp, 1, function(x) {  return(mrnas_p[mrnas_p$sys_id==x['sys_id'], 'avg_count_pro'])   }  )
        pred_temp$avg_count_rna = apply(pred_temp, 1, function(x) {  return(mrnas_p[mrnas_p$sys_id==x['sys_id'], 'avg_count_rna'])   }  )
        pred_temp$pred_pro = pred_temp$avg_count_rna * pred_temp$pred
        
        pred_temp$clr = '#FF5577'
        
      } else {
        pred_temp = read.table(sprintf('./ml_output/ml_predictions_%s_%s_%s.csv', learner,fs,i),header=F,sep=',', stringsAsFactors=F)
        pred_temp = cbind(as.data.frame(test_ids[1:951]), pred_temp, stringsAsFactors=F) # ORDER STAYED THE SAME
        colnames(pred_temp) = c('sys_id','real','pred')
        
        # DE-SCALE THE DATA
        descaled_real = (pred_temp$real * sd(mrnas_p[mrnas_p$sys_id %in% train_ids, 'pro_per_rna'])) + mean(mrnas_p[mrnas_p$sys_id %in% train_ids, 'pro_per_rna'])
        m = lm(mrnas_p[mrnas_p$sys_id %in% test_ids[1:951], 'pro_per_rna'] ~ descaled_real)
        descaled_real = descaled_real * m$coefficients[2] + m$coefficients[1] # should match the real rp_ratio in the mrnas_p data set
        pred_temp$real = descaled_real
        
        descaled_pred = (pred_temp$pred * sd(mrnas_p[mrnas_p$sys_id %in% train_ids, 'pro_per_rna'])) + mean(mrnas_p[mrnas_p$sys_id %in% train_ids, 'pro_per_rna'])
        descaled_pred = descaled_pred * m$coefficients[2] + m$coefficients[1] # De-scale using only the training data!
        pred_temp$pred = as.numeric(as.character(descaled_pred))
        
        pred_temp$clr = '#006666'
      }
      
      pred_temp$avg_count_pro = apply(pred_temp, 1, function(x) {  return(mrnas_p[mrnas_p$sys_id==x['sys_id'], 'avg_count_pro'])   }  )
      pred_temp$avg_count_rna = apply(pred_temp, 1, function(x) {  return(mrnas_p[mrnas_p$sys_id==x['sys_id'], 'avg_count_rna'])   }  )
      pred_temp$pred_pro = pred_temp$avg_count_rna * pred_temp$pred
      
      pred_temp$log_avg_count_pro = log(pred_temp$avg_count_pro)
      pred_temp$log_pred_pro = log(pred_temp$pred_pro)
      
      r2s = rbind(r2s, data.frame(data_set=fs, test_set=i, r2=CalcR2(pred_temp$log_avg_count_pro, pred_temp$log_pred_pro)))
    }
  }
  
  t.test(r2s[r2s$data_set=='x0_110110100','r2'],r2s[r2s$data_set=='x0_110110100-0100','r2'])$p.value
  t.test(r2s[r2s$data_set=='x0_110110100','r2'],r2s[r2s$data_set=='x0_110110100-0010','r2'])$p.value
  t.test(r2s[r2s$data_set=='x0_110110100','r2'],r2s[r2s$data_set=='x0_110110100-0001','r2'])$p.value
  t.test(r2s[r2s$data_set=='x0_110110100','r2'],r2s[r2s$data_set=='x0_110110100-1000','r2'])$p.value
  t.test(r2s[r2s$data_set=='x0_110110100','r2'],r2s[r2s$data_set=='x0_110110100-1111','r2'])$p.value
  
  
  for(fs in unique(r2s$data_set)){
    p = t.test(r2s[r2s$data_set=='x0_110110100','r2'],r2s[r2s$data_set==fs,'r2'])$p.value
    fs = PrettifyFeatureSetNames(fs)
    print(paste(fs,"--- PVALUE: ", p))
  } 
}

CreateHeatmap()
CreateRpRatioPrediction()
CreateProteinAbundancePredictions()
CreateFeatureVsRPRatio()
CreateFeatureImportance()
CreateFeatureImportanceAFS()
CreatePermutationTest()
CreateRpkmComparison()






