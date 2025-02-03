read96wMaps = function(folder='./input/maps/') {
  lapply(list.files(folder, pattern='library_plate', full.names=T), fread) %>%
    rbindlist()
}


registerQuadrants = function(form='384w', plate_number = NULL, bio_rep = NULL,
  tech_rep = NULL) {

  if(!form %in% c('384w', '1536w')) stop("\nUnknown plate format. Should be '384w' or '1536'")
  if(is.null(plate_number)) stop('\nWas expecting plate numbers!')
  if(length(plate_number)%%4!=0) stop('\nWas expecting plates for all four quadrants!')

  if(is.null(bio_rep)) bio_rep=rep(1, length(plate_number))
  if(is.null(tech_rep)) tech_rep=rep(1, length(plate_number))

  if(length(plate_number)!=length(bio_rep)) stop('\nNumber of bioreps does not match the number of plates!')
  if(length(plate_number)!=length(tech_rep)) stop('\nNumber of techreps does not match the number of plates!')

  n = length(plate_number)/4
  mat = matrix(plate_number, nrow=4)
  mat = cbind(1:4, mat)

  plates = paste0('plt_', seq_len(n))
  colnames(mat) = c('qrt', plates)

  out = as.data.frame(mat) %>%
    pivot_longer(!contains('qrt'), names_to='to', values_to='from') %>%
    arrange(to, qrt) %>%
    # tidy up 
    mutate(to = gsub('plt_', '', to)) %>%
    mutate(across(everything(), as.numeric))
  
  out = cbind(out, bio_rep, tech_rep)

  if(form=='384w') {
    names(out) = c('qrt384', 'plt384', 'plt96', 'biorep96', 'techrep96')
    if(length(unique(out$plt384))==1) out = subset(out, select=-plt384)
  } else {
    names(out) = c('qrt1536', 'plt1536', 'plt384', 'biorep384', 'techrep384')
    if(length(unique(out$plt1536))==1) out = subset(out, select=-plt1536)
  }

  out
}



add384RowAndCol = function(dat) {

  # add 384w coordinates 
  rowcol384 = rbind(
    data.table(
      qrt384 = 1,
      row96    = rep(LETTERS[1:8], each=12),
      col96    = rep(1:12, 8), 
      row384   = rep(LETTERS[seq(1, 16, 2)], each=12),
      col384   = rep(seq(1, 24, 2), 8)
    ),
    
    data.table(
      qrt384 = 2,
      row96    = rep(LETTERS[1:8], each=12),
      col96    = rep(1:12, 8), 
      row384   = rep(LETTERS[seq(1, 16, 2)], each=12),
      col384   = rep(seq(2, 24, 2), 8)
    ),
    
    data.table(
      qrt384 = 3,
      row96    = rep(LETTERS[1:8], each=12),
      col96    = rep(1:12, 8), 
      row384   = rep(LETTERS[seq(2, 16, 2)], each=12),
      col384   = rep(seq(1, 24, 2), 8)
    ),
  
    data.table(
      qrt384 = 4,
      row96    = rep(LETTERS[1:8], each=12),
      col96    = rep(1:12, 8), 
      row384   = rep(LETTERS[seq(2, 16, 2)], each=12),
      col384   = rep(seq(2, 24, 2), 8)
    )
  )
  left_join(dat, rowcol384, by=c('qrt384'), relationship='many-to-many')
}

add1536RowAndCol = function(dat) {
  # add 1536w coordinates 
  rowcol1536 = rbind(
    data.table(
      qrt1536  = 1, 
      row384   = rep(LETTERS[1:16], each=24),
      col384   = rep(seq(1, 24), 16),
      row1536  = rep(seq(1, 32, 2), each=24),
      col1536  = rep(seq(1, 48, 2), 16)
    ),
    
    data.table(
      qrt1536  = 2, 
      row384   = rep(LETTERS[1:16], each=24),
      col384   = rep(seq(1, 24), 16),
      row1536  = rep(seq(1, 32, 2), each=24),
      col1536  = rep(seq(2, 48, 2), 16)
    ),
    
    data.table(
      qrt1536 = 3, 
      row384   = rep(LETTERS[1:16], each=24),
      col384   = rep(seq(1, 24), 16),
      row1536  = rep(seq(2, 32, 2), each=24),
      col1536  = rep(seq(1, 48, 2), 16)
    ),
    
    data.table(
      qrt1536 = 4, 
      row384   = rep(LETTERS[1:16], each=24),
      col384   = rep(seq(1, 24), 16),
      row1536  = rep(seq(2, 32, 2), each=24),
      col1536  = rep(seq(2, 48, 2), 16)
    )
  )
  left_join(dat, rowcol1536, by=c('qrt1536', 'row384', 'col384'), relationship='many-to-many')
}

loadIrisFiles = function(folder) {
  name = list.files(folder, pattern='.iris$')
  path = list.files(folder, pattern='.iris$', full.names=T)
  out = lapply(path, fread)
  names(out) = name
  rbindlist(out, idcol='file')
}


plotReplicateCorrelation = function(dat) {

  # remove factor levels with all NA's ------
  percentage_of_nas = dat %>% 
    group_by(cond) %>%
    select(contains('rep')) %>%
    summarize(across(everything(), ~mean(is.na(.)))) 
  
  percentage_of_nas[percentage_of_nas==1] = NA
  cols_to_remove = percentage_of_nas %>% select(where(~any(is.na(.)))) %>%
    select(contains('rep')) %>%
    colnames()

  dat = select(dat, !cols_to_remove)
  # ------------------------------------

  # plot replicate correlation
  n_plots = ncol(dat) - 4

  has_zero = apply(dat[, 4:ncol(dat)] == 0, 1, sum, na.rm=T) != 0

  dat[!has_zero] %>%
    #filter(`rep1` > 0 & `rep2` > 0)  %>%

    mutate(across(is.numeric, ~log10(.x + thr))) %>%
  {
    ggpairs(., columns=4:ncol(.), aes(color=cond, alpha=0.3), 
      lower = list(size = 0.5/n_plots),
      upper = list(continuous = wrap('cor', size = 6/sqrt(n_plots))),
    ) +
    xlab(expression(log[10]*(opacity))) +
    ylab(expression(log[10]*(opacity))) +
    scale_color_brewer(palette='Dark2') +
    theme_bw()
  }
}


plotHC = function(dat) {
  # plot hierachical clustering

  # remove factor levels with all NA's ------
  percentage_of_nas = dat %>% 
    group_by(cond) %>%
    select(contains('rep')) %>%
    summarize(across(everything(), ~mean(is.na(.)))) 
  
  percentage_of_nas[percentage_of_nas==1] = NA
  cols_to_remove = percentage_of_nas %>% select(where(~any(is.na(.)))) %>%
    select(contains('rep')) %>%
    colnames()

  dat = select(dat, !cols_to_remove)
  # ------------------------------------

  # make data wider still 
  dat_hc = dat %>%
    pivot_longer(cols=4:ncol(.), names_to = 'rep') %>%
    pivot_wider(
      id_cols=genename,
      names_from = c(cond, numb, rep),
      values_from = value
    ) %>%
    setDT()

  m = dat_hc[, -1] %>% as.matrix()
  
  rownames(m) = dat_hc[['genename']]
  p = colnames(m)
  
  meta = data.frame(foo = p) %>%
    separate_wider_delim(
      foo, '_', names = c('cond', 'conc', 'rep')) %>% 
    data.frame(., row.names = p)
  
  #mat = m %>% cor(use='complete.obs') 
  #
  #o = rownames(mat)
  #hc = hclust(as.dist(1 - mat))  # as.dist: use the correlation distance
  #mat = mat[hc$order, hc$order]
  #mat[lower.tri(mat)] = NA
  #mat = mat[o, o]

  getBrewerColors = function(var_from_meta, palette_name) {
    varname = unique(meta[[var_from_meta]])
    colorRampPalette(RColorBrewer::brewer.pal(n=8, palette_name))(length(varname)) %>%
      setNames(varname)
  }
  
  annCols = list(
    cond = getBrewerColors('cond', 'Dark2'),
    conc = getBrewerColors('conc', 'Accent'),
    rep = getBrewerColors('rep', 'Paired')
  )

  #m[is.na(m)] = 0 - thr

  pheatmap(
    log10(m + thr),  # assume multiplicative errors
    annotation_col=select(meta, rep, conc, cond),
    annotation_colors = annCols,
    show_colnames=T,
    fontsize_row=2, 
    fontsize_col=2, 
    #clustering_method='single',
    na_col = '#FFFFFF',
    heatmap_legend_param = list(title = expression(log[10]*(opacity)))
  )
}


getResultsFromLinearModel = function(dat, folder, plot_pca_qc=T) {

  x = 'biorep_all'

  dat = dat %>%
    filter(numb %in% c('10010201', '10010500201', '1001001', '100105001')) %>%
    mutate(cond = gsub('Spectet|Spectetamp', '', cond))


  dat_long = dat %>% 
    # make sure to retain values_drop_na or the plotMDS will not work
    pivot_longer(contains('rep'), names_to='rep', values_to = 'opacity', values_drop_na = T) %>%
    mutate(rep = gsub('ep', '', rep)) %>%
    group_by(cond, numb, rep) %>%
    mutate(
      fitness = opacity/median(opacity[genename=='gfp'])
    )

  # Remove the ones you have information only from one condition
  # Why? lmFit can report back logFC NA, but method='robust' will trip 'rlm' up
  only_zeros_in_a_cond = dat_long %>%
    group_by(cond, genename) %>%
    summarize(mu = mean(fitness, na.rm=T)) %>%
    filter(mu == 0) %>%
    pull(genename)

  dat_long = dat_long %>% filter(!genename %in% only_zeros_in_a_cond)

  fitness = dat_long %>%
    filter(genename!='gfp') %>%
    select(-opacity) %>%
    pivot_wider(
      id_cols=genename,
      names_from = c(cond, numb, rep),
      values_from = fitness
    ) %>%
    setDT() 

  m = fitness[, -1] %>% as.matrix() %>% log2()
  genes = fitness[['genename']]
  rownames(m) = genes
  f = data.frame(genes)
  rownames(f) = genes
  p = colnames(m)
  
  meta = data.frame(foo = p) %>%
    separate_wider_delim(
      foo, '_', names = c('cond', 'conc', 'rep')) %>% 
    data.frame(., row.names = p) %>%
    mutate(rep = gsub('ep', '', rep))
  
  eset = ExpressionSet(m, AnnotatedDataFrame(meta), AnnotatedDataFrame(f))
 
  if(plot_pca_qc) {
    pdf(paste0('./output/', folder, '/qc_', x, '_pca.pdf'), h=5, w=10)
    par(mfrow=c(1, 2), pty='s')
    plotMDS(eset, labels=pData(eset)[, 'rep'], cex=0.75)
    plotMDS(eset, labels=pData(eset)[, 'cond'], cex=0.75)
    dev.off()
  }
  
  design = model.matrix(~ 0 + cond, data=pData(eset))
  colnames(design) = gsub('cond', '', colnames(design))
  cm = makeContrasts(ara = AraIPTG - IPTG, levels=design)

  fit = lmFit(eset, method='robust', design, maxit=10000)
  fit2 = contrasts.fit(fit, cm)
  fit2 = eBayes(fit2)
  res = topTable(fit2, genelist=fit2$genes, number = Inf, sort.by = 'none') %>%
    janitor::clean_names()

  res %>% {
    ggplot(., aes(log_fc, -log10(adj_p_val))) +
      geom_point(pch=21, col='grey30', fill='grey60') +
      ggrepel::geom_text_repel(
        data = filter(., adj_p_val < 0.01 & abs(log_fc) > 0.2 | abs(log_fc)>0.3 ),
        aes(label=genes)
      ) +
    labs(x=expression(log[2]*'FC'), y=expression(-log[10]*'(adj p-val)')) +
    theme_cowplot() 
  }

  ggsave(paste0('./output/', folder, '/res_volcano_', x, '.pdf'))

  fwrite(res, paste0('./output/', folder, '/res_', x, '.csv'))
}

