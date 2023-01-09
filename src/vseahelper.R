vsea <- function(af.thr,cadd.thr,epitype,mut, weighted = 0, mode = 'normal', gene.list = NULL){
  
  if(mode == 'normal'){
    
    # low.freq <- info.expanded %>% filter(AF_CTRL <= af.thr,
    #                                      CADD >= cadd.thr, GROUPS == epitype) %>% 
    #   dplyr::select(SYMBOL,HGVSP,ESTIMATES)
    
    low.freq <- info.expanded %>% filter(AF_CTRL <= af.thr, GROUPS == epitype) %>% 
      dplyr::select(SYMBOL,HGVSP,ESTIMATES)
    
    low.freq <- low.freq[!duplicated(low.freq[c("SYMBOL","HGVSP")]),]
    
    high.freq <- info.expanded %>% filter(AF_CASE > AF_CTRL,
                                          PVAL <= 0.05,
                                          MUTATION == mut, 
                                          GROUPS == epitype) %>%
      dplyr::select(SYMBOL,HGVSP)
    
    high.freq <- high.freq[!duplicated(high.freq[c("SYMBOL","HGVSP")]),]  
    
  } else if(mode == 'penalized') {
    
    low.freq <- info.expanded %>% filter(AF_CTRL <= af.thr,
                                         CADD >= cadd.thr, GROUPS == epitype) %>%
      mutate(ESTIMATES = log10(ESTIMATES)*-log10(PVAL)) %>%
      dplyr::select(SYMBOL,HGVSP,ESTIMATES)
    
    low.freq <- low.freq[!duplicated(low.freq[c("SYMBOL","HGVSP")]),]
    
    high.freq <- info.expanded %>% filter(AF_CASE > AF_CTRL,
                                          MUTATION == mut, GROUPS == epitype) %>%
      dplyr::select(SYMBOL,HGVSP)
    
    high.freq <- high.freq[!duplicated(high.freq[c("SYMBOL","HGVSP")]),]
  } else if (mode == 'filtered'){

      low.freq <- info.expanded %>% filter(AF_CTRL <= af.thr,
                                         CADD >= cadd.thr, GROUPS == epitype) %>% 
      dplyr::select(SYMBOL,HGVSP,ESTIMATES)
    
    low.freq <- low.freq[!duplicated(low.freq[c("SYMBOL","HGVSP")]),]
    
    high.freq <- info.expanded %>% filter(AF_CASE > AF_CTRL, ESTIMATES >= 2, PVAL <= 0.05,
                                          MUTATION == mut, GROUPS == epitype) %>%
      dplyr::select(SYMBOL,HGVSP)
    
    high.freq <- high.freq[!duplicated(high.freq[c("SYMBOL","HGVSP")]),]
  } else if (mode == 'modified'){
    
    low.freq <- info.expanded %>% filter(AF_CTRL <= af.thr,AN_CTRL > 0,
                                         CADD >= cadd.thr, GROUPS == epitype, SYMBOL %in% gene.list) %>% 
      dplyr::select(SYMBOL,HGVSP,ESTIMATES)
    
    low.freq <- low.freq[!duplicated(low.freq[c("SYMBOL","HGVSP")]),]
    
    high.freq <- info.expanded %>% filter(AF_CASE > AF_CTRL, ESTIMATES >= 2, PVAL <= 0.05,
                                          MUTATION == mut, GROUPS == epitype, SYMBOL %in% gene.list) %>%
      dplyr::select(SYMBOL,HGVSP)
    
    high.freq <- high.freq[!duplicated(high.freq[c("SYMBOL","HGVSP")]),]
    
  } else if (mode == 'cadd'){
    
    low.freq <- info.expanded %>% filter(AF_CTRL <= af.thr, GROUPS == epitype) %>%
      mutate(ESTIMATES = CADD) %>%
      dplyr::select(SYMBOL,HGVSP,ESTIMATES)
    
    low.freq <- low.freq[!duplicated(low.freq[c("SYMBOL","HGVSP")]),]
    
    high.freq <- info.expanded %>% filter(AF_CASE > AF_CTRL, ESTIMATES > 1, PVAL <= 0.05,
                                          MUTATION == mut, GROUPS == epitype) %>%
      dplyr::select(SYMBOL,HGVSP)
    
    high.freq <- high.freq[!duplicated(high.freq[c("SYMBOL","HGVSP")]),]
    
  } else if (mode == 'gsea'){
    
    low.freq <- info.expanded %>% 
      dplyr::select(SYMBOL,HGVSP,ESTIMATES)
    
    low.freq <- low.freq[!duplicated(low.freq[c("SYMBOL","HGVSP")]),]
    
    high.freq <- info.expanded %>% filter(AF_CTRL <= af.thr,AN_CTRL > 0,
                                          CADD >= cadd.thr,
                                          AF_CASE > AF_CTRL, ESTIMATES >= 2, 
                                          PVAL <= 0.05,
                                          AN_CTRL > 0,
                                          MUTATION == mut, 
                                          GROUPS == epitype, 
                                          SYMBOL %in% gene.list) %>%
      dplyr::select(SYMBOL,HGVSP)
    
    high.freq <- high.freq[!duplicated(high.freq[c("SYMBOL","HGVSP")]),]
    
  }
  
  
  
  
  low.freqnames <- as.data.table(low.freq)[, toString(ESTIMATES), by = list(SYMBOL)] %>% as.data.frame
  low.freq <- as.data.table(low.freq)[, toString(HGVSP), by = list(SYMBOL)] %>% as.data.frame
  high.freq <- as.data.table(high.freq)[, toString(HGVSP), by = list(SYMBOL)] %>% as.data.frame
  
  colnames(low.freq) <- c("GENES","LF_VARIANTS")
  colnames(low.freqnames) <- c("GENES","LF_ESTIMATES")
  colnames(high.freq) <- c("GENES","HF_VARIANTS")
  
  variant.summary <- NULL
  variant.summary <- inner_join(low.freq,low.freqnames, "GENES")
  variant.summary <- inner_join(variant.summary,high.freq, "GENES")
  
  
  variant.summary <- variant.summary %>% mutate(
    LF_VARIANTS = as.list(strsplit(LF_VARIANTS,",")),
    HF_VARIANTS = as.list(strsplit(HF_VARIANTS,",")),
    LF_ESTIMATES = as.list(strsplit(LF_ESTIMATES,",")),
    LF_ESTIMATES_COUNT = sapply(LF_ESTIMATES,length),
    LF_COUNT = sapply(LF_VARIANTS,length),
    HF_COUNT = sapply(HF_VARIANTS,length)
  ) %>% dplyr::select(GENES,LF_COUNT,HF_COUNT,everything())
  
  
  names(variant.summary$LF_VARIANTS) <- variant.summary$GENES
  names(variant.summary$HF_VARIANTS) <- variant.summary$GENES
  names(variant.summary$LF_ESTIMATES) <- variant.summary$GENES
  
  VSEA <- NULL
  for (variant in names(variant.summary$LF_VARIANTS)) {
    var.vector <- variant.summary$LF_ESTIMATES[[variant]] %>% as.numeric
    names(var.vector) <- variant.summary$LF_VARIANTS[[variant]]
    tmp <- cbind(TERM = variant,gsea.compute(var.vector,
                                             variant.summary$HF_VARIANTS[[variant]],
                                             w = weighted, nperm = 1000))

    VSEA <- rbind(VSEA,tmp)
   
  }
  
  
  qvalues <- NULL
  
  
  for (GENE in VSEA$TERM) {
    dataset <- VSEA %>% filter(TERM == GENE)
    qvalue <- gsea.fdr(dataset$NES,
                    str_split(dataset$PERM,",") %>% unlist %>% as.numeric,
                    VSEA$NES)
    qvalues <- c(qvalues,qvalue)
    
  }
  
  VSEA$FDR <- qvalues
  
  return(VSEA)
}
