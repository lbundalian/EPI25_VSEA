# Create geneset ----------------------------------------------------------
# this is the list of all genes involved per Epilepsy group 
signature.gene<- list(DEE = VSEA_List_MODIFIED$DEE_NONSYNONYMOUS$TERM,
                      NAFE = VSEA_List_MODIFIED$NAFE_NONSYNONYMOUS$TERM,
                      GGE = VSEA_List_MODIFIED$GGE_NONSYNONYMOUS$TERM)


# this is the list of genes relagted to specific paths
gobp <- new("GeneSet", "data/c5.go.bp.v7.5.1.symbols.gmt", )
gobp <- gobp@geneset




# Create gene set  --------------------------------------------------------
godf <- data.frame(PATH = NULL,GENES = NULL)
for (i in names(gobp)) {
  print(i)
  print(gobp[[i]] %>% paste0(.,collapse = ","))
  godf <- rbind(godf,data.frame(PATH = i,
                                GENES = gobp[[i]] %>% unlist %>%
                                  paste0(.,collapse = ",")))
}


godf <- separate_rows(godf, all_of(c("GENES")), sep=",")


brain.detected <- brain.genes.df %>% filter(`RNA brain regional distribution` %in% c('Detected in all','Detected in many') )


gene.set <- NULL
gene.set <- list(
  'ION CHANNEL AND INTERACTORS' = godf %>% filter(PATH %like% "_ION") %>%
    select(GENES) %>% unlist %>% unique,
  'NEUROTRANSMITTER REGULATORS' = godf %>% filter(PATH %like% "NEUROTRANSMI") %>%
    select(GENES) %>% unlist %>% unique,
  'NERVOUS SYSTEM DEVELOPMENT' = c(godf %>% filter(PATH %like% "NERVOUS_SYSTEM_DEV") %>%
    select(GENES) %>% unlist %>% unique,godf %>% filter(PATH %like% "BRAIN_DEV") %>%
      select(GENES) %>% unlist %>% unique) %>% unique,
  'SYNAPTIC FUNCTION MODULATORS' = godf %>% filter(PATH %like% "SYNAP") %>%
    select(GENES) %>% unlist %>% unique,
  'GABAERGIC SIGNALING MODULATORS' = godf %>% filter(PATH %like% "GABA") %>%
    select(GENES) %>% unlist %>% unique,
  'GLUTAMATERGIC SIGNALING MODULATORS' = godf %>% filter(PATH %like% "GLUTAMAT") %>%
    select(GENES) %>% unlist %>% unique,
  'BRAIN EXPRESSED' = brain.detected$Gene %>% unlist
  
)


fet.res <- NULL
for (s in names(signature.gene)) {
  for (g in names(gene.set)) {
    tmp <- enrich.fet(
      pathway = gene.set[[g]],
      overlap = length(signature.gene[[s]][signature.gene[[s]] %in% gene.set[[g]]]),
      signature = signature.gene[[s]],
      backpop = 26588
    )
    fet.res <- rbind(fet.res,data.frame(GROUP = s,PATH = g, 
                                        PVAL = tmp$p.value,
                                        HITS = paste(signature.gene[[s]][signature.gene[[s]] %in% gene.set[[g]]],collapse = ","),
                                        ESTIMATE = tmp$estimate,
                                        UPPER = tmp$conf.int[2],
                                        LOWER = tmp$conf.int[1]))
  }
}

fet.res$FDR <- p.adjust(fet.res$PVAL, method = "fdr")
View(fet.res)
fet.res <- fet.res %>% mutate(SIG = case_when(
  FDR <= 0.01 ~ "**",
  FDR > 0.01 & FDR <= 0.05 ~ "*",
  TRUE ~ ""
))



png("CUSTOM_ORA.png", width = 42, height = 21, units = 'cm', res = 300)
ggplot(data = fet.res %>% filter(PVAL < 1) , aes(y = PATH, x = log10(ESTIMATE),
                       xmin = log10(LOWER+0.001),
                       xmax = log10(UPPER+0.001))) +
  geom_point(aes(col = PATH),size = 10) +
  geom_errorbarh(height = .1, aes(col = PATH)) +
  geom_vline(xintercept = 0, color='black', linetype='dashed', alpha=.8) +
  geom_text_repel( aes( label = SIG), direction = "y",size = 8,
                   point.padding = unit(8, "points"), nudge_y = 0.055, color = "red") +
  facet_wrap(~GROUP) +
  xlim(-3, 3) +
  ggtitle("") +
  xlab("Log Enrichment Scores") +
  ylab("") +
  theme_bw() +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(2), hjust = 0.5),
        axis.title = element_text(size = rel(2)),
        axis.text.x = element_text(angle = 90,hjust = 0.95,vjust = 0.2, size =rel(2)),
        axis.text.y = element_text(size = rel(2),face = "bold")) +
  scale_fill_brewer(palette = "Set1")

dev.off()




# DEE
dee.gene.path <- expandInfo("HITS",fet.res %>% filter(PVAL != 1,GROUP == "DEE"))
dee.gene.path <- dee.gene.path %>%
  group_by(HITS) %>%
  do(.,data.frame(PATH = paste(.$PATH,collapse=","),stringsAsFactors=FALSE)) %>%
  as.data.frame() 
colnames(dee.gene.path) <- c("GENES","CATEGORIES")

dee.gene.sig <- as.data.frame(signature.gene$DEE)
colnames(dee.gene.sig) <- c("GENES")
dee.gene.path <- full_join(dee.gene.path,dee.gene.sig,by = "GENES")
dee.gene.path <- dee.gene.path %>% mutate(
  CATEGORIES = case_when(
    is.na(CATEGORIES) ~ "UNCATEGORIZED",
    TRUE ~ CATEGORIES
    )
)


gge.gene.path <- expandInfo("HITS",fet.res %>% filter(PVAL != 1,GROUP == "GGE"))
gge.gene.path <- gge.gene.path %>%
  group_by(HITS) %>%
  do(.,data.frame(PATH = paste(.$PATH,collapse=","),stringsAsFactors=FALSE)) %>%
  as.data.frame() 
colnames(gge.gene.path) <- c("GENES","CATEGORIES")

gge.gene.sig <- as.data.frame(signature.gene$GGE)
colnames(gge.gene.sig) <- c("GENES")
gge.gene.path <- full_join(gge.gene.path,gge.gene.sig,by = "GENES")
gge.gene.path <- gge.gene.path %>% mutate(
  CATEGORIES = case_when(
    is.na(CATEGORIES) ~ "UNCATEGORIZED",
    TRUE ~ CATEGORIES
  )
)

nafe.gene.path <- expandInfo("HITS",fet.res %>% filter(PVAL != 1,GROUP == "NAFE"))
nafe.gene.path <- nafe.gene.path %>%
  group_by(HITS) %>%
  do(.,data.frame(PATH = paste(.$PATH,collapse=","),stringsAsFactors=FALSE)) %>%
  as.data.frame() 
colnames(nafe.gene.path) <- c("GENES","CATEGORIES")

nafe.gene.sig <- as.data.frame(signature.gene$NAFE)
colnames(nafe.gene.sig) <- c("GENES")
nafe.gene.path <- full_join(nafe.gene.path,nafe.gene.sig,by = "GENES")
nafe.gene.path <- nafe.gene.path %>% mutate(
  CATEGORIES = case_when(
    is.na(CATEGORIES) ~ "UNCATEGORIZED",
    TRUE ~ CATEGORIES
  )
)

epi.enrichment.list <- list(
  NAFE = nafe.gene.path,
  DEE = dee.gene.path,
  GGE = gge.gene.path
)


View(fet.res)
sort.fet.res <- with(fet.res, fet.res[order(GROUP, PATH, HITS) , ])

View(sort.fet.res)

custom_table <- fet.res %>% select(GROUP,PATH,HITS)
custom_table_wide <- reshape(custom_table, idvar = "PATH", timevar = "GROUP", direction = "wide")
colnames(custom_table_wide) <- c("PATH","DEE","NAFE","GGE")



for_sort <- expandInfo("HITS",fet.res)
View(for_sort)

custom_table_sorted <- sort.fet.res %>% select(GROUP,PATH,HITS)
custom_table_wide_sorted <- reshape(custom_table_sorted, idvar = "PATH", timevar = "GROUP", direction = "wide")
colnames(custom_table_wide_sorted) <- c("PATH","DEE","NAFE","GGE")
View(custom_table_wide_sorted)




openxlsx::write.xlsx(custom_table_wide,"Table_1_FUNCTIONAL_CLASS.xlsx")

