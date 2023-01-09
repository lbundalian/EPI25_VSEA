
get.estimates <- function(gene.list = NULL){
  
  dee.variants <- separate_rows(dee.filtered %>% filter(SYMBOL %in% gene.list) %>% dplyr::select(SYMBOL,LEAD.VAR), all_of(c("LEAD.VAR")), sep=",")
  dee.variants$GROUPS <- "DEE"
  colnames(dee.variants) <- c("SYMBOL","HGVSP","GROUPS")
  
  gge.variants <- separate_rows(gge.filtered  %>% filter(SYMBOL %in% gene.list) %>% dplyr::select(SYMBOL,LEAD.VAR), all_of(c("LEAD.VAR")), sep=",")
  gge.variants$GROUPS <- "GGE"
  colnames(gge.variants) <- c("SYMBOL","HGVSP","GROUPS")
  
  nafe.variants <- separate_rows(nafe.filtered  %>% filter(SYMBOL %in% gene.list) %>% dplyr::select(SYMBOL,LEAD.VAR), all_of(c("LEAD.VAR")), sep=",")
  nafe.variants$GROUPS <- "NAFE"
  colnames(nafe.variants) <- c("SYMBOL","HGVSP","GROUPS")
  
  
  
  dee.variants <- dee.variants %>% mutate(
    SYMBOL = trimws(SYMBOL),
    HGVSP = trimws(HGVSP),
    GROUPS = trimws(GROUPS)
  )
  
  gge.variants <- gge.variants %>% mutate(
    SYMBOL = trimws(SYMBOL),
    HGVSP = trimws(HGVSP),
    GROUPS = trimws(GROUPS)
  )
  
  nafe.variants <- nafe.variants %>% mutate(
    SYMBOL = trimws(SYMBOL),
    HGVSP = trimws(HGVSP),
    GROUPS = trimws(GROUPS)
  )
  
  
  
  tmp <- info.expanded  %>% mutate(
    SYMBOL = trimws(SYMBOL),
    HGVSP = trimws(HGVSP),
    GROUPS = trimws(GROUPS)
  )
  
  
  dee.variants <- inner_join(dee.variants,tmp, by = c("SYMBOL","HGVSP","GROUPS")) 
  nafe.variants <- inner_join(nafe.variants,tmp, by = c("SYMBOL","HGVSP","GROUPS"))
  gge.variants <- inner_join(gge.variants,tmp, by = c("SYMBOL","HGVSP","GROUPS"))
  
  all.variants <- rbind(gge.variants,nafe.variants,dee.variants)
  
  
}

dee.filtered <- Leading_Modified$DEE_NONSYNONYMOUS %>% select(SYMBOL,HGVSP)
colnames(dee.filtered) <- c("SYMBOL","LEAD.VAR")

gge.filtered <- Leading_Modified$GGE_NONSYNONYMOUS %>% select(SYMBOL,HGVSP)
colnames(gge.filtered) <- c("SYMBOL","LEAD.VAR")

nafe.filtered <- Leading_Modified$NAFE_NONSYNONYMOUS %>% select(SYMBOL,HGVSP)
colnames(nafe.filtered) <- c("SYMBOL","LEAD.VAR")


var.list <- NULL
for (s in names(gene.set)) {
  var.list[[s]]  <- get.estimates(gene.set[[s]])
}


for (v in names(var.list)) {
  var.list[[v]] <- var.list[[v]] %>% mutate(PATH = v)
}

combined.var <- NULL
for (v in names(var.list)) {
  combined.var <- rbind(combined.var,var.list[[v]])
}

#source("src/packages.R")
combined.var.filtered <- combined.var[!duplicated(combined.var[,c('SYMBOL','HGVSP','GROUPS','PATH')]),] 

get_sig <- function(group,path){
  sig <- NULL
  sig <- fet.res %>% filter(GROUP == group, PATH == path) %>% select(SIG) %>% unlist
  return(sig)
}


combined.var.filtered$SIG <- mapply(get_sig,combined.var.filtered$GROUPS,combined.var.filtered$PATH)

p <- ggplot(combined.var.filtered, aes(x = PATH, y = ESTIMATES)) + 
  geom_boxplot(aes(fill = PATH),outlier.size = 1,outlier.shape = 1) + 
  geom_jitter(size = 1,shape = 1) +
  facet_wrap(~GROUPS, scales="free_x") +
  theme_bw() +
  xlab("Gene Set") + 
  ylab("Odds Ratio") + 
  coord_flip() +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(2), hjust = 0.5),
        axis.title = element_text(size = rel(2)),
        axis.text.x = element_text(angle = 90,hjust = 0.95,vjust = 0.2, size =rel(2)),
        axis.text.y = element_text(size = rel(2),face = "bold"))



png("graph/CUSTOM_ORA_BOX.png", width = 42, height = 21, units = 'cm', res = 300)
print(p)
dev.off()




p.rescaled <- ggplot(combined.var.filtered, aes(x = PATH, y = ESTIMATES)) + 
  geom_boxplot(aes(fill = PATH),outlier.size = 1,outlier.shape = 1) + 
  geom_jitter(size = 1,shape = 1) +
  facet_wrap(~GROUPS) +
  theme_bw() +
  xlab("Gene Set") + 
  ylab("Odds Ratio") + 
  coord_flip() +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(2), hjust = 0.5),
        axis.title = element_text(size = rel(2)),
        axis.text.x = element_text(angle = 90,hjust = 0.95,vjust = 0.2, size =rel(2)),
        axis.text.y = element_text(size = rel(2),face = "bold"))


png("CUSTOM_ORA_BOX_RESCALED.png", width = 42, height = 21, units = 'cm', res = 300)
print(p.rescaled)
dev.off()
