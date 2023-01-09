

# Use to separate multiple values in a single cell ------------------------
expandInfo = function (columns_selected, dataset){
  require(tidyr)
  df <- separate_rows(dataset, all_of(columns_selected), sep=",")
}



# Use to combine multiple values in a single cell -------------------------
collateInfo = function (columns_selected, dataset){
  for (cols in columns_selected) {
    print(cols)
    collated <- as.matrix(dataset[[cols]]) %>% apply(.,1,get_count)
    print(collated)
    dataset[[cols]] <- collated  
    
  }
  return(dataset)
}


# For getting the sum of a vector listed in a single cell -----------------
get_count = function(count_in_string){
  count_vector <- as.numeric(unlist(stri_split(count_in_string,fixed=',')))
  return(sum(count_vector))
}



# For extracting only the needed columns ----------------------------------
extractInfo = function (columns_selected,dataset){
  dataset <- dataset %>% select(all_of(columns_selected))
  return(dataset)
}

buildGDB = function(gene_vector){
  gene.db <- getBM(filters= "ensembl_gene_id", 
                        attributes = c("ensembl_gene_id","hgnc_symbol"),
                        values = gene_vector ,
                        mart = useDataset("hsapiens_gene_ensembl", useMart("ensembl")))
  colnames(gene.db ) <- c("GENE_ID","SYMBOL")
  return(gene.db)
}


bindGDB = function(dataset, genedb){
  dataset <- inner_join(dataset,genedb,by="GENE_ID")  
}


calculateAF <- function(dataset){
  
  dataset <- dataset %>% mutate(AF_CASE = AC_CASE/AN_CASE,
                        AF_CTRL = AC_CTRL/AN_CTRL)
  return(dataset)
}


calculateOR <- function(dataset){
  
  AC_CASE <- dataset$AC_CASE %>% as.numeric
  AN_CASE <- dataset$AN_CASE %>% as.numeric
  AC_CTRL <- dataset$AC_CTRL %>% as.numeric
  AN_CTRL <- dataset$AN_CTRL %>% as.numeric
  
  odds <- fmsb::oddsratio(AC_CASE,(AN_CASE-AC_CASE),AC_CTRL,(AN_CTRL-AC_CTRL))
  oddresults <- data.frame(pval = odds$p.value,estimate = odds$estimate)
  colnames(oddresults) <- c('PVAL','ESTIMATES')
  return(oddresults)
  
}


# Classify each variants as risk or protective factor ---------------------
classify <- function(dataset){
  dataset <- dataset %>% filter(ESTIMATES != Inf,PVAL != Inf)
  dataset <- dataset %>% mutate(FACTORS = 
                                                case_when(
                                                  ESTIMATES > 1 ~ "RISK",
                                                  ESTIMATES < 1 ~ "PROTECTIVE",
                                                  TRUE ~ "NONE"
                                                ))
  return(dataset)
}




# Collapsed the variants per genes ----------------------------------------
collapsed_variants <- function(variants = NULL, 
                               mutation = "SYNONYMOUS", 
                               group = "DEE",
                               ac_thr = 0){
  if(!is.na(group)){
    variants_filtered <- variants %>% filter(MUTATION == mutation, 
                                             AC_CTRL >= ac_thr,
                                             GROUPS == group)
  } else {
    
    variants_filtered <- variants %>% filter(MUTATION == mutation, 
                                             AC_CTRL >= ac_thr)
    
  }

  variants_collapsed <- table(variants_filtered$SYMBOL, variants_filtered$FACTORS) %>% as.data.frame.matrix
  variants_collapsed <- cbind(variants_collapsed,
                              calc_risk(variants_collapsed$RISK, 
                                        variants_collapsed$PROTECTIVE,
                                        rownames(variants_collapsed)))
  variants_collapsed$RISK.PADJ <- p.adjust(variants_collapsed$RISK.PVAL, method = "BH")
  
  
  
  return(variants_collapsed)
}


calc_risk <- function(case, control, genes){
  
  AC_CASE <- case %>% as.numeric %>% unlist
  AC_CTRL <- control %>% as.numeric %>% unlist
  Alleles <- genes
  riskdata <- data.frame(AC_CASE,AC_CTRL) %>% as.matrix
  dimnames(riskdata) <- list(ALLELES = Alleles, Outcome = c("Case","Control"))
  risks <- riskratio(riskdata, rev = "c")
  
  # relation 
  riskresults <- data.frame(pval = risks$p.value[,"chi.square"]
                            ,estimate = risks$measure[,"estimate"])
  
  colnames(riskresults) <- c('RISK.PVAL','RISK.ESTIMATES')
  return(riskresults)
}






'%notin%' <-  Negate('%in%')


qqplotter_uniform<-function(pvalues, 
                            should.thin=T, thin.obs.places=2, thin.exp.places=2, 
                            title = "QQ Plot",
                            remarks = "QQ Plot",
                            xlab=expression(paste("Expected (",-log[10], " p-value)")),
                            ylab=expression(paste("Observed (",-log[10], " p-value)")), 
                            draw.conf=TRUE, conf.points=1000, conf.col="lightgray", conf.alpha=.05,
                            already.transformed=FALSE, pch=20, aspect="iso", prepanel=prepanel.qqunif,
                            par.settings=list(superpose.symbol=list(pch=pch)), ...) {
  
  
  
  pvalues <- pvalues[!is.na(pvalues)]
  grp<-NULL
  n<-1
  exp.x<-c()
  if(is.list(pvalues)) {
    nn<-sapply(pvalues, length)
    rs<-cumsum(nn)
    re<-rs-nn+1
    n<-min(nn)
    if (!is.null(names(pvalues))) {
      grp=factor(rep(names(pvalues), nn), levels=names(pvalues))
      names(pvalues)<-NULL
    } else {
      grp=factor(rep(1:length(pvalues), nn))
    }
    pvo<-pvalues
    pvalues<-numeric(sum(nn))
    exp.x<-numeric(sum(nn))
    for(i in 1:length(pvo)) {
      if (!already.transformed) {
        pvalues[rs[i]:re[i]] <- -log10(pvo[[i]])
        exp.x[rs[i]:re[i]] <- -log10((rank(pvo[[i]], ties.method="first")-.5)/nn[i])
      } else {
        pvalues[rs[i]:re[i]] <- pvo[[i]]
        exp.x[rs[i]:re[i]] <- -log10((nn[i]+1-rank(pvo[[i]], ties.method="first")-.5)/(nn[i]+1))
      }
    }
  } else {
    n <- length(pvalues)+1
    if (!already.transformed) {
      exp.x <- -log10((rank(pvalues, ties.method="first")-.5)/n)
      pvalues <- -log10(pvalues)
    } else {
      exp.x <- -log10((n-rank(pvalues, ties.method="first")-.5)/n)
    }
  }
  
  
  #this is a helper function to draw the confidence interval
  panel.qqconf<-function(n, conf.points=1000, conf.col="gray", conf.alpha=.05, ...) {
    require(grid)
    conf.points = min(conf.points, n-1);
    mpts<-matrix(nrow=conf.points*2, ncol=2)
    for(i in seq(from=1, to=conf.points)) {
      mpts[i,1]<- -log10((i-.5)/n)
      mpts[i,2]<- -log10(qbeta(1-conf.alpha/2, i, n-i))
      mpts[conf.points*2+1-i,1]<- -log10((i-.5)/n)
      mpts[conf.points*2+1-i,2]<- -log10(qbeta(conf.alpha/2, i, n-i))
    }
    grid.polygon(x=mpts[,1],y=mpts[,2], gp=gpar(fill=conf.col, lty=0), default.units="native")
  }
  
  
  
  
  prepanel.qqunif= function(x,y,...) {
    A = list()
    A$xlim = range(x, y)*1.02
    A$xlim[1]=0
    A$ylim = A$xlim
    return(A)
  }  
  #draw the plot
  p_xlab = xlab
  p_ylab = ylab
  p_title = title
  p_remarks = remarks
  
  #hist(exp.x, main = "Theoretical P Values")
  #hist(pvalues, main = "Actual P Values")
  xyplot(pvalues~exp.x, groups=grp, xlab=p_xlab, ylab=p_ylab, aspect=aspect,main = p_title,
         sub = p_remarks,
         prepanel=prepanel, scales=list(axs="i"), pch=pch,
         panel = function(x, y, ...) {
           if (draw.conf) {
             panel.qqconf(n, conf.points=conf.points, 
                          conf.col=conf.col, conf.alpha=conf.alpha)
           };
           
           
           # take the genes with the lowes pvlaues
           panel.xyplot(x,y, ...);
           #xlabels <- pvalues %>% tail(.,10)
           #panel.text(xlabels,xlabels,labels = names(xlabels), pos = 3)
           panel.abline(0,1);
         }, par.settings=par.settings, ...
  )
  
  
  
}


qqplotter_lambda <- function(values){
  
  values <- values[!is.na(values)]  
  chisq <- qchisq(values, 1, lower.tail = F)
  lambda <- median(chisq) / qchisq(0.5,1)
  return(lambda)
  
}


# test if non synonymous risk factor variants are over presented in the gene set of
# non synonymous variants
hypergeometric_test <- function(x){

  ns.variants.risk <- x[1]
  ns.all.variants <- x[2]
  s.variants.risk <- x[3]
  s.all.variants <- x[4]
  # Build Contingency Table
  dat <- data.frame(
      "non_synonymous" = c(ns.variants.risk, ns.all.variants),
      "synonymous" = c(s.variants.risk, s.all.variants),
      row.names = c("RISK", "NOTRISK"),
      stringsAsFactors = FALSE
    )
  colnames(dat) <- c("NS", "S")
    
  result <- fisher.test(dat, alternative = "greater")
    
  return (result$p.value)
  
  
}



# reannotation of variants ------------------------------------------------
annotate <- function(dataset){
  dataset <- dataset %>% filter(ESTIMATES != Inf,PVAL != Inf)
  dataset <- dataset %>% mutate(ANNOTATIONS = 
                                  case_when(
                                    CONSEQUENCE == "loss of function" ~ "PTV",
                                    CONSEQUENCE == "synonymous" ~ "SYNONYMOUS",
                                    CONSEQUENCE %in% c("missense","splice region","inframe indel") & POLYPHEN == "benign" ~ "BENIGN_MISSENSE",
                                    CONSEQUENCE %in% c("missense","splice region","inframe indel") & POLYPHEN %in% c("probably_damaging","possibly_damaging")  & MPC < 2~ "DAMAGING_MISSENSE",
                                    CONSEQUENCE %in% c("missense","splice region","inframe indel") & POLYPHEN %in% c("probably_damaging","possibly_damaging") & MPC >= 2 ~ "DAMAGING_MISSENSE_MPC",
                                    TRUE ~ "NONE"
                                  ))
  return(dataset)
}



perform_fisher <- function(x){
  
  case_exposed <- x[1]
  case_unexposed <- x[2]
  ctrl_exposed <- x[3]
  ctrl_unexposed <- x[4]
  print(x[5])
  # Build Contingency Table
  dat <- data.frame(
    "allele_yes" = c(case_exposed, case_unexposed),
    "allele_no" = c(ctrl_exposed, ctrl_unexposed),
    row.names = c("WAL", "WOAL"),
    stringsAsFactors = FALSE
  )
  colnames(dat) <- c("CASE", "CTRL")
  
  fet <- fisher.test(dat, conf.int = TRUE)
  fetres <- as.data.frame(pval = fet$p.value, or = fet$estimate, fet$conf.int)
  
  return (fetres)
}



get_overlap <- function(x,y){
  z <- intersect(x,y)
  return(length(z))
}


enrich.fet <- function(pathway,overlap,signature, backpop = 26588 ){
  
  contable <- matrix(c(
    dg = length(signature) - overlap,
    dr = overlap,
    ng = backpop - length(signature) - length(pathway) + overlap,
    nr = length(pathway) - overlap),2,2,dimnames = list(c("GREEN","RED"),c("DRAWN","NOT DRAWN")))
  
  print(contable)
  
  fet <- fisher.test(contable,alt="less")
  return(fet)
}


enrichment_test <-  function(overlap, hf, pop, lf,tail){
  p.val <- phyper( q=overlap-1, # number of variants common in control and case - 1          
                   m=hf, # number of all high frequency variants across the population         
                   n=pop-hf, # number of non high frequency variants across the population
                   k=lf, # number of controls with low frequency       
                   lower.tail=tail)
  score <- (overlap/lf)/(hf/pop)
  
  return(data.frame(p = p.val, es = score))
}









# ggplot(sig_variants %>% filter(!is.na(CONSEQUENCE.GROUPED)), 
#        aes(ODD.ESTIMATES, colour = CONSEQUENCE.GROUPED)) +
#   stat_ecdf() +
#   ggtitle("ECDF of OR") +
#   xlab("Odd Ratio") + 
#   ylab("Cumulative Sum") + 
#   theme(plot.title = element_text(size = rel(1.5), hjust = 0.5), 
#         axis.title = element_text(size = rel(1.25)),
#         axis.text.x=element_text(angle=90,hjust=0.95,vjust=0.2)) +
#   labs(color = "CONSEQUENCE")
# 
# 
# 
# 
# ggplot(sig_variants %>% filter(!is.na(CONSEQUENCE.GROUPED)), aes(x = ODD.ESTIMATES, fill = CONSEQUENCE.GROUPED)) +                       # Draw overlaying histogram
#   geom_density()+ 
#   scale_x_log10() +
#   ggtitle("Density Plot of OR") +
#   xlab("Odd Ratio") + 
#   ylab("Density") + 
#   theme(plot.title = element_text(size = rel(1.5), hjust = 0.5), 
#         axis.title = element_text(size = rel(1.25)),
#         axis.text.x=element_text(angle=90,hjust=0.95,vjust=0.2)) +
#   labs(color = "CONSEQUENCE")
# 
# 
# ggplot(sig_variants %>% filter(!is.na(CONSEQUENCE.GROUPED)), aes(x = ODD.ESTIMATES, fill = CONSEQUENCE.GROUPED)) +                       # Draw overlaying histogram
#   geom_histogram() +
#   scale_x_log10() +
#   ggtitle("Histogram Plot of OR") +
#   xlab("Odd Ratio") + 
#   ylab("Frequency") + 
#   theme(plot.title = element_text(size = rel(1.5), hjust = 0.5), 
#         axis.title = element_text(size = rel(1.25)),
#         axis.text.x=element_text(angle=90,hjust=0.95,vjust=0.2)) +
#   labs(color = "CONSEQUENCE")
# 

enrichGO <-  function(gene_vector,category = "BP",thr = 1, score = "std"){
  gene_df <- data.frame(SYMBOL = names(gene_vector),PADJ = gene_vector) 
  
  gene_names <- bitr(names(gene_vector), fromType = "SYMBOL",
                     toType = c("ENTREZID"),
                     OrgDb = org.Hs.eg.db) 
  
  gene_set <- inner_join(gene_df,gene_names,by="SYMBOL") %>% dplyr::select(c(ENTREZID,PADJ)) 
    
  genes <- gene_set$ENTREZID
  gene_list <-gene_set$PADJ
  names(gene_list) <- genes
  enrichment <- gseGO(geneList = sort(gene_list, decreasing = TRUE),
                      OrgDb        = org.Hs.eg.db,
                      ont          = category,
                      minGSSize    = 100,
                      maxGSSize    = 500,
                      pvalueCutoff = thr,
                      pAdjustMethod = "fdr",
                      verbose      = FALSE,
                      eps = 0,
                      scoreType = score)
  return(setReadable(enrichment,"org.Hs.eg.db","ENTREZID")@result)
}


enrichKP <-  function(gene_vector,thr = 1){
  gene_set <- as.data.frame(gene_vector) %>% rownames_to_column("SYMBOL")
  colnames(gene_set) <- c("SYMBOL","PADJ")
  gene_names <- bitr(names(gene_vector), fromType = "SYMBOL",
                     toType = c("ENTREZID"),
                     OrgDb = org.Hs.eg.db)
  #print(gene_names)
  gene_set <- inner_join(gene_set,gene_names,by="SYMBOL") %>% dplyr::select(c(ENTREZID,PADJ)) 
  genes <- gene_set$ENTREZID
  gene_list <-gene_set$PADJ
  names(gene_list) <- genes
  enrichment <- enrichKEGG(gene         = genes,
                           organism     = 'hsa',
                           pvalueCutoff = thr)
}


enrichDisease <- function(gene_vector){
  
  geneList <- bitr(names(gene_vector), fromType = "SYMBOL",
                     toType = c("ENTREZID"),
                     OrgDb = org.Hs.eg.db)
  print(geneList)
  do <- enrichDO(gene          = geneList %>% unlist,
                 ont           = "DO",
                 pvalueCutoff  = 0.05,
                 pAdjustMethod = "BH",
                 universe      = names(gene_vector),
                 minGSSize     = 5,
                 maxGSSize     = 500,
                 qvalueCutoff  = 0.05,
                 readable      = FALSE) 
}


vectorize <- function(genes,values){
  v <- values
  names(v) <- genes
  return(sort(v, decreasing = TRUE))
}



# Enrichment Analysis by unweighted KS ------------------------------------

gsea.perm <- function(score,vector,nperm = 1000, sample.size, w = 0){
  
  set.seed(123)
  
  vector <- sort(vector,decreasing = TRUE)
  
  es.perm <- NULL
  nes <- 0
  p.val <- 0
  zs <- 0
  vector.length <- length(vector)
  
  for (n in c(1:nperm)) {
    
    # Compute ES
    counter <- c(rep(FALSE,vector.length))
    
    col <- randperm(vector.length,sample.size)
    counter[col] <- TRUE
    hits <- cumsum(counter*(vector^w))
    miss <- cumsum((1-counter)*(vector^w))
    hits.normalized <- hits/length(counter[counter == TRUE])
    miss.normalized <- miss/length(counter[counter == FALSE])
    es <- max(hits.normalized - miss.normalized)
    es.perm <- c(es.perm,es)
  } 
  
  p.val <- length(es.perm[es.perm >= score])/nperm
  zs <- (score - mean(es.perm))/std(es.perm)
  nes <- score/abs(mean(es.perm[es.perm > 0]))
  
  
  return(data.frame(PVAL = p.val, NES = nes, Z = zs, 
                    PERM = es.perm  %>% toString))
}

# ranked list is the variant
gsea.compute <- function(ranked.list, gene.set, nperm = 1000, w = 0){
  
  ranked.list <- sort(ranked.list,decreasing = TRUE)
  
  counter <- names(ranked.list) %in% (gene.set %>% unlist)
  
  enriched.genes <- names(ranked.list[counter])
  
  hits <- cumsum(counter*(ranked.list^w))
  miss <- cumsum((1-counter)*(ranked.list^w))
  hits.normalized <- hits/length(counter[counter == TRUE])
  miss.normalized <- miss/length(counter[counter == FALSE])
  
  es.all <- hits.normalized - miss.normalized
  es.index <- which.max(es.all)
  es <- max(hits.normalized - miss.normalized)
  
  if(length(es.index) != 0){
    lead.genes <- names(ranked.list[1:es.index])
    leading.edge <- names(ranked.list)[names(ranked.list) %in% enriched.genes & names(ranked.list) %in% lead.genes] %>% toString  
  } else {
    leading.edge <- NA
  }
  
  perm.results <- gsea.perm(es, ranked.list, nperm, length(counter[counter == TRUE]), w = w)
  
  es.scores <- as.numeric(es.all)

  
  return(cbind(data.frame(ES = es,Lead = leading.edge, 
                          Scores = es.all %>% toString,
                          VI = names(ranked.list) %>% toString,
                          VB = gene.set %>% unlist %>% toString),
               perm.results)) 
}


gsea.fdr <- function(NES.score,ES.perm,NES.obs, nperm = 1000){
  
  NES.perm <- ES.perm[ES.perm > 0]/abs(mean(ES.perm[ES.perm > 0]))
  perm <- length(NES.perm[NES.perm >= NES.score & NES.perm >= 0])/length(NES.perm)
  
  obs <- length(NES.obs[NES.obs >= NES.score & NES.obs >= 0])/length(NES.obs)
  
  q <- perm/obs
  
  return(q)
}







enrichReactome <-  function(gene_vector,thr = 1){
  gene_set <- as.data.frame(gene_vector) %>% rownames_to_column("SYMBOL")
  colnames(gene_set) <- c("SYMBOL","PADJ")
  gene_names <- bitr(names(gene_vector), fromType = "SYMBOL",
                     toType = c("ENTREZID"),
                     OrgDb = org.Hs.eg.db)
  #print(gene_names)
  gene_set <- inner_join(gene_set,gene_names,by="SYMBOL") %>% dplyr::select(c(ENTREZID,PADJ)) 
  genes <- gene_set$ENTREZID
  gene_list <-gene_set$PADJ
  names(gene_list) <- genes
  enrichment <- enrichPathway(gene = genes, pvalueCutoff = thr, readable = TRUE)
  return(enrichment)
}



#  Custom Enrichment ------------------------------------------------------
initializeGSEA = function(name = "Hallmark",pathway.db = "Data/h.all.v7.4.symbols.gmt" ){
  db <- gmtPathways(pathway.db)
  return(db)
}

performGSEA = function(name = "Hallmark",perm =10000, db, genevector){
  result <- fgsea(pathways = db, stats=sort(genevector, decreasing = T), 
                  nPermSimple = perm, scoreType = "pos", gseaParam = 0)
  result <- result %>% arrange(desc(NES))
}




# For Over representation -------------------------------------------------
calc.score <- function(hits,test,pathway,population = 23000){
  score <- 0
  score <- (hits/test)/(pathway/population)
  return(score)
}

znorm <- function(x, mean,sd){
  z <- 0
  z <- (log10(x+0.001) - mean)/sd
  return(z)
}


znorm2 <- function(x, mean,sd){
  z <- 0
  z <- (x- mean)/sd
  return(z)
}


calc.nes <- function(score, absolute.mean){
  nes <- 0
  nes <- score/absolute.mean
  return(nes)
}


