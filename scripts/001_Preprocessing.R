# Load packages required for the analysis ---------------------------------
# source("src/packages.R")
# source("src/helper.R")

# Load the VCF ------------------------------------------------------------
if (!file.exists("results/epi25_info.csv")) {
  vcf <- read.vcf("data/epi25.vcf.bgz",split.info = TRUE, split.samples = TRUE)
  info <- vcf$vcf
  colnames(info) <- casefold(colnames(info),upper = TRUE)
  write.csv(info,"results/epi25_info.csv")  
} else {
  info <- vroom("results/epi25_info.csv") %>% dplyr::select(-c("...1")) 
}



# Preprocess the information from the VCF file ----------------------------
if(!file.exists("results/epi25_expanded.csv")){
  info.expanded <- expandInfo(c("GROUPS","AC_CASE","AN_CASE","AC_CTRL","AN_CTRL"),
                              info)
  write.csv(info.expanded,"results/epi25_expanded.csv")
} else {
  info.expanded <- vroom("results/epi25_expanded.csv") %>% dplyr::select(-c("...1")) 
}

if(!file.exists("results/epi25_collated.csv")){
  info.collated <- collateInfo(c("AC_CASE","AN_CASE","AC_CTRL","AN_CTRL"), info)  
  write.csv(info.collated,"results/epi25_collated.csv")
} else {
  info.collated <- vroom("results/epi25_collated.csv") %>% dplyr::select(-c("...1"))
}

selected_cols <- c("CHROM","POS","REF","ALT","GENE_ID","CONSEQUENCE", "HGVSP","HGVSC",
                   "VARIANT.ID","CADD","POLYPHEN","MPC","GROUPS","AC_CASE","AN_CASE","AC_CTRL","AN_CTRL")

info.expanded <- extractInfo(selected_cols,info.expanded)
info.collated <- extractInfo(selected_cols,info.collated)


# Acquire the gene 'SYMBOL' list  -----------------------------------------
if(!file.exists("results/genes.csv")){
  
  gene.list <- buildGDB(info$GENE_ID)
  write.csv(gene.list,"results/genes.csv")
  
} else {
  gene.list <- vroom("results/genes.csv") %>% dplyr::select(-c("...1"))
}

# Bind the genes database to the dataset ----------------------------------
info.expanded <- bindGDB(info.expanded, gene.list)
info.collated <- bindGDB(info.collated, gene.list)



# Modify the dataset ------------------------------------------------------

info.expanded <- info.expanded %>% 
  mutate(CADD = as.numeric(CADD),
         AC_CASE = as.numeric(AC_CASE),
         AN_CASE = as.numeric(AN_CASE),
         AC_CTRL = as.numeric(AC_CTRL),
         AN_CTRL = as.numeric(AN_CTRL)) %>%
  dplyr::select(SYMBOL,GENE_ID, everything()) %>%
  filter(CONSEQUENCE %in% c('missense','synonymous','loss of function','splice region'),
         AC_CTRL > 0,
         AC_CASE > 0, 
         HGVSP != '.')


info.collated <- info.collated %>% 
  mutate(CADD = as.numeric(CADD),
         AC_CASE = as.numeric(AC_CASE),
         AN_CASE = as.numeric(AN_CASE),
         AC_CTRL = as.numeric(AC_CTRL),
         AN_CTRL = as.numeric(AN_CTRL)) %>%
  dplyr::select(SYMBOL,GENE_ID, everything()) %>%
  filter(CONSEQUENCE %in% c('missense','synonymous','loss of function','splice region'),
         AC_CTRL > 0,
         AC_CASE > 0, 
         HGVSP != '.')


# Calculate Allele Frequency ----------------------------------------------
info.expanded <- calculateAF(info.expanded)
info.collated <- calculateAF(info.collated)


# Calculate Odds Ratio ----------------------------------------------------
info.expanded <- cbind(info.expanded,calculateOR(info.expanded))
info.collated <- cbind(info.collated,calculateOR(info.collated))


# Classify the variants as protective or risk factor ----------------------
info.expanded <- classify(info.expanded)
info.collated <- classify(info.collated)


# Add groups based from CONSEQUENCE ---------------------------------------
info.expanded <- info.expanded %>% mutate(MUTATION =
                                    case_when(
                                      CONSEQUENCE == "synonymous" ~ "SYNONYMOUS",
                                      TRUE ~ "NONSYNONYMOUS"
                                      ))

info.collated <- info.collated %>% mutate(MUTATION =
                                            case_when(
                                              CONSEQUENCE == "synonymous" ~ "SYNONYMOUS",
                                              TRUE ~ "NONSYNONYMOUS"
                                            ))



# Add annotation based from CONSEQUENCE ---------------------------------------
info.expanded <- annotate(info.expanded)
info.collated <- annotate(info.collated)


# Read the GO database ----------------------------------------------------
GO.database <- readRDS("data/gene_ontology_db.Rds")



# Get brain genes ---------------------------------------------------------
brain.genes <- vroom("data/brain_gtex_all.tsv") %>% 
  select(Gene) %>% unlist


brain.genes.df <- vroom("data/brain_gtex_all.tsv") 


