proc <- function(show = TRUE){
  
  cat(paste0("___Computing Enrichment in ",EPI.type," ",MUTATION,"___\n"))
  cat("_________________Parameters__________________\n")
  cat(paste0("Epilepys Group: ",EPI.type,"\n"))
  cat(paste0("Mutation:       ",MUTATION,"\n"))
  cat(paste0("CADD Threshold: ",CADD.threshold,"\n"))
  cat(paste0("AF Threshold:   ",AF.threshold,"\n"))
  cat(paste0("Enrichment mode: ",VSEA.mode,"\n"))
  
}

# MODIFIED ---------------------------------------------------------------

EPI.type <- "DEE"
MUTATION <- "NONSYNONYMOUS"
CADD.threshold <- 20
AF.threshold <- 0.001
VSEA.mode <- "modified"

proc(TRUE)

VSEA_DEE_NONSYNONYMOUS <- vsea(AF.threshold,CADD.threshold,EPI.type,MUTATION,mode = VSEA.mode, gene.list = qualified.genes)


EPI.type <- "DEE"
MUTATION <- "SYNONYMOUS"
CADD.threshold <- 20
AF.threshold <- 0.001

proc(TRUE)

VSEA_DEE_SYNONYMOUS <- vsea(AF.threshold,CADD.threshold,EPI.type,MUTATION,mode = VSEA.mode, gene.list = qualified.genes)


EPI.type <- "GGE"
MUTATION <- "NONSYNONYMOUS"
CADD.threshold <- 20
AF.threshold <- 0.001


proc(TRUE)

VSEA_GGE_NONSYNONYMOUS <- vsea(AF.threshold,CADD.threshold,EPI.type,MUTATION,mode = VSEA.mode, gene.list = qualified.genes)


EPI.type <- "GGE"
MUTATION <- "SYNONYMOUS"
CADD.threshold <- 20
AF.threshold <- 0.001


proc(TRUE)


VSEA_GGE_SYNONYMOUS <- vsea(AF.threshold,CADD.threshold,EPI.type,MUTATION,mode = VSEA.mode, gene.list = qualified.genes)

EPI.type <- "NAFE"
MUTATION <- "NONSYNONYMOUS"
CADD.threshold <- 20
AF.threshold <- 0.001

proc(TRUE)

VSEA_NAFE_NONSYNONYMOUS <- vsea(AF.threshold,CADD.threshold,EPI.type,MUTATION,mode = VSEA.mode, gene.list = qualified.genes)



EPI.type <- "NAFE"
MUTATION <- "SYNONYMOUS"
CADD.threshold <- 20
AF.threshold <- 0.001

proc(TRUE)

VSEA_NAFE_SYNONYMOUS <- vsea(AF.threshold,CADD.threshold,EPI.type,MUTATION,mode = VSEA.mode, gene.list = qualified.genes)


VSEA_List_MODIFIED <- list(
  DEE_NONSYNONYMOUS = VSEA_DEE_NONSYNONYMOUS %>% 
    filter(Z >= 1.96,FDR <= 0.05) %>% select(TERM,ES,PVAL,NES,Z,FDR, Lead
),
  DEE_SYNONYMOUS = VSEA_DEE_SYNONYMOUS %>% 
    filter(Z >= 1.96,FDR <= 0.05)%>% select(TERM,ES,PVAL,NES,Z,FDR, Lead
    ),
  NAFE_NONSYNONYMOUS = VSEA_NAFE_NONSYNONYMOUS %>% 
    filter(Z >= 1.96,FDR <= 0.05)%>% select(TERM,ES,PVAL,NES,Z,FDR, Lead
    ),
  NAFE_SYNONYMOUS = VSEA_NAFE_SYNONYMOUS %>% 
    filter(Z >= 1.96,FDR <= 0.05)%>% select(TERM,ES,PVAL,NES,Z,FDR, Lead
    ),
  GGE_NONSYNONYMOUS = VSEA_GGE_NONSYNONYMOUS %>% 
    filter(Z >= 1.96,FDR <= 0.05)%>% select(TERM,ES,PVAL,NES,Z,FDR, Lead
    ),
  GGE_SYNONYMOUS = VSEA_GGE_SYNONYMOUS %>% 
    filter(Z >= 1.96,FDR <= 0.05) %>% select(TERM,ES,PVAL,NES,Z,FDR, Lead
    )
)

openxlsx::write.xlsx(VSEA_List_MODIFIED,"Supplementary_Table_1_VSEA.xlsx")



# Counting of Enriched genes ----------------------------------------------



# Z threshold
cat(paste0("DEE NON SYNONYMOUS:",VSEA_DEE_NONSYNONYMOUS %>% 
             filter(Z >= 1.96,FDR <= 0.05) %>% nrow), "\n")
cat(paste0("DEE SYNONYMOUS:",VSEA_DEE_SYNONYMOUS %>% 
             filter(Z >= 1.96,FDR <= 0.05) %>% nrow), "\n")

cat(paste0("GGE NON SYNONYMOUS:",VSEA_GGE_NONSYNONYMOUS %>% 
             filter(Z >= 1.96,FDR <= 0.05) %>% nrow), "\n")
cat(paste0("GGE SYNONYMOUS:",VSEA_GGE_SYNONYMOUS %>% 
             filter(Z >= 1.96,FDR <= 0.05) %>% nrow), "\n")

cat(paste0("NAFE NON SYNONYMOUS:",VSEA_NAFE_NONSYNONYMOUS %>% 
             filter(Z >= 1.96,FDR <= 0.05) %>% nrow), "\n")
cat(paste0("NAFE SYNONYMOUS:",VSEA_NAFE_SYNONYMOUS %>% 
             filter(Z >= 1.96,FDR <= 0.05) %>% nrow), "\n")

