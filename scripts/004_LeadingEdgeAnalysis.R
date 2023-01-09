

# Helper Functions --------------------------------------------------------

get_leading <- function(VSEA_LIST){
  
  VSEA_LEADING <- NULL
  for (item in names(VSEA_LIST)) {
    
    VSEA_LEADING[[item]] <- core.enrichment(VSEA_LIST[[item]],item %>% gsub("_.*", "", .))  
    
  }
  
  VSEA_LEADING_LIST <- list(
    DEE_SYNONYMOUS = VSEA_LEADING$DEE_SYNONYMOUS,
    DEE_NONSYNONYMOUS = VSEA_LEADING$DEE_NONSYNONYMOUS,
    NAFE_SYNONYMOUS = VSEA_LEADING$NAFE_SYNONYMOUS,
    NAFE_NONSYNONYMOUS = VSEA_LEADING$NAFE_NONSYNONYMOUS,
    GGE_SYNONYMOUS = VSEA_LEADING$GGE_SYNONYMOUS,
    GGE_NONSYNONYMOUS = VSEA_LEADING$GGE_NONSYNONYMOUS  
  )
  return(VSEA_LEADING_LIST)
}




# Leading Edge Function ---------------------------------------------------
core.enrichment <- function(dataset,epitype){
  VSEA_LEADING <- dataset %>% 
    filter(FDR <= 0.05, Z >= 1.96) %>% 
    select(TERM,Lead)
  
  print(epitype)
  print(VSEA_LEADING)
  VSEA_LEADING <- separate_rows(VSEA_LEADING, all_of(c("Lead")), sep=",")
  
  
  colnames(VSEA_LEADING) <- c("SYMBOL","HGVSP")
  
  VSEA_LEADING <- VSEA_LEADING %>% mutate(HGVSP = trimws(HGVSP))
  INFO_LEADING <- info.expanded %>% mutate(HGVSP = trimws(HGVSP))
  
  VSEA_LEADING <- inner_join(VSEA_LEADING, INFO_LEADING, by = c("SYMBOL","HGVSP")) %>% 
    filter(GROUPS == epitype)
  

  return(VSEA_LEADING)
}



# Usage -------------------------------------------------------------------

Leading_Modified <- get_leading(VSEA_List_MODIFIED)


for (g in names(Leading_Modified)) {
  print(Leading_Modified[[g]] %>% nrow)
}

openxlsx::write.xlsx(Leading_Modified,"VSEA_Leading_Modified.xlsx")

