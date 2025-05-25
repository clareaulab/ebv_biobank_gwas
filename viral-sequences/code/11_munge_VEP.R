library(data.table)
library(dplyr)

process_one <- function(one){
  dt <- fread(one, skip = "#", header = TRUE) %>% data.frame()
  head(dt)
  dt %>% filter(Consequence == "missense_variant")
  dt$Location <- as.numeric(gsub("NC_007605.1:", "", dt$Location ))
  dt$raw <- stringr::str_split_fixed(dt$X.Uploaded_variation, "_", 4)[,4]
  dt <- dt %>% mutate(REF = substr(raw, 1, 1), ALT = substr(raw, 3, 3))
  dt[,c("Location","Feature","Consequence","Amino_acids", "Codons", "REF", "ALT")] 
}
dtf <- lapply(list.files("../data/veps/", full.names = TRUE), process_one) %>% rbindlist() %>% data.frame()
saveRDS(dtf, file = "../output/full_EBV_annotation.rds")
