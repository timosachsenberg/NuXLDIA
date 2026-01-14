library(tidyverse)
library(tidyr)
library(dplyr) 
library(stringr)
library(base)

#NuXl output for pep is different, no peak annotations  

#ALL FOR UV 
setwd("D:\\ECOLI_2025\\Ecoli_UVXL_DDA_singleshots_2026\\003-TextExporter-out")

file_list_UV <- list.files(pattern = "*.peptides\\.unknown$", full.names = TRUE)


# Read all files and store them in a list of data frames
df_list <- lapply(file_list_UV, read_tsv)

# Combine all data frames into one
pep_brp <- do.call(rbind, df_list)

new_colnames <- c(
  "RT" = 1,
  "CCS" = 29,
  "nuxl_score" = 68,
  "peak_ann" = 16
)


# Rename columns using the named vector
colnames(pep_brp)[new_colnames] <- names(new_colnames)

pep_brp <- pep_brp %>% mutate(IM = CCS*2.6867811*(10^25)/(charge*3/16*10000*(10^20)*((2*pi/(((mz*charge)*28.0134/(mz*charge+28.0134)*1.66053904*(10^(-27)))*1.38064852*(10^(-23))*(305)))^0.5)*1.6021766208*(10^(-19))))



pep_brpj <- pep_brp %>%  mutate (id = paste(pep_brp$sequence,pep_brp$charge,sep="_"))

setwd("D:\\ECOLI_2025\\Ecoli_UVXL_DDA_singleshots_2026\\library_shots")

pep_brpj$peak_ann <-  gsub ('"','',pep_brpj$peak_ann, fixed=T)

pep_brpsel <- pep_brpj %>% dplyr::select (nuxl_score,score, id,charge,RT,IM,sequence,accessions,peak_ann,mz)

pep_brpsel$score <- as.numeric(pep_brpsel$score)
#select lines with minimal score
pep_brpsel_min <- pep_brpsel %>% 
  group_by(id) %>% 
  slice_min(score, n = 1)
#exclude redundancy 
pep_brpsel_min <- pep_brpsel_min %>% group_by(id) %>% top_n(1, nuxl_score)

pep_brpsel_min$peak_ann <- str_replace_all(pep_brpsel_min$peak_ann, "\\|", ";")
pep_brpsel_min$peak_ann <- str_replace_all(pep_brpsel_min$peak_ann, ",", "|")


uvpselall<-pep_brpsel_min

#exclude redundancy 
uvpselall <- uvpselall %>% group_by(id) %>% top_n(1, nuxl_score)
uvpselall <- uvpselall %>% group_by(id) %>% top_n(1, RT)

uvpselall <- uvpselall %>% mutate (stripped = sequence)

uvpselall$stripped <-  gsub ('(Oxidation)','',uvpselall$stripped, fixed=T) 
#try loop for binding spectral library from nuxl table

#loop on short df (subset big table) - it fucking works!

#tab541s <- tab541msel[1:5,] #subset for faster work, now use tjsel 
#and its derivatives

uvpselall <- uvpselall[,-c(1,2)]

#dfnew <- df2[1,] #just first line in the starting db
#dfnew<-dfnew[,-6]
#dfnew<-dfnew%>% mutate(PrecursorMz <- NA)

dfnew<-read_tsv('template_lib.tsv')

for(i in 1:nrow(uvpselall)) {       # for-loop over rows
  
  datal <- as.character(uvpselall[i,7]) 
  
  rowsl <- strsplit(datal, ';')[[1]]
  
  
  dfl<-tibble(rowsl) 
  df3<-dfl %>%  separate_wider_delim(rowsl, "|", names = 
                                       c("FragmentMz", "RelativeIntensity",
                                         "FragmentCharge","annotation"))
  
  
  df3 <- df3 %>% mutate(ModifiedPeptide=as.matrix(uvpselall[i,5])
                        )
  
  df3 <- df3 %>% mutate(StrippedPeptide=as.matrix(uvpselall[i,9])
  )
  
  df3 <- df3 %>% mutate(PrecursorCharge=as.matrix(uvpselall[i,2])
                        )
  
  df3 <- df3 %>% mutate(RT=as.matrix(uvpselall[i,3])
                        )
  
  df3 <- df3 %>% mutate(IonMobility=as.matrix(uvpselall[i,4])
                        )
  
  df3 <- df3 %>% mutate(ProteinID=as.matrix(uvpselall[i,6])
                        ) 
  
  df3 <- df3 %>% mutate(PrecursorMz=as.matrix(uvpselall[i,8])
  )
  
  dfnew <- rbind(dfnew,df3)
  
}




#leave b and y ions only
dfnew2<-dfnew %>% filter( grepl('b',annotation) | grepl('y',annotation)  ) 
#create fragment type from annotation
dfnew2$FragmentType <- ifelse(grepl("b", dfnew2$annotation), "b", "y")

#delete first original extra row
dfnew2<-dfnew2[-1,]

colnames(dfnew2)[8] = "FragmentAnnotation"
colnames(dfnew2) [3] = "Tr_recalibrated"

write.table(dfnew2, file = "UVECO_PE_lib_shot130925.tsv", row.names=FALSE, sep="\t")

#run XL library
uv_xlfilt_pep <- rbind (dfnew1, dfnew2)

write.table(uv_xlfilt_pep, file = "UVECO_SHOT_SPLIBRARY_130126.tsv", row.names=FALSE, sep="\t")

