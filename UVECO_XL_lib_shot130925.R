library(tidyverse)
library(tidyr)
library(dplyr); library(readr); library(ggplot2); library(ggpubr)
library (arsenal)
library(stringr)
library(base)
library(data.table)
library(stringr)

#ALL FOR UV 
setwd("D:\\ECOLI_2025\\Ecoli_UVXL_DDA_singleshots_2026\\003-TextExporter-out")

file_list_UV <- list.files(pattern = "*.XLs\\.unknown$", full.names = TRUE)

# Read all files and store them in a list of data frames
df_list_UV <- lapply(file_list_UV, read_tsv)

# Combine all data frames into one
xl_UV <- do.call(rbind, df_list_UV)




colnames(xl_UV)[1]  <- "RT"
colnames(xl_UV)[36]  <- "nucleic"
colnames(xl_UV)[54]  <- "XL"
colnames(xl_UV)[29]  <- "CCS"
colnames(xl_UV)[46]  <- "best_loc"
colnames(xl_UV)[47]  <- "loc"
colnames(xl_UV)[48]  <- "locscore"
colnames (xl_UV) [68] <- "nuxl_score"

colnames (xl_UV) [16] <- "peak_ann"


xl_brpj <- xl_UV %>%  mutate (id = paste(xl_UV$sequence,xl_UV$nucleic,xl_UV$charge,sep="_"))

setwd('D:\\ECOLI_2025\\Ecoli_UVXL_DDA_singleshots_2026\\library_shots')


#back to the library
xl_brpj$peak_ann <-  gsub ('"','',xl_brpj$peak_ann, fixed=T)

xl_brp_sel <- xl_brpj %>% dplyr::select (id,charge,RT,CCS,
                                         sequence,accessions,nucleic,peak_ann,loc,mz,nuxl_score,
                                         locscore)


xl_brp_sel$peak_ann <- str_replace_all(xl_brp_sel$peak_ann, "\\|", ";")
xl_brp_sel$peak_ann <- str_replace_all(xl_brp_sel$peak_ann, ",", "|")


shotbrp_sel <- xl_brp_sel

shotbrp_sel <- shotbrp_sel %>% filter(locscore>0) #4041


tj_max <- shotbrp_sel %>% group_by(id) %>% top_n(1, nuxl_score) 
tj_max <-  distinct(tj_max) #519

tjmax <- tj_max %>% group_by(id) %>% top_n(1, locscore) 

tkmax <- ungroup (tjmax) %>% select (sequence) %>% distinct () #207


tj_maxl <- shotbrp_sel %>% group_by(id) %>% top_n(1, locscore) 
tj_maxl <- tj_maxl %>% group_by(id) %>% top_n(1, nuxl_score) 
tj_maxl <-  distinct(tj_maxl) #509



tj_n <- shotbrp_sel %>% group_by(id) %>% summarise(n=n())
tj_n$n <- as.numeric(tj_n$n)


#do not forget mz
tj_sel <- tj_maxl %>% dplyr::select (id,charge,RT,CCS,sequence,accessions,nucleic,peak_ann,loc,mz)

#CCS must be converted back to IM! 


tj_sel <- left_join(tj_sel, tj_n, by=c("id"="id") ) 


uv_sel <- tj_sel %>% mutate(IM = CCS*2.6867811*(10^25)/(charge*3/16*10000*(10^20)*((2*pi/(((mz*charge)*28.0134/(mz*charge+28.0134)*1.66053904*(10^(-27)))*1.38064852*(10^(-23))*(305)))^0.5)*1.6021766208*(10^(-19))))


#introducing required format of ModifiedSequence

uv_sel <- uv_sel %>% mutate (local = loc+1)


#exclude hits with CSM=1  
#uv_sel <- uv_sel %>% filter (n > 1) #

#exclude hits with charge 4 and 5
#uv_sel <- uv_sel %>% filter (charge < 4) #303 vs 556

uv_sel <- uv_sel[,-11]



#insert parentheses to nucleic to insert it further as modification


uv_sel$nucleic = gsub('.*^','(', uv_sel$nucleic)

uv_sel$par = ')' 

uv_sel$nucleic <- paste(uv_sel$nucleic,uv_sel$par)
uv_sel$nucleic = gsub(' ','', uv_sel$nucleic)
uv_sel <- uv_sel[,-13]
uv_sell <- uv_sel[,-4]
#uv_sell$local = gsub('0','1', uv_sell$local) #modification must be at the residue
#divide uv_sel based on sequence - including and excluding modifications

uv_nomod <- uv_sell %>% filter (
  !grepl("(",sequence, fixed=T )
  
) 

uv_nomod <- uv_nomod %>% mutate (mod=NA)


uv_modo <- uv_sell %>% filter (
  grepl("(O",sequence, fixed=T )
) 
uv_modo$mod = NA



fun_insert <- function(x, pos, insert) {       # Create own function
  gsub(paste0("^(.{", pos, "})(.*)$"),
       paste0("\\1", insert, "\\2"),
       x)
}





for(i in 1:nrow(uv_nomod)) {
  
  uv_nomod[i,12] <- fun_insert(x = uv_nomod[i,4],    # Apply own function
                               pos = uv_nomod[i,11], 
                               insert = uv_nomod[i,6])
}

#it works, then work with tj_modo and tj_modc - existing modifications

#Oxidation M
for(i in 1:nrow(uv_modo)) {
  
  uv_modo[i,12] <- fun_insert(x = uv_modo[i,4],    # Apply own function
                              pos = uv_modo[i,11], 
                              insert = uv_modo[i,6])
}


uv_modo_ok <- uv_modo %>% filter (
  grepl("M(Oxidation)",mod, fixed=T )
)  

uv_modo_nok <- setdiff(uv_modo, uv_modo_ok)


for(i in 1:nrow(uv_modo_nok)) {
  
  uv_modo_nok[i,12] <- fun_insert(x = uv_modo_nok[i,4],    # Apply own function
                              pos = uv_modo_nok[i,11] + 11, 
                              insert = uv_modo_nok[i,6])
  }


#MANUALLY CORRECT 3 LINES 
uv_modo_nok [1,12] <- "ATLGEVGNAEH(U)M(Oxidation)LR" 
uv_modo_nok [4,12] <- "ATLGEVGNAEH(C-H3N1)M(Oxidation)LR" 
uv_modo_nok [5,12] <- "ATLGEVGNAEH(CG-H3N1-H2O1)M(Oxidation)LR" 


uv_nomod$local<-as.numeric(uv_nomod$local)
uv_modo_ok$local<-as.numeric(uv_modo_ok$local)
uv_modo_nok$local<-as.numeric(uv_modo_nok$local)

uvlib <- rbind (uv_nomod, uv_modo_ok,uv_modo_nok)



uvl <- uvlib

uvl$sequence <-  gsub ('(Oxidation)','',uvl$sequence, fixed=T)     
#tjl$sequence <-  gsub ('(Carbamidomethyl)','',tjl$sequence, fixed=T)     


#try loop for binding spectral library from nuxl table

#loop on short df (subset big table) - it fucking works!

#tab541s <- tab541msel[1:5,] #subset for faster work, now use tjsel 
#and its derivatives

uvl <- uvl %>% filter(!grepl("DECOY", accessions))

#dfnew <- df2[1,] #just first line in the starting db
#dfnew<-dfnew[,-6]
#dfnew<-dfnew%>% mutate(PrecursorMz <- NA)

dfnew<-read_tsv('template_lib.tsv')

for(i in 1:nrow(uvl)) {       # for-loop over rows
  
  datal <- as.character(uvl[i,7]) 
  
  rowsl <- strsplit(datal, ';')[[1]]
  
  
  dfl<-tibble(rowsl) 
  df3<-dfl %>%  separate_wider_delim(rowsl, "|", names = 
                                       c("FragmentMz", "RelativeIntensity",
                                         "FragmentCharge","annotation"))
  
  
  df3 <- df3 %>% mutate(ModifiedPeptide=as.matrix(uvl[i,12])
  )
  
  df3 <- df3 %>% mutate(StrippedPeptide=as.matrix(uvl[i,4])
  )
  
  df3 <- df3 %>% mutate(PrecursorCharge=as.matrix(uvl[i,2])
  )
  
  df3 <- df3 %>% mutate(RT=as.matrix(uvl[i,3])
  )
  
  df3 <- df3 %>% mutate(IonMobility=as.matrix(uvl[i,10])
  )
  
  df3 <- df3 %>% mutate(ProteinID=as.matrix(uvl[i,5])
  ) 
  
  df3 <- df3 %>% mutate(PrecursorMz=as.matrix(uvl[i,9])
  )
  
  dfnew <- rbind(dfnew,df3)
  
}



#leave b and y ions only
dfnew1<-dfnew %>% filter( grepl('b',annotation) | grepl('y',annotation)  ) 
#create fragment type from annotation
dfnew1$FragmentType <- ifelse(grepl("b", dfnew1$annotation), "b", "y")

#delete first original extra row
dfnew1<-dfnew1[-1,]

colnames(dfnew1)[8] = "FragmentAnnotation"
colnames(dfnew1) [3] = "Tr_recalibrated"

write.table(dfnew1, file = "UVECO_XL_lib_shot130925.tsv", row.names=FALSE, sep="\t")

#add library for peptides

#peplibdf <- read_tsv("nucuv_pbrplib_150425.tsv") 

#libmay <- rbind(dfnew1,peplibdf)
#write.table(libmay, file = "ne_uvxl_lib060525.tsv", row.names=FALSE, sep="\t")
