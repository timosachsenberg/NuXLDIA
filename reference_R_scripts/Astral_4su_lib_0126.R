library(tidyverse)
library(tidyr)
library(dplyr); library(readr); library(ggplot2); library(ggpubr)
library (arsenal)
library(stringr)
library(base)
library(data.table)
library(stringr)
library(stringi)


################################################################
#CROSSLINK LIBRARY



setwd("D:\\Riboastral_processing\\Rosa_human_2025\\splib2026\\003-TextExporter-out")

file_list <- list.files(pattern = "*.XLs\\.unknown$", full.names = TRUE)

# Read all files and store them in a list of data frames
df_list <- lapply(file_list, read_tsv)

# Combine all data frames into one
xl_s1 <- do.call(rbind, df_list)




colnames(xl_s1)[1]  <- "RT"
colnames(xl_s1)[35]  <- "nucleic"
#colnames(xl_s1)[49]  <- "XL"
#colnames(xl_s1)[28]  <- "CCS"
colnames(xl_s1)[42]  <- "best_loc"
colnames(xl_s1)[43]  <- "loc"
colnames(xl_s1)[44]  <- "locscore"
colnames (xl_s1) [64] <- "nuxl_score"

colnames (xl_s1) [16] <- "peak_ann"


xl_s1j <- xl_s1 %>%  mutate (id = paste(xl_s1$sequence,xl_s1$nucleic,xl_s1$charge,sep="_"))

setwd("D:\\Riboastral_processing\\Rosa_human_2025\\splib2026")

xl_s1j$peak_ann <-  gsub ('"','',xl_s1j$peak_ann, fixed=T)

xl_s1_sel <- xl_s1j %>% dplyr::select (id,charge,RT,sequence,accessions,nucleic,peak_ann,loc,mz,nuxl_score,locscore)

#xl_s1_sel <- xl_s1_sel %>%  mutate(origin="s1")

#Tune peak annotations as for manual production of lists 

xl_s1_sel$peak_ann <- str_replace_all(xl_s1_sel$peak_ann, "\\|", ";")
xl_s1_sel$peak_ann <- str_replace_all(xl_s1_sel$peak_ann, ",", "|")



#s1_sel<- xl_s1_sel #5745


s1_sel <- xl_s1_sel %>% filter(locscore>0) #4124

tj_maxl <- s1_sel %>% group_by(id) %>% top_n(1, locscore) 
tj_maxl <- tj_maxl %>% group_by(id) %>% top_n(1, nuxl_score) 
tj_maxl <-  distinct(tj_maxl) #1237


tj_n <- s1_sel %>% group_by(id) %>% summarise(n=n())
tj_n$n <- as.numeric(tj_n$n)

#do not forget mz
tj_sel <- tj_maxl %>% dplyr::select (id,charge,RT,sequence,accessions,nucleic,peak_ann,loc,mz)

tj_sel <- left_join(tj_sel, tj_n, by=c("id"="id") )  #1237


#introducing required format of ModifiedSequence

uv_sel <- tj_sel %>% mutate (local = loc+1)

#remove n - needed for understanding where we are by numbers 
uv_sel <- uv_sel[,-10]
#insert parentheses to nucleic to insert it further as modification


uv_sel$nucleic = gsub('.*^','(', uv_sel$nucleic)

uv_sel$par = ')' 

uv_sel$nucleic <- paste(uv_sel$nucleic,uv_sel$par)
uv_sel$nucleic = gsub(' ','', uv_sel$nucleic)
uv_sel <- uv_sel[,-11]

#select sequences w/o oxidation M
uv_nomod <- uv_sel %>% filter (
  !grepl("(",sequence, fixed=T )
  ) 

uv_nomod <- uv_nomod %>% mutate (mod=NA)

fun_insert <- function(x, pos, insert) {       # Create own function
  gsub(paste0("^(.{", pos, "})(.*)$"),
       paste0("\\1", insert, "\\2"),
       x)
  }





for(i in 1:nrow(uv_nomod)) {
  
  uv_nomod[i,11] <- fun_insert(x = uv_nomod[i,4],    # Apply own function
                               pos = uv_nomod[i,10], 
                               insert = uv_nomod[i,6])
  }



#oxidations M present

uv_modo <- uv_sel %>% filter (
  grepl("(O",sequence, fixed=T )
) 

uv_modo$ox <- str_count(uv_modo$sequence, fixed("Oxidation"))


uv_modo1 <- uv_modo %>% filter (ox == 1)

uv_modo1$seq <- uv_modo1$sequence 

uv_modo1$sequence <- gsub ("\\(Oxidation\\)","", uv_modo1$sequence)

#put positions of O to the $oxpos

positions <- str_locate_all(uv_modo1$seq, fixed("O"))

# Store as comma-separated string of positions
uv_modo1$oxpos <- sapply(positions, function(x) {
  if(nrow(x) == 0) {
    NA  # or "" for empty string
  } else {
    paste(x[, "start"], collapse = ",")
  }
})


uv_modo1$oxpos <- as.numeric(uv_modo1$oxpos)

uv_modo1$oxpos <- uv_modo1$oxpos - 2 

uv_modo1$mod = NA






for(i in 1:nrow(uv_modo1)) {
  
  uv_modo1[i,14] <- fun_insert(x = uv_modo1[i,4],    # Apply own function
                               pos = uv_modo1[i,10], 
                               insert = uv_modo1[i,6])
  }


uv_modo1$nucount <- nchar(uv_modo1$mod) - nchar(uv_modo1$sequence)

#stringi

uv_modo1$mod_new <- mapply(function(mod, oxpos, local, nucount) {
  if(is.na(oxpos) || oxpos == "" || is.na(local)) return(mod)
  
  positions <- as.numeric(strsplit(as.character(oxpos), ",")[[1]])
  positions <- sort(positions, decreasing = TRUE)
  
  for(pos in positions) {
    # Calculate position to insert AFTER
    if(pos <= local) {
      insert_after <- pos  # We'll insert after this position
    } else {
      insert_after <- pos + nucount
    }
    
    # Insert AFTER the character at insert_after position
    mod <- stri_sub_replace(mod,
                            from = insert_after + 1,
                            to = insert_after,
                            replacement = "(Oxidation)")
  }
  return(mod)
}, uv_modo1$mod, uv_modo1$oxpos, uv_modo1$local, uv_modo1$nucount, SIMPLIFY = TRUE)


uv_modo1_ready <- uv_modo1 %>% dplyr::select (-ox,-mod,-seq,-oxpos,-nucount)

colnames(uv_modo1_ready) [11] = "mod"

# 2 Oxidations M

uv_modo2 <- uv_modo %>% filter (ox == 2)

uv_modo2$seq <- uv_modo2$sequence 

uv_modo2$sequence <- gsub ("\\(Oxidation\\)","", uv_modo2$sequence)

#put positions of O to the $oxpos

positions2 <- str_locate_all(uv_modo2$seq, fixed("O"))

# Store as comma-separated string of positions
uv_modo2$oxpos <- sapply(positions2, function(x) {
  if(nrow(x) == 0) {
    NA  # or "" for empty string
  } else {
    paste(x[, "start"], collapse = ",")
  }
})

uv_modo2 <- uv_modo2 %>%
  separate(oxpos, into = c("oxpos1", "oxpos2"), sep = ",", remove = TRUE)


uv_modo2$oxpos1 <- as.numeric(uv_modo2$oxpos1)

uv_modo2$oxpos2 <- as.numeric(uv_modo2$oxpos2)

uv_modo2$oxpos1 <- uv_modo2$oxpos1 - 2 

uv_modo2$oxpos2 <- uv_modo2$oxpos2 - 2 - 11 

uv_modo2$mod = NA


for(i in 1:nrow(uv_modo2)) {
  
  uv_modo2[i,15] <- fun_insert(x = uv_modo2[i,4],    # Apply own function
                               pos = uv_modo2[i,10], 
                               insert = uv_modo2[i,6])
  }


uv_modo2$nucount <- nchar(uv_modo2$mod) - nchar(uv_modo2$sequence)

#stringi
uv_modo2$mod_final <- mapply(function(mod, oxpos1, oxpos2, local, nucount) {
  if((is.na(oxpos1) || oxpos1 == "") && (is.na(oxpos2) || oxpos2 == "")) return(mod)
  
  # Combine all oxidation positions
  positions <- c()
  if(!is.na(oxpos1) && oxpos1 != "") {
    positions <- c(positions, as.numeric(strsplit(as.character(oxpos1), ",")[[1]]))
  }
  if(!is.na(oxpos2) && oxpos2 != "") {
    positions <- c(positions, as.numeric(strsplit(as.character(oxpos2), ",")[[1]]))
  }
  
  # Sort in descending order for correct insertion
  positions <- sort(positions, decreasing = TRUE)
  
  # Process all positions
  for(pos in positions) {
    # Calculate position to insert AFTER
    if(pos <= local) {
      insert_after <- pos
    } else {
      insert_after <- pos + nucount
    }
    
    # Insert "(Oxidation)" AFTER the character at insert_after position
    mod <- stri_sub_replace(mod,
                            from = insert_after + 1,
                            to = insert_after,
                            replacement = "(Oxidation)")
  }
  return(mod)
}, uv_modo2$mod, uv_modo2$oxpos1, uv_modo2$oxpos2, uv_modo2$local, uv_modo2$nucount, SIMPLIFY = TRUE)



uv_modo2_ready <- uv_modo2 %>% dplyr::select (-ox,-mod,-seq,-oxpos1, -oxpos2, -nucount)

colnames(uv_modo2_ready) [11] = "mod"



uvl <- rbind (uv_nomod, uv_modo1_ready, uv_modo2_ready)


uvl$sequence <-  gsub ('(Oxidation)','',uvl$sequence, fixed=T)     

dfnew<-read_tsv('template_lib.tsv')

dfnew<-dfnew[,-4]

for(i in 1:nrow(uvl)) {       # for-loop over rows
  
  datal <- as.character(uvl[i,7]) 
  
  rowsl <- strsplit(datal, ';')[[1]]
  
  
  dfl<-tibble(rowsl) 
  df3<-dfl %>%  separate_wider_delim(rowsl, "|", names = 
                                       c("FragmentMz", "RelativeIntensity",
                                         "FragmentCharge","annotation"))
  
  
  df3 <- df3 %>% mutate(ModifiedPeptide=as.matrix(uvl[i,11])
  )
  
  df3 <- df3 %>% mutate(StrippedPeptide=as.matrix(uvl[i,4])
  )
  
  df3 <- df3 %>% mutate(PrecursorCharge=as.matrix(uvl[i,2])
  )
  
  df3 <- df3 %>% mutate(RT=as.matrix(uvl[i,3])
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

colnames(dfnew1)[7] = "FragmentAnnotation"
colnames(dfnew1) [3] = "Tr_recalibrated"

write.table(dfnew1, file = "pacce_xl_astr220126.tsv", row.names=FALSE, sep="\t")

######################################################################
#PEPTIDE LIBRARY

setwd("D:\\Riboastral_processing\\Rosa_human_2025\\splib2026\\003-TextExporter-out")

file_list <- list.files(pattern = "*.peptides\\.unknown$", full.names = TRUE)

# Read all files and store them in a list of data frames
df_list <- lapply(file_list, read_tsv)

# Combine all data frames into one
uva <- do.call(rbind, df_list)


colnames(uva)[1]  <- "RT"
colnames(uva)[35]  <- "nucleic"
#colnames(uva)[49]  <- "XL"
#colnames(uva)[28]  <- "CCS"
colnames(uva)[42]  <- "best_loc"
colnames(uva)[43]  <- "loc"
colnames(uva)[44]  <- "locscore"
colnames (uva) [64] <- "nuxl_score"

colnames (uva) [16] <- "peak_ann"


uvaj <- uva %>%  mutate (id = paste(uva$sequence,uva$charge,sep="_"))

setwd("D:\\Riboastral_processing\\Rosa_human_2025\\splib2026")

uvaj$peak_ann <-  gsub ('"','',uvaj$peak_ann, fixed=T)


uvaj$peak_ann <- str_replace_all(uvaj$peak_ann, "\\|", ";")
uvaj$peak_ann <- str_replace_all(uvaj$peak_ann, ",", "|")

#select maximal NuXlscore at the same sequence, prec_charge, nucleic
tj_max <- uvaj %>% group_by(id) %>% top_n(1, nuxl_score) 
tj_max <-  distinct(tj_max)
tj_max <- tj_max %>% group_by(id) %>% top_n(-1, score) #6913

tj_n <- uvaj %>% group_by(id) %>% summarise(n=n())
tj_n$n <- as.numeric(tj_n$n)




#do not forget mz
tj_sel <- tj_max %>% dplyr::select (id,charge,RT,sequence,accessions,nucleic,peak_ann,loc,mz)


tj_sel <- left_join(tj_sel, tj_n, by=c("id"="id") ) 



#exclude hits with CSM=1 
#uv_sel <- tj_sel %>% filter (n > 1) #5615

#exclude hits with charge 5
#uv_sel <- uv_sel %>% filter (charge < 5) #795

#remove n
uv_sel <- uv_sel[,-10]
#insert parentheses to nucleic to insert it 




#insert parentheses to nucleic to insert it further as modification


#uv_sel$nucleic = gsub('.*^','(', uv_sel$nucleic)

#uv_sel$par = ')' 

#uv_sel$nucleic <- paste(uv_sel$nucleic,uv_sel$par)
#uv_sel$nucleic = gsub(' ','', uv_sel$nucleic)
#uv_sel <- uv_sel[,-11]
#uv_sell <- uv_sel[,-4]
#uv_sell <- uv_sel
#uv_sell$local = gsub('0','1', uv_sell$local) #modification must be at the residue
#divide uv_sel based on sequence - including and excluding modifications

#select sequences w/o oxidation M
#uv_nomod <- uv_sell %>% filter (
#  !grepl("(",sequence, fixed=T ))


#uv_nomod$mod = NA

#tj_modoc <- tjsel %>% filter (
# grepl("(O",sequence, fixed=T) &
#  grepl("(C",sequence, fixed=T)
#)
#tj_modoc$mod = NA
#select lines with oxidation M
#uv_modo <- uv_sell %>% filter (
#  grepl("(O",sequence, fixed=T )
#) 
#uv_modo$mod = NA
#tj_modo <- setdiff(tj_modo,tj_modoc)


#tj_modc <- tjsel %>% filter (
# grepl("(C",sequence, fixed=T )
#) 
#tj_modc$mod = NA
#tj_modc <- setdiff(tj_modc,tj_modoc)


#function inserting modification to the library

#fun_insert <- function(x, pos, insert) {       # Create own function
#  gsub(paste0("^(.{", pos, "})(.*)$"),
#       paste0("\\1", insert, "\\2"),
#       x)




#my_string_new1 <- fun_insert(x = tj_nomod$sequence,    # Apply own function
#                            pos = tj_nomod$local, 
#                           insert = tj_nomod$nucleic)
#my1<-as.data.frame(my_string_new1)  
#tj_nomod$mod <- my1 #в таком виде запоминается только local из первой строчки (2)
#и дальше лупит по ним везде как будто local=2. нужен цикл?      


#for(i in 1:nrow(uv_nomod)) {

#  uv_nomod[i,11] <- fun_insert(x = uv_nomod[i,4],    # Apply own function
#                               pos = uv_nomod[i,10], 
#                               insert = uv_nomod[i,6])
#}

#it works, then work with tj_modo and tj_modc - existing modifications

#Oxidation M
#for(i in 1:nrow(uv_modo)) {

#  uv_modo[i,11] <- fun_insert(x = uv_modo[i,4],    # Apply own function
#                              pos = uv_modo[i,10], 
#                              insert = uv_modo[i,6])
#}



uvl <- uv_sel

uvl<- uvl %>% mutate (Modified.Sequence = sequence)


uvl$sequence <-  gsub ('(Oxidation)','',uvl$sequence, fixed=T)     
#tjl$sequence <-  gsub ('(Carbamidomethyl)','',tjl$sequence, fixed=T)     


#try loop for binding spectral library from nuxl table

#loop on short df (subset big table) - it fucking works!

#tab541s <- tab541msel[1:5,] #subset for faster work, now use tjsel 
#and its derivatives



#dfnew <- df2[1,] #just first line in the starting db
#dfnew<-dfnew[,-6]
#dfnew<-dfnew%>% mutate(PrecursorMz <- NA)

dfnew<-read_tsv('template_lib.tsv')

dfnew<-dfnew[,-4]

for(i in 1:nrow(uvl)) {       # for-loop over rows
  
  datal <- as.character(uvl[i,7]) 
  
  rowsl <- strsplit(datal, ';')[[1]]
  
  
  dfl<-tibble(rowsl) 
  df3<-dfl %>%  separate_wider_delim(rowsl, "|", names = 
                                       c("FragmentMz", "RelativeIntensity",
                                         "FragmentCharge","annotation"))
  
  
  df3 <- df3 %>% mutate(ModifiedPeptide=as.matrix(uvl[i,10])
  )
  
  df3 <- df3 %>% mutate(StrippedPeptide=as.matrix(uvl[i,4])
  )
  
  df3 <- df3 %>% mutate(PrecursorCharge=as.matrix(uvl[i,2])
  )
  
  df3 <- df3 %>% mutate(RT=as.matrix(uvl[i,3])
  )
  
  df3 <- df3 %>% mutate(ProteinID=as.matrix(uvl[i,5])
  ) 
  
  df3 <- df3 %>% mutate(PrecursorMz=as.matrix(uvl[i,9])
  )
  
  dfnew <- rbind(dfnew,df3)
  
}




#leave b and y ions only
dfnew2<-dfnew %>% filter( grepl('b',annotation) | grepl('y',annotation)  ) 
#create fragment type from annotation
dfnew2$FragmentType <- ifelse(grepl("b", dfnew2$annotation), "b", "y")

#delete first original extra row
dfnew2<-dfnew2[-1,]

colnames(dfnew2)[7] = "FragmentAnnotation"
colnames(dfnew2) [3] = "Tr_recalibrated"

write.table(dfnew2, file = "pacce_peplib220126.tsv", row.names=FALSE, sep="\t")

#run XL library
paccelib <- rbind (dfnew1, dfnew2)

write.table(paccelib, file = "pacce_SPECTRAL_LIB_220125.tsv", row.names=FALSE, sep="\t")






#############################
#for compatibility with dia-nn 2+


paccelib_irt <- paccelib %>% mutate (predicted_iRT = Tr_recalibrated) 

#rename columns, add fragment numbers



#add fragment numbers
paccelib_irt <- paccelib_irt  %>%
  extract(FragmentAnnotation, 
          into = "FragmentSeriesNumber", 
          regex = ".*?(\\d+).*", 
          remove = FALSE) %>%
  mutate(FragmentSeriesNumber = as.numeric(FragmentSeriesNumber))

paccelib_irt$FragmenLossType <- ifelse(grepl("'", paccelib_irt$FragmentAnnotation),
                                       "unknown", NA)

colnames (paccelib_irt) [1] = "ModifiedPeptideSequence" 

# Replace (C...), (U...), (A...), (G...), (S...) with [...]
paccelib_irt$ModifiedPeptideSequence <- gsub(
  "\\((C[^)]*)\\)", 
  "[\\1]", 
  paccelib_irt$ModifiedPeptideSequence
)
paccelib_irt$ModifiedPeptideSequence <- gsub(
  "\\((U[^)]*)\\)", 
  "[\\1]", 
  paccelib_irt$ModifiedPeptideSequence
)
paccelib_irt$ModifiedPeptideSequence <- gsub(
  "\\((A[^)]*)\\)", 
  "[\\1]", 
  paccelib_irt$ModifiedPeptideSequence
)
paccelib_irt$ModifiedPeptideSequence <- gsub(
  "\\((G[^)]*)\\)", 
  "[\\1]", 
  paccelib_irt$ModifiedPeptideSequence
)

paccelib_irt$ModifiedPeptideSequence <- gsub(
  "\\((S[^)]*)\\)", 
  "[\\1]", 
  paccelib_irt$ModifiedPeptideSequence
)

paccelib_irt$ModifiedPeptideSequence <- gsub(
  "\\(Oxidation\\)", 
  "\\(UniMod:35\\)", 
  paccelib_irt$ModifiedPeptideSequence
)

colnames (paccelib_irt) [3] = "AverageExperimentalRetentionTime"

#colnames (paccelib_irt) [4] = "PrecursorIonMobility"

colnames (paccelib_irt) [4] = "PeptideSequence"

colnames (paccelib_irt) [6] = "ProteinId"

colnames (paccelib_irt) [7] = "Annotation"

colnames (paccelib_irt) [9] = "ProductMz"

colnames (paccelib_irt) [10] = "LibraryIntensity"

colnames (paccelib_irt) [13] = "NormalizedRetentionTime"



write.table(paccelib_irt, file = "PACCE_LIB_DIANN2_220126.tsv", row.names=FALSE, 
            sep="\t", na = "")






