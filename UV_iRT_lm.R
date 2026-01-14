library(tidyverse)
library(readr)
library(ggplot2)
library(ggpubr)
library(stringr)
library(base)

setwd("D:\\ECOLI_2025\\Ecoli_UVXL_DDA_singleshots_2026\\library_shots")

helapsm <- read.table("D:\\ECOLI_2025\\Ecoli_UVXL_DDA_singleshots_2026\\library_shots\\HeLa_fragout\\psm.tsv",
                     header=T, sep = "\t")

helalib <- read.table ("D:\\ECOLI_2025\\Ecoli_UVXL_DDA_singleshots_2026\\library_shots\\HeLa_fragout\\library.tsv", 
                       header=T, sep = "\t")

helapsm$modpep <- ifelse(
  !is.na(helapsm$Modified.Peptide) & helapsm$Modified.Peptide != "", 
  helapsm$Modified.Peptide, 
  helapsm$Peptide
)

helapsm <- helapsm %>%  mutate (
  id = paste(helapsm$modpep, helapsm$Charge,sep="_")
                                )

helalib <- helalib %>%  mutate (
  id = paste(helalib$ModifiedPeptideSequence, helalib$PrecursorCharge,sep="_")
)


helapsms <- helapsm %>% dplyr::select (id, Retention) 




helapsms <- helapsms %>% group_by (id) %>% summarize (
  RT = mean(Retention)
                                                  )

#remove modified peptides

helapsms <- helapsms %>% filter (
  !grepl ("\\[",id)
                                 )


helalibs1 <- helalib %>% dplyr::select (id, NormalizedRetentionTime) %>% distinct()


helalibs1 <- helalibs1 %>% filter (
  !grepl ("UniMod",id)
)


helart <- helalib %>% dplyr::select (id, NormalizedRetentionTime, 
                                       AverageExperimentalRetentionTime) %>% distinct()



helart1 <- left_join (helapsms, helalibs1, by = "id") 

helart1 <- helart1 %>% drop_na()


ggplot(helart1, aes(x = RT, y = NormalizedRetentionTime)) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue") +
  labs(title = "",
       x = "RT",
       y = "iRT") +
  theme_minimal()

colnames (helart1) [3] = "iRT"


colnames (helart) [2] = "iRT"
colnames (helart) [3] = "RT"

ggplot(helart, aes(x = RT, y = iRT)) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue") +
  labs(title = "",
       x = "RT",
       y = "iRT") +
  theme_minimal()

ggplot(helart, aes(x = RT, y = iRT)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", color = "blue", se = TRUE) +
  labs(title = "",
       x = "RT",
       y = "iRT") +
  theme_minimal() +
  scale_x_continuous(breaks = seq(round(min(helart$RT)/50)*50, 
                                  round(max(helart$RT)/50)*50, 
                                  by = 50)) +
  theme(panel.grid.major = element_line(color = "gray85", linewidth = 0.5),
        panel.grid.minor = element_blank())  # Remove minor grid lines


# Install and load required package

library(segmented)

# Fit initial linear model
fit_lm <- lm(iRT ~ RT, data = helart)

# Fit segmented regression with 2 breakpoints (3 segments)
fit_seg <- segmented(fit_lm, seg.Z = ~RT, npsi = 3)

# Summary of the model
summary(fit_seg)

# Get breakpoints
breakpoints <- fit_seg$psi[, 2]
print(paste("Breakpoints at x =", breakpoints[1], "and x =", breakpoints[2],
            "and x =", breakpoints[3]))

# Plot the results
plot(helart$RT, helart$iRT, pch = 16, col = "gray",
     xlab = "RT", ylab = "iRT", main = "Piecewise Linear Regression")
plot(fit_seg, add = TRUE, col = "red", lwd = 2)
points(fit_seg, col = "blue", pch = 19)

# Function to predict new y values
predict_piecewise <- function(new_x, segmented_model) {
  # Use the predict method for segmented objects
  # The predictor variable in your model is 'RT', not 'x'
  predictions <- predict(segmented_model, newdata = data.frame(RT = new_x))
  return(predictions)
}

# Example: Predict for new x values
new_x <- c(810, 1400, 2106, 2200, 2716)
new_y <- predict_piecewise(new_x, fit_seg)

# Print the predictions
print("Predictions:")
print(data.frame(RT = new_x, predicted_iRT = new_y))

# Get the piecewise linear formula
get_piecewise_formula <- function(segmented_model) {
  coefs <- coef(segmented_model)
  breakpoints <- segmented_model$psi[, 2]
  
  # Intercept and slopes for each segment
  intercept <- coefs[1]
  slope1 <- coefs[2]  # First segment slope
  slope2 <- coefs[2] + coefs[3]  # Second segment slope
  slope3 <- coefs[2] + coefs[3] + coefs[4]  # Third segment slope
  
  formula_text <- paste0(
    "y = \n",
    "  if x < ", breakpoints[1], ": ", round(intercept, 3), " + ", round(slope1, 3), " * x\n",
    "  if ", breakpoints[1], " ≤ x < ", breakpoints[2], ": ", 
    round(intercept - slope1*breakpoints[1] + slope2*breakpoints[1], 3), 
    " + ", round(slope2, 3), " * x\n",
    "  if x ≥ ", breakpoints[2], ": ",
    round(intercept - slope1*breakpoints[1] + slope2*breakpoints[1] - 
            slope2*breakpoints[2] + slope3*breakpoints[2], 3),
    " + ", round(slope3, 3), " * x"
  )
  
  return(formula_text)
}

# Display the formula
cat(get_piecewise_formula(fit_seg))


#let's reshape the UV library 

ecouvlib_rt <- read.table ("D:\\ECOLI_2025\\Ecoli_UVXL_DDA_singleshots_2026\\library_shots\\UVECO_SHOT_SPLIBRARY_130126.tsv", 
                       header=T, sep = "\t")

#predict with segmented lib
lib_x <- ecouvlib_rt$Tr_recalibrated
lib_y <- predict_piecewise(lib_x, fit_seg)

iRT_lib <- data.frame(RT = lib_x, predicted_iRT = lib_y) %>% distinct()

ecouvlib_irt <- left_join (
  ecouvlib_rt, iRT_lib, by = c("Tr_recalibrated" = "RT")
  ) 

#rename columns, add fragment numbers

#add fragment numbers
ecouvlib_irt <- ecouvlib_irt  %>%
    extract(FragmentAnnotation, 
          into = "FragmentSeriesNumber", 
          regex = ".*?(\\d+).*", 
          remove = FALSE) %>%
  mutate(FragmentSeriesNumber = as.numeric(FragmentSeriesNumber))

ecouvlib_irt$FragmenLossType <- ifelse(grepl("'", ecouvlib_irt$FragmentAnnotation),
                                       "unknown", NA)

colnames (ecouvlib_irt) [1] = "ModifiedPeptideSequence" 

# Replace (C...), (U...), (A...), (G...) with [...]
ecouvlib_irt$ModifiedPeptideSequence <- gsub(
  "\\((C[^)]*)\\)", 
  "[\\1]", 
  ecouvlib_irt$ModifiedPeptideSequence
)
ecouvlib_irt$ModifiedPeptideSequence <- gsub(
  "\\((U[^)]*)\\)", 
  "[\\1]", 
  ecouvlib_irt$ModifiedPeptideSequence
)
ecouvlib_irt$ModifiedPeptideSequence <- gsub(
  "\\((A[^)]*)\\)", 
  "[\\1]", 
  ecouvlib_irt$ModifiedPeptideSequence
)
ecouvlib_irt$ModifiedPeptideSequence <- gsub(
  "\\((G[^)]*)\\)", 
  "[\\1]", 
  ecouvlib_irt$ModifiedPeptideSequence
)

ecouvlib_irt$ModifiedPeptideSequence <- gsub(
  "\\(Oxidation\\)", 
  "\\(UniMod:35\\)", 
  ecouvlib_irt$ModifiedPeptideSequence
)

colnames (ecouvlib_irt) [3] = "AverageExperimentalRetentionTime"

colnames (ecouvlib_irt) [4] = "PrecursorIonMobility"

colnames (ecouvlib_irt) [5] = "PeptideSequence"

colnames (ecouvlib_irt) [7] = "ProteinId"

colnames (ecouvlib_irt) [8] = "Annotation"

colnames (ecouvlib_irt) [10] = "ProductMz"

colnames (ecouvlib_irt) [11] = "LibraryIntensity"

colnames (ecouvlib_irt) [14] = "NormalizedRetentionTime"



write.table(ecouvlib_irt, file = "UVECO_SHOT_LIB_iRT_130126.tsv", row.names=FALSE, sep="\t")

